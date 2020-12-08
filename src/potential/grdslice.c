/*--------------------------------------------------------------------
 *	Copyright (c) 1991-2020 by the GMT Team (https://www.generic-mapping-tools.org/team.html)
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation; version 3 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	Contact info: www.generic-mapping-tools.org
 *--------------------------------------------------------------------*/
/* grdslice reads a grid and detects isolated peaks by using
 * slice contouring.  We only contour interior, closed contours.
 * We contour from maximum value and go downwards.  As the contour value
 * drops we will chop off the tops of smaller and smaller peaks.  The
 * first time we do this we initialize a new peak object and then the
 * contour rings further down are appended via pointers to the down slice.
 * As the contours widen it is possible that more than one peak share
 * the same, broader base slices.  We may then traverse the list of peak
 * and their slices to do analysis.  Each slice stores its area and mean
 * location.
 * The usual input grids will in Mercator projection. However, for testing
 * with synthetic data, input grids can be in either km or geographic coordinates.
 * THe original application of this tool was to detect seamounts from grids
 * of gravity or vertical gradient anomalies [see ref list in docs].
 *
 * Author:	Paul Wessel and Seung-Sep Kim
 * Date:	1-DEC-2020 (original date 5-DEC-2006 and partly based on grdcontour)
 * Version:	2.0	Revised for GMT 6
 */

#include "gmt_dev.h"

#define THIS_MODULE_CLASSIC_NAME	"grdslice"
#define THIS_MODULE_MODERN_NAME	"grdslice"
#define THIS_MODULE_LIB		"potential"
#define THIS_MODULE_PURPOSE	"Detect isolated peaks by contour-slicing a grid"
#define THIS_MODULE_KEYS	"<G{"
#define THIS_MODULE_NEEDS	"g"
#define THIS_MODULE_OPTIONS "-:RVf"

struct GRDSLICE_CTRL {
	struct GMT_CONTOUR contour;
	struct GRDSLICE_In {
		bool active;
		char *file;
	} In;
	struct GRDSLICE_A {	/* -A<cutoff> */
		bool active;
		double cutoff;
	} A;
	struct GRDSLICE_C {	/* -C<cont_int> */
		bool active;
		double interval;
	} C;
	struct GRDSLICE_D {	/* -D<dumpfile> */
		bool active;
		char *file;
	} D;
	struct GRDSLICE_L {	/* -L<Low/high> */
		bool active;
		double low, high;
	} L;
	struct GRDSLICE_S {	/* -S<smooth> */
		bool active;
		unsigned int value;
	} S;
	struct GRDSLICE_T {  /* -T<bottom_level>/<area_cutoff> */
		bool active;
		double blevel, acutoff;
	} T;
	struct GRDSLICE_Z {	/* -Z[+s<fact>][+o<shift>] */
		bool active;
		double scale, offset;
	} Z;
};

struct GRDSLICE_SLICE {	/* Hold each contour slice information */
	int n;				/* Number of points in this slice polygon */
	int id;				/* Unique ID number */
	int shared;			/* Number of peaks that share this slice */
	double z;			/* Z-value (contour) of this slice */
	double *x, *y;			/* The array of Mercator (x,y) coordinates */
	double xmin, xmax, ymin, ymax;	/* Extreme Mercator coordinates of polygon */
	double x_mean, y_mean;		/* Mean coordinate as approximation of center point */
	double area;			/* Area of polygon in km^2 */
	double azimuth;			/* Azimuth of major axis of approximate ellipse with same area */
	double major, minor;		/* Length of axes (in km) of approximate ellipse */
	double fit;			/* 0-100% of how well an ellipse explains the shape of contour (100% is perfect) */
	struct GRDSLICE_SLICE *next;		/* Pointer to next slice in same contour level */
	struct GRDSLICE_SLICE *down;		/* Pointer to the next slice down in this stack of slices */
};

struct GRDSLICE_PEAK {	/* Hold start of peak and linked list of slices */
	int id;				/* Unique ID number */
	double x, y;			/* Mean Mercator coordinate of peak */
	double z;			/* Z-value (contour) of this peak */
	struct GRDSLICE_SLICE *start;		/* Pointer to top slice in stack */
};

static void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct GRDSLICE_CTRL *C;
	
	C = gmt_M_memory (GMT, NULL, 1, struct GRDSLICE_CTRL);
	
	/* Initialize values whose defaults are not 0/false/NULL */
	
	C->D.file = strdup ("PEAK");
	C->L.low = 0;
	C->L.high = DBL_MAX;
	C->Z.scale = 1.0;
	
	return (C);
}

void Free_Ctrl (struct GMT_CTRL *GMT, struct GRDSLICE_CTRL *C) {	/* Deallocate control structure */
	gmt_M_free (GMT, C);
}

static double grdslice_contour_area_merc (double x[], double y[], uint64_t n) {
	/* Computes contour area and reverses polygon if needed so always CCW */
	uint64_t i, j;
	double area, xold, yold;
	
	/* Sign will be +ve if polygon is CW, negative if CCW */
	
	area = yold = 0.0;
	xold = x[n-1];	/* We know this is a closed contour already */
	yold = y[n-1];
	
	for (i = 0; i < n; i++) {
		area += (xold - x[i]) * (yold + y[i]);
		xold = x[i];
		yold = y[i];
	}
	if (area > 0.0) {	/* Must reverse */
		for (i = 0, j = n-1; i < n/2; i++, j--) {
			gmt_M_double_swap (x[i], x[j]);
			gmt_M_double_swap (y[i], y[j]);
		}
	}
	return (0.5 * fabs (area));
}

static int usage (struct GMTAPI_CTRL *API, int level) {
	const char *name = gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);

	GMT_Message (API, GMT_TIME_NONE, "usage: %s <grid> -C<cont_int> [-A<cutoff>] [-D<prefix>] [-L<Low/high>]\n", name);
	GMT_Message (API, GMT_TIME_NONE, "\t[%s] [-S<smooth>] [-T<bottom_level>/<area_cutoff>] [%s] [-Z[+s<scale>][+o<offset>]]\n\n", GMT_Rgeo_OPT, GMT_V_OPT);

	if (level == GMT_SYNOPSIS) return (GMT_MODULE_SYNOPSIS);

	GMT_Message (API, GMT_TIME_NONE, "\t<grid> is the grid file to be sliced.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-C Contours slice interval.\n");
	GMT_Message (API, GMT_TIME_NONE, "\n\tOPTIONS:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-A Skip peaks that are within <cutoff> km from a larger peak [no skipping].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-D Sets filename prefix for output files (<prefix>_slices.txt and <prefix>_pos.txt) [PEAK].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-L Limit contours to this range [Default is -L0/inf].\n");
	GMT_Option (API, "R");
	GMT_Message (API, GMT_TIME_NONE, "\t-S Smooth contours by interpolation at approximately gridsize/<smooth> intervals.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-T Specify the bottom level of contouring <bottom_level> and ignore\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   contours whose area are less than <area_cutoff> in km^2\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   With -D, it will write <prefix>_bottom.txt and <prefix>_indices.txt.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   With -L, the bottom must be equal or larger than the low value.\n");
	GMT_Option (API, "V");
	GMT_Message (API, GMT_TIME_NONE, "\t-Z Subtract <shift> (via +o<shift> [0]) and multiply data by <fact> (via +s<fact> [1]).\n");
	GMT_Option (API, "f,.");

	return (GMT_MODULE_USAGE);
}

static int parse (struct GMT_CTRL *GMT, struct GRDSLICE_CTRL *Ctrl, struct GMT_OPTION *options) {
	unsigned int n_errors = 0, n_files = 0, pos = 0;
	int j;
	char p[GMT_LEN64] = {""};
	struct GMT_OPTION *opt = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;

	for (opt = options; opt; opt = opt->next) {	/* Process all the options given */

		switch (opt->option) {

			case '<':	/* Input file (only one is accepted) */
				if (n_files++ > 0) { n_errors++; continue; }
				Ctrl->In.active = true;
				if (opt->arg[0]) Ctrl->In.file = strdup (opt->arg);
				if (GMT_Get_FilePath (GMT->parent, GMT_IS_GRID, GMT_IN, GMT_FILE_REMOTE, &(Ctrl->In.file))) n_errors++;
				break;

			/* Processes program-specific parameters */

			case 'A':
				Ctrl->A.active = true;
				Ctrl->A.cutoff = atof (opt->arg);
				break;
			case 'C':
				Ctrl->C.active = true;
				Ctrl->C.interval = atof (opt->arg);
				break;
			case 'D':
				Ctrl->D.active = true;
				gmt_M_str_free (Ctrl->D.file);
				Ctrl->D.file = strdup (opt->arg);
				break;
			case 'L':
				Ctrl->L.active = true;
				sscanf (opt->arg, "%lf/%lf", &Ctrl->L.low, &Ctrl->L.high);
				break;
			case 'S':
				Ctrl->S.active = true;
				j = atoi (opt->arg);
				n_errors += gmt_M_check_condition (GMT, j < 0, "Option -S: Smooth_factor must be > 0\n");
				Ctrl->S.value = j;
				break;
			case 'T':
				Ctrl->T.active = true;
				sscanf (opt->arg, "%lf/%lf", &Ctrl->T.blevel, &Ctrl->T.acutoff);
				break;
			case 'Z':
				Ctrl->Z.active = true;
				while (gmt_getmodopt (GMT, 'Z', opt->arg, "so", &pos, p, &n_errors) && n_errors == 0) {
					switch (p[0]) {
						case 's':	Ctrl->Z.scale  = atof (&p[1]);	break;
						case 'o':	Ctrl->Z.offset = atof (&p[1]);	break;
						default: 	/* These are caught in gmt_getmodopt so break is just for Coverity */
							break;
					}
				}
				break;
			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, n_files != 1 || Ctrl->In.file == NULL, "Must specify a single grid file\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->C.interval <= 0.0, "Option -C: Must specify positive contour interval\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->L.low > Ctrl->L.high, "Option -L: lower limit > upper!\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->T.active && Ctrl->T.blevel < Ctrl->L.low, "Option -T: Bottom level < lower limit!\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->Z.scale == 0.0, "Option -Z: factor must be nonzero\n");

	return (n_errors ? GMT_PARSE_ERROR : GMT_NOERROR);

}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

EXTERN_MSC int GMT_grdslice (void *V_API, int mode, void *args) {
	bool begin = true, inside, *skip = NULL, is_mercator = false;

	int error, n_skipped = 0, kp_id = 0, cs;

	unsigned int c, n_contours, n_edges, nrots, M = 2, n_inside, n_peaks = 0, n_slices = 0, *edge = NULL;

	uint64_t ij, n, i, n_alloc = GMT_CHUNK;
	int64_t ns;

	char file[PATH_MAX] = {""}, bfile[PATH_MAX] = {""}, kfile[PATH_MAX] = {""};
	
	FILE *fp = NULL, *tp = NULL, *kp = NULL;

	double aspect, cval, min, max, small, scale = 1.0, minor, major, sa, ca, sin_a, cos_a, area, dr, a;
	double small_x, small_y, lon, lat = 0.0, min_area, max_area, xr, yr, r, r_fit, rms, merc_x0 = 0.0, merc_y0 = 0.0;
	double wesn[4], A[4], EigenValue[2], EigenVector[4], out[9], work1[2], work2[2], *x = NULL, *y = NULL, *contour = NULL;

	struct GMT_GRID *G = NULL, *G_orig = NULL;
	struct GRDSLICE_SLICE *this_slice = NULL, *poly = NULL, **slice = NULL, **last = NULL;
	struct GRDSLICE_PEAK **peak = NULL, *this_peak = NULL;
	struct GRDSLICE_CTRL *Ctrl = NULL;
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) return (usage (API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if ((error = gmt_report_usage (API, options, 0, usage)) != GMT_NOERROR) bailout (error);	/* Give usage if requested */

	/* Parse the command-line arguments */

	if ((GMT = gmt_init_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_KEYS, THIS_MODULE_NEEDS, NULL, &options, &GMT_cpy)) == NULL) bailout (API->error); /* Save current state */
	if (GMT_Parse_Common (API, THIS_MODULE_OPTIONS, options)) Return (API->error);
	Ctrl = New_Ctrl (GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse (GMT, Ctrl, options)) != GMT_NOERROR) Return (error);

	/*---------------------------- This is the grdslice main code ----------------------------*/
	
	GMT_Report (API, GMT_MSG_INFORMATION, "Allocate memory and read data file\n");

	if (!strcmp (Ctrl->In.file,  "=")) {
		GMT_Report (API, GMT_MSG_ERROR, "Piping of grids not supported!\n");
		Return (EXIT_FAILURE);
	}

	gmt_M_memcpy (wesn, GMT->common.R.wesn, 4, double);	/* Current -R setting, if any */

	if ((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_ONLY, NULL, Ctrl->In.file, NULL)) == NULL) {
		Return (API->error);
	}

	if (gmt_M_is_subset (GMT, G->header, wesn)) {	/* Subset requested; make sure wesn matches header spacing */
		if ((error = gmt_M_err_fail (GMT, gmt_adjust_loose_wesn (GMT, wesn, G->header), "")))
			Return (error);
	}

	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_DATA_ONLY, wesn, Ctrl->In.file, G) == NULL) {
		Return (API->error);	/* Get subset */
	}

	if (strstr (G->header->remark, "Spherical Mercator Projected with -Jm1")) {
		/* Grid coordinates are in spherical Mercator units that we must invert before output */
		double wesn_m[4];
		/* Select plain Mercator on a sphere with -Jm1 -R0/360/-lat/+lat */
		GMT->current.setting.proj_ellipsoid = gmt_get_ellipsoid (GMT, "Sphere");
		GMT->current.proj.units_pr_degree = true;
		GMT->current.proj.pars[0] = 180.0;
		GMT->current.proj.pars[1] = 0.0;
		GMT->current.proj.pars[2] = 1.0;
		GMT->current.proj.projection = GMT->current.proj.projection_GMT = GMT_MERCATOR;
		gmt_set_geographic (GMT, GMT_IN);
		GMT->common.J.active = true;

		wesn_m[XLO] = GMT_IMG_MINLON;		wesn_m[XHI] = GMT_IMG_MAXLON;
		wesn_m[YLO] = GMT_IMG_MINLAT_80;	wesn_m[YHI] = GMT_IMG_MAXLAT_80;
		if (gmt_map_setup (GMT, wesn_m)) Return (GMT_PROJECTION_ERROR);
		gmt_geo_to_xy (GMT, 0.0, 0.0, &merc_x0, &merc_y0);
		is_mercator = true;
	}
	
	n_edges = G->header->n_rows * (urint (ceil (G->header->n_columns / 16.0)));
	edge = gmt_M_memory (GMT, NULL, n_edges, unsigned int);	/* Bit flags used to keep track of contours */

	if (!(Ctrl->Z.scale == 1.0 && Ctrl->Z.offset == 0.0)) {	/* Must transform z grid */
		GMT_Report (API, GMT_MSG_INFORMATION, "Subtracting %g and multiplying grid by %g\n", Ctrl->Z.offset, Ctrl->Z.scale);
		G->header->z_min = (G->header->z_min - Ctrl->Z.offset) * Ctrl->Z.scale;
		G->header->z_max = (G->header->z_max - Ctrl->Z.offset) * Ctrl->Z.scale;
		if (Ctrl->Z.scale < 0.0) gmt_M_double_swap (G->header->z_min, G->header->z_max);
		/* Since gmt_scale_and_offset_f applies z' = z * scale + offset we must adjust Z.offset first: */
		Ctrl->Z.offset *= Ctrl->Z.scale;
		gmt_scale_and_offset_f (GMT, G->data, G->header->size, Ctrl->Z.scale, -Ctrl->Z.offset);
	}
	if (Ctrl->L.low  > G->header->z_min) G->header->z_min = Ctrl->L.low;	/* Possibly clip the z range */
	if (Ctrl->L.high < G->header->z_max) G->header->z_max = Ctrl->L.high;

	small = Ctrl->C.interval * 1.0e-6;

	small_x = 0.01 * G->header->inc[GMT_X];	small_y = 0.01 * G->header->inc[GMT_Y];
	min_area = DBL_MAX;	max_area = -DBL_MAX;

	n_alloc = GMT_CHUNK;
	contour = gmt_M_memory (GMT, NULL, n_alloc, double);

	/* Set up contour intervals automatically from Ctrl->C.interval */
	
	min = floor (G->header->z_min / Ctrl->C.interval) * Ctrl->C.interval; if (min < G->header->z_min) min += Ctrl->C.interval;
	max = ceil  (G->header->z_max / Ctrl->C.interval) * Ctrl->C.interval; if (max > G->header->z_max) max -= Ctrl->C.interval;
	for (cs = irint (min/Ctrl->C.interval), n_contours = 0; cs <= irint (max/Ctrl->C.interval); cs++, n_contours++) {
		if (n_contours == n_alloc) {
			n_alloc += GMT_CHUNK;
			contour = gmt_M_memory (GMT, contour, n_alloc, double);
		}
		contour[n_contours] = cs * Ctrl->C.interval;
	}
	contour = gmt_M_memory (GMT, contour, n_contours, double);

	/* Because we are doing single-precision, we cannot subtract incrementally but must start with the
	 * original grid values and subtract the current contour value. */
	 
	if ((G_orig = GMT_Duplicate_Data (API, GMT_IS_GRID, GMT_DUPLICATE_DATA, G)) == NULL) {
		gmt_M_free (GMT, contour);
		Return (GMT_RUNTIME_ERROR);
	}

	/* Get initial memory allocations for slices and peaks */
	slice = gmt_M_memory (GMT, NULL, n_contours, struct GRDSLICE_SLICE *);
	last  = gmt_M_memory (GMT, NULL, n_contours, struct GRDSLICE_SLICE *);
	peak   = gmt_M_memory (GMT, NULL, n_alloc,    struct GRDSLICE_PEAK *);

	for (cs = (int)(n_contours-1); cs >= 0; cs--) {	/* For each contour value cval but starting from the top */
		c = (unsigned int)cs;

		/* Reset markers and set up new zero-contour */

		cval = contour[c];

		/* New approach to avoid round-off */

		for (ij = 0; ij < G->header->size; ij++) {
			G->data[ij] = G_orig->data[ij] - (gmt_grdfloat)cval;		/* If there are NaNs they will remain NaNs */
			if (G->data[ij] == 0.0) G->data[ij] += (gmt_grdfloat)small;	  /* There will be no actual zero-values, just -ve and +ve values */
		}

		while ((ns = gmt_contours (GMT, G, Ctrl->S.value, GMT->current.setting.interpolant, 0, edge, &begin, &x, &y)) > 0) {
			n = (uint64_t)ns;
			if (!(fabs (x[0] - x[n-1]) < small_x && fabs (y[0] - y[n-1]) < small_y)) {	/* Not a closed contour */
				gmt_M_free (GMT, x);
				gmt_M_free (GMT, y);
				continue;
			}
			x[n-1] = x[0];	y[n-1] = y[0];	/* Close the contour for sure */
			
			/* Allocate structure for the new contour slice */
			
			this_slice    = gmt_M_memory (GMT, NULL, 1, struct GRDSLICE_SLICE);
			this_slice->x = gmt_M_memory (GMT, NULL, n, double);
			this_slice->y = gmt_M_memory (GMT, NULL, n, double);
			
			/* Get information about this slice and store it in structure */
			
			this_slice->xmin = this_slice->ymin = +DBL_MAX;
			this_slice->xmax = this_slice->ymax = -DBL_MAX;
			this_slice->n = n;
			this_slice->z = cval;
			area = grdslice_contour_area_merc (x, y, n);	/* Area in Mercator^2 units; this step also ensures the contour orientation is CCW */
			gmt_M_memcpy (this_slice->x, x, n, double);	/* Copy polygon coordinates */
			gmt_M_memcpy (this_slice->y, y, n, double);
			for (i = 0; i < n; i++) {	/* Determine extreme coordinate values */
				if (x[i] < this_slice->xmin) this_slice->xmin = this_slice->x[i];
				if (x[i] > this_slice->xmax) this_slice->xmax = this_slice->x[i];
				if (y[i] < this_slice->ymin) this_slice->ymin = this_slice->y[i];
				if (y[i] > this_slice->ymax) this_slice->ymax = this_slice->y[i];
				/* Determine mean location (centroid) by getting sum */
				this_slice->x_mean += this_slice->x[i];
				this_slice->y_mean += this_slice->y[i];
			}
			/* Normalize to get mean position */
			this_slice->x_mean /= n;
			this_slice->y_mean /= n;

			if (is_mercator) { /* Mercator grid; adjust area for Mercator scale */
				gmt_xy_to_geo (GMT, &lon, &lat, this_slice->x_mean + merc_x0, this_slice->y_mean + merc_y0);
				scale = GMT->current.proj.DIST_KM_PR_DEG * cosd (lat);
			}
			else if (gmt_M_is_geographic (GMT, GMT_IN)) {	/* Geographic */
					lat = this_slice->y_mean;
					scale = GMT->current.proj.DIST_KM_PR_DEG * cosd (lat);
			}

			this_slice->area = area * (scale * scale);	/* Now in km^2 unless for Cartesian grids */

			GMT_Report (API, GMT_MSG_DEBUG, "Area = %g with scale = %g for z = %g [lat = %g]\n", area, scale, this_slice->z, lat);
			
			/* Find orientation of major/minor axes and aspect ratio */
			A[0] = A[1] = A[2] = A[3] = 0.0;
			for (i = 0; i < n; i++) {	/* Build dot-product 2x2 matrix */
				x[i] = this_slice->x[i] - this_slice->x_mean;
				y[i] = this_slice->y[i] - this_slice->y_mean;
				A[0] += x[i] * x[i];
				A[1] += x[i] * y[i];
				A[3] += y[i] * y[i];
			} 
			A[2] = A[1];
			if (gmt_jacobi (GMT, A, M, M, EigenValue, EigenVector, work1, work2, &nrots)) {	/* Solve eigen-system A = EigenVector * EigenValue * EigenVector^T */
				GMT_Report (API, GMT_MSG_WARNING, "Eigenvalue routine failed to converge in 50 sweeps.\n");
			}
			aspect = sqrt (EigenValue[0] / EigenValue[1]);	/* Major/minor aspect ratio */
			this_slice->azimuth = atan2 (EigenVector[1], EigenVector[0]) * R2D;	/* Actual, not azimuth yet - just angle CCW from horizontal */
			sincosd (this_slice->azimuth, &sa, &ca);
			this_slice->azimuth = 90.0 - this_slice->azimuth;		/* Now it is a proper azimuth */
			if (this_slice->azimuth < 0.0) this_slice->azimuth += 360.0;
			
			this_slice->minor = sqrt (this_slice->area / (aspect * M_PI));	/* Ellipse axes in km */
			this_slice->major = this_slice->minor * aspect;
			minor = sqrt (area / (aspect * M_PI));		/* Axes in map units */
			major = minor * aspect;
			
			/* Determine a measure of fit by calculating 100.0 * (1 - rms(delta_r) / major) */
			rms = 0.0;
			for (i = 0; i < n; i++) {	/* For each point */
				xr = x[i] * ca - y[i] * sa;	/* Rotate vector into eigen-system axes */
				yr = x[i] * sa + y[i] * ca;
				a = atan2 (yr, xr);		/* Angle of this vector in radians */
				r = hypot (xr, yr);		/* Magnitude of vector */
				sincos (a, &sin_a, &cos_a);
				xr = major * cos_a;		/* Ellipse model prediction */
				yr = minor * sin_a;
				r_fit = hypot (xr, yr);		/* Model radius(angle */
				dr = r - r_fit;			/* Misfit */
				rms += dr * dr;
			}
			this_slice->fit = 100.0 * (1 - sqrt (rms / n) / major);	/* Our fit parameter */
			
			/* Update information of min/max area */
			if (this_slice->area > max_area) max_area = this_slice->area;
			if (this_slice->area < min_area) min_area = this_slice->area;
			
			gmt_M_free (GMT, x);	gmt_M_free (GMT, y);	/* Free original memory returned by gmt_contours */
			
			n_slices++;
			
			/* Append this slice to end of slice array */
			
			if (last[c]) {	/* We already had at least one slice at this level - append it at end of the list */
				last[c]->next = this_slice;
				last[c] = last[c]->next;	/* Update the pointer to the last slice in this list */
			}
			else	/* First time - initialize slice and last to point to this slice */
				slice[c] = last[c] = this_slice;

			/* Determine if any previous (higher) level contours are inside the new (lower) contour */
			
			n_inside = 0;			/* Not inside anything so far */
			if (c < (n_contours-1) && slice[c+1]) {	/* There is a previous contour level and it has contours - must check for nesting */
				poly = slice[c+1];		/* First contour polygon for the previous contour level */
				while (poly) {			/* As long as there are more polygons */
					inside = true;		/* Set to true initially but this is usually reversed by one of two tests below: */
					if (poly->x[0] < this_slice->xmin || poly->x[0] > this_slice->xmax || poly->y[0] < this_slice->ymin || poly->y[0] > this_slice->ymax) {
						inside = false;		/* Outside polygon's extreme x/y-range */
					}
					else if (gmt_non_zero_winding (GMT, poly->x[0], poly->y[0], this_slice->x, this_slice->y, this_slice->n) < 2) {
						inside = false;		/* Inside rectangle, but still outside the polygon (returns 0 = outside, 1 on line (should not happen), 2 inside) */
					}
					if (inside) {	/* Here we know that the previous contour slice is inside the new one so we point to the new slice as "downstream" from the old one */
						poly->down = this_slice;
						poly->shared++;	/* Number of peaks sharing this slice */
						n_inside++;
					}

					poly = poly->next;			/* Go to next polygon in the list of contours at the previous level */
				}
			}
			if (n_inside == 0) {	/* No previous contours contained by this contour - initialize a new peak location at the center of this slice */
				this_peak = gmt_M_memory (GMT, NULL, 1, struct GRDSLICE_PEAK);
				this_peak->x = this_slice->x_mean;	/* Just use mean location for now - perhaps later choose the grid maximum */
				this_peak->y = this_slice->y_mean;
				this_peak->z = this_slice->z;		/* Likewise - this might eventually be the local high value in the grid */
				this_peak->start = this_slice;		/* This is the top slice in the stack below this point */
				this_peak->id = n_peaks;			/* This is the unique peak ID */
				peak[n_peaks] = this_peak;			/* Add peak to peak array */
				n_peaks++;
				if (n_peaks == n_alloc) {	/* Get more memory */
					n_alloc += GMT_CHUNK;
					peak = gmt_M_memory (GMT, peak, n_alloc, struct GRDSLICE_PEAK *);
				}
			}

		}
		GMT_Report (API, GMT_MSG_INFORMATION, "Tracing the %8.2f contour: # of slices: %6d # of peaks: %6d\n", cval, n_slices, n_peaks);
	}
	
	peak = gmt_M_memory (GMT, peak, n_peaks, struct GRDSLICE_PEAK *);

	if (Ctrl->A.active) {	/* Calculate distances between peaks */
		int which;
		double dist;
		skip = gmt_M_memory (GMT, NULL, n_peaks, bool);
		for (c = 0; c < n_peaks; c++) {	/* For each peak */
			for (i = c+1; i < n_peaks; i++) {	/* For each other peak */
				dist = 0.001 * gmt_great_circle_dist_meter (GMT, peak[c]->start->x_mean, peak[c]->start->y_mean, peak[i]->start->x_mean, peak[i]->start->y_mean);
				if (dist < Ctrl->A.cutoff) {	/* These two peaks are too close */
					which = (peak[c]->start->z >= peak[i]->start->z) ? i : c;
					skip[which] = true;
				}
			}
		}
	}
	
	sprintf (file, "%s_pos.txt", Ctrl->D.file);
	if ((fp = fopen (file, "w")) == NULL) {
		GMT_Report (API, GMT_MSG_ERROR, "Unable to create file %s\n", file);
		Return (EXIT_FAILURE);
	}
	for (c = n_skipped = 0; c < n_peaks; c++) {	/* For each peak */
		if (Ctrl->A.active && skip[c]) {
			n_skipped++;
		}
		else {
			fprintf (fp, "> -Z%g -L%d x y z id area major minor azimuth fit \n", peak[c]->z, peak[c]->id);
			poly = peak[c]->start;		/* First contour polygon originating form this peak */
			while (poly) {			/* As long as there are more polygons at this level */
				out[0] = poly->x_mean;	out[1] = poly->y_mean;	out[2] = poly->z;
				out[3] = peak[c]->id; out[4] = poly->area; out[5] = poly->major;
				out[6] = poly->minor; out[7] = poly->azimuth; out[8] = poly->fit;
				if (is_mercator) gmt_xy_to_geo (GMT, &out[0], &out[1], out[0] + merc_x0, out[1] + merc_y0);	/* Get lon, lat */
				GMT->current.io.output (GMT, fp, 9, out, NULL);
				poly = poly->down;		/* Go to next polygon in this contour list */
			}
		}
		gmt_M_free (GMT, peak[c]);		/* Free memory as we go */
	}
	fclose (fp);
	
	sprintf (file, "%s_slices.txt", Ctrl->D.file);
	if ((fp = fopen (file, "w")) == NULL) {
		GMT_Report (API, GMT_MSG_ERROR, "Unable to create file %s\n", file);
		Return (EXIT_FAILURE);
	}
	
	if (Ctrl->T.active) {
		sprintf (bfile, "%s_bottom.txt", Ctrl->D.file);
		if ((tp = fopen (bfile, "w")) == NULL) {
			GMT_Report (API, GMT_MSG_ERROR, "Unable to create file %s\n", bfile);
			Return (EXIT_FAILURE);
		}
		sprintf (kfile, "%s_indices.txt", Ctrl->D.file);
		if ((kp = fopen (kfile, "w")) == NULL) {
			GMT_Report (API, GMT_MSG_ERROR, "Unable to create file %s\n", kfile);
			Return (EXIT_FAILURE);
		}
	}
	
	for (c = 0; c < n_contours; c++) {	/* For each contour value cval */
		poly = slice[c];		/* First contour polygon at this contour level */
		while (poly) {			/* As long as there are more polygons at this level */
			/* Output a multisegment header and body */
			out[0] = poly->x_mean;	out[1] = poly->y_mean;
			if (is_mercator) gmt_xy_to_geo (GMT, &out[0], &out[1], out[0] + merc_x0, out[1] + merc_y0);	/* Get lon, lat */
			fprintf (fp, "> -Z%g -L%g -N%d -S%g/%g/%g/%g/%g/%g\n", poly->area, poly->z, poly->shared, out[0], out[1], poly->azimuth, poly->major, poly->minor, poly->fit);
			if (Ctrl->T.active && doubleAlmostEqualZero (poly->z, Ctrl->T.blevel) && poly->area >= Ctrl->T.acutoff) {
				kp_id++;
				fprintf (tp, "> %d -Z%g -L%g -N%d -S%g/%g/%g/%g/%g/%g\n", kp_id, poly->area, poly->z, poly->shared, out[0], out[1], poly->azimuth, poly->major, poly->minor, poly->fit);
				/* indices file format: lon lat x y id */
				gmt_xy_to_geo (GMT, &out[0], &out[1], out[0] + merc_x0, out[1] + merc_y0);
				fprintf (kp, "%g\t%g\t%g\t%g\t%d\n", out[0], out[1], poly->x_mean, poly->y_mean, kp_id);
			}
			for (i = 0; i < poly->n; i++) {
				out[0] = poly->x[i];	out[1] = poly->y[i];	out[2] = poly->z;
				if (is_mercator) gmt_xy_to_geo (GMT, &out[0], &out[1], out[0] + merc_x0, out[1] + merc_y0);	/* Get lon, lat */
				GMT->current.io.output (GMT, fp, 3, out, NULL);
				if (Ctrl->T.active && doubleAlmostEqualZero (poly->z, Ctrl->T.blevel) && poly->area >= Ctrl->T.acutoff)
					GMT->current.io.output (GMT, tp, 3, out, NULL);
			}
			this_slice = poly;
			poly = poly->next;		/* Go to next polygon in this contour list */
            
			gmt_M_free (GMT, this_slice->x);
			gmt_M_free (GMT, this_slice->y);
			gmt_M_free (GMT, this_slice);	/* Free memory as we go */
		}
	}
    fclose (fp);

	if (Ctrl->T.active) {
		fclose (tp);
		fclose (kp);
		if (kp_id == 0) {
			GMT_Report (API, GMT_MSG_WARNING, "No peak bottoms or indices written to %s and %s despite -T being set\n", bfile, kfile);
			gmt_remove_file (GMT, bfile);
			gmt_remove_file (GMT, kfile);
		}
	}

	gmt_M_free (GMT, edge);
	gmt_M_free (GMT, contour);
	gmt_M_free (GMT, slice);
	gmt_M_free (GMT, last);
	gmt_M_free (GMT, peak);
	if (Ctrl->A.active) gmt_M_free (GMT, skip);

	GMT_Report (API, GMT_MSG_INFORMATION, "Done, min/max areas: %g %g\n", min_area, max_area);
	GMT_Report (API, GMT_MSG_INFORMATION, "Data written to %s_slices.txt and %s_pos.txt\n", Ctrl->D.file, Ctrl->D.file);

	Return (GMT_NOERROR);
}
