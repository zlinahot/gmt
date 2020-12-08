.. index:: ! grdslice
.. include:: ../module_supplements_purpose.rst_

********
grdslice
********

|grdslice_purpose|

Synopsis
--------

.. include:: ../../common_SYN_OPTs.rst_

**gmt grdslice**
*grid*
|-C|\ *interval*
[ |-A| ]
[ |-D|\ [*prefix*] ]
[ |-L|\ *low/high* ]
[ |-S|\ *smoothfactor* ]
[ |-T|\ *bottom_level*/*area_cutoff* ]
[ |SYN_OPT-V| ]
[ |-Z|\ [**+s**\ *factor*][**+o**\ *shift*] ]
[ |SYN_OPT-f| ]
[ |SYN_OPT--| ]

|No-spaces|

Description
-----------

**grdslice** is used to determine local *peaks* in the grid.  To do so, it reads a binary 2-D
grid and contour-slices the surface using the selected constant contour interval,
starting from the max value and moving down to lower and lower values.  As new contours are found,
their centers are considered the location of a new peak, and we trace the hierarchy of how
each new slice belongs to groups of previously detected peaks.  Slice shapes are equated to ellipses
and their major and major axes and orientations are computed. Information are written to two separate
data files.

Required Arguments
------------------

*grid*
    2-D gridded data set to be sliced. (See GRID FILE FORMATS below). We automatically detect if this
    grid is geographic or in spherical Mercator units (such as created by :doc:`img2grd </supplements/img/img2grd>`.
    If plain Cartesian you may need to use modifiers to convert horizontal distance units to km.

.. _-C:

**-C**\ *interval*
    Contour *interval* for slicing.  The smaller, the more detailed the results and the longer
    the execution time.

Optional Arguments
------------------

.. _-A:

**-A**
    Examine the peak locations found and skip those that fall within *cutoff* km from another
    peak with a larger amplitude [no skipping].

.. _-D:

**-D**\ *prefix*
    Set the filename *prefix* for the two output files (*prefix*\ _slices.txt and *prefix*\ _pos.txt) [PEAK].

.. _-L:

**-L**\ *low/high*
    Limit range: Only do the contour slicing within the specified range of contours [use he entire range].

.. _-S:

**-S**\ *smoothfactor*
    Resample the contour lines at roughly every (gridbox_size/*smoothfactor*) interval [no resampling].

.. _-T:

**-T**\ *bottom_level*/*area_cutoff*
    Specify the bottom level of contouring as *bottom_level* and ignore contours whose area are less
    than *area_cutoff* in km^2 (or user units for Cartesian grids).  With **-D**, we will write
    *prefix*\ _bottom.txt and *prefix*\ _indices.txt files as well. If **-L** is used, the *bottom_level*
    must be equal or larger than the *low* value.

.. _-Z:

**-Z**\ [**+s**\ *factor*][**+o**\ *shift*]
    Use to subtract *shift* from the data and multiply the results by
    *factor* before contouring starts [1/0]. (Numbers in **-C**, **-L**
    and **-T** refer to values *after* this scaling has occurred.)

.. |Add_-f| unicode:: 0x20 .. just an invisible code
.. include:: ../../explain_-f.rst_

.. |Add_-V| unicode:: 0x20 .. just an invisible code
..  include:: ../../explain_-V.rst_

.. include:: ../../explain_colon.rst_

.. include:: ../../explain_help.rst_

Examples
--------

To slice the grid gravity_smts.grd every 0.1 mGal and write results to files called X_*, try:

::

    grdslice gravity_smts.grd -C0.1 -DX -V

References
----------

Kim, S.-S., and P. Wessel (2011), New global seamount census from altimetry-derived gravity data,
*Geophys. J. Int., 186*, 615-631,
`http://dx.doi.org/10.1111/j.1365-246X.2011.05076.x <http://dx.doi.org/10.1111/j.1365-246X.2011.05076.x>`_.

Wessel, P., D. T. Sandwell, and S.-S. Kim (2010), The global seamount census, *Oceanography, 23(1)*, 24-33,
`http://dx.doi.org/10.5670/oceanog.2010.60 <http://dx.doi.org/10.5670/oceanog.2010.60>`_.

See Also
--------

:doc:`gmt.conf </gmt.conf>`, :doc:`gmt </gmt>`,
:doc:`grdcontour </grdcontour>`,
:doc:`grdmath </grdmath>`,
:doc:`grdseamount </supplements/potential/grdseamount>`,
:doc:`img2grd </supplements/img/img2grd>`
