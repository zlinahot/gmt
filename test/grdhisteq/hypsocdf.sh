gmt begin hypsocdf ps
	gmt grdhisteq @earth_relief_30m -C11 -D -Fcdf.txt > levels
	gmt plot cdf.txt -R-10000/6000/0/1.1 -JX18c -Bxaf+l"Depth" -Byafg+lCDF -BWSrt -W1p
	awk '{printf "> \n%s 0\n%s 1.1\n", $1, $1}' levels | gmt plot -W0.5p
gmt end show

