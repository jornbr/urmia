These files were provided by Eisa Bozorgzadeh by e-mail on 15/10/2014,
with the following comments:

----

You may find it interesting that the bathymetry of Uremia Lake has been carried out already. Here are following informations: 
1- Uremia lake Bathymetry DEM 
2- Excel sheet including Volume-Area- Depth 
The Bathometry file was derived from the bathymetry operation in 2011. 
Please know that I have no Idea about the reliability of bathymetry file. You may need to check this file to find how reliable it is. 
I look forward having your valuable comments on the presented data. 
Best Regards 
Eisa 

----

The Python script parse_bathymetry.py parses the bathymetry, computes, plots
and saves bottom level-surface area data (hypsograph.dat), and plots the bathymetry.

Issues with the provided bathymetry:

(1) It is unknown what reference level (datum) was used for the elevation values.
The maximum elevation value is 1278.0 m in the bathymetry, while the maximum surface elevation
of Lake Urmia in the time series provided by Dr. Karimi in Iran exceeds 1278.4 m. Either the bathymetry
does not extend far enough in  the horizontal, and as a result excludes the higher elevations
that at times have been lake-covered, or the current bathymetry and the time series use
different references for "mean sea level".

(2) If you plot a histogram of the elevations in the bathymetry file, you see obvious peaks at
whole meter levels (e.g., 1270.0, 1271.0, 1272.0, etc. - se BoxModel/interpolated_area_dist.png)
This is very likely an artifact of the way the bathymetry was produced. E.g., it could hint at 
interpolation from contour lines.

Jorn Bruggeman, 6/1/2015