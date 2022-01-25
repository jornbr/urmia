Box model of Lake Urmia, describing the water budget and salinity in the lake.
It considers the bathymetry implicitly through a hypsograph [cumulative distribution
of surface area as function of bottom elevation], and also considers salt precipitation
and its effect on bottom morphology.

boxmodel.py is the main script to run the model. The first lines in this
file contain configuration settings. The script uses the following inputs:

- UrLa.dat: hypsograph and time series of runoff, evaporation, precipitation, water level,
  provided by Dr. Karimi during our stay in Iran.
- bathymetry.dat: high resolution hypsograph produced from bathymetry provided by Eisa.
  see Data/Bathymetry
- water_budget.pickle: water budget data derived from UrLa.dat, stored in binary
  "pickle" format to enable quick reads by boxmodel.py. The script that generates this file
  from UrLa.dat is parse_water_budget.py.

The script boxmodel.py generates several plots (*.png) and can also create animation stills (animate/*.png)
if the "animate" option is set to True in boxmodel.py.

Sample animations derived from these are provided as movie_bathymetry1.avi and movie_bathymetry2.avi.
The former uses the low-res hypsograph from UrLa.dat; the latter the high-res data from hypsograph.dat.
Note that the salt layer that appears late in these animations is equivalent to *solid salt*, that is, it
describes the volume taken up by pure salt. The observed salt layer will be thicker because it is porous,
and contains a mixture of solid salt and water.

To run, you need Python with the numpy and matplotlib modules. On Windows, it is easiest to download
a ready-made Python distribution such as the Enthought Python Distribution (https://www.enthought.com/products/epd/)
or Anaconda (http://continuum.io).