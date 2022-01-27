from osgeo import gdal,osr
import pyproj
import numpy
import pylab

# Open bathymetry file [uses GDAL]
ds = gdal.Open('Bathimetry.img')

# Read first variable
band = ds.GetRasterBand(1)

# Convert to array with negative [invalid] elevation values masked out.
levels = band.ReadAsArray()
levels = numpy.ma.array(levels,mask=levels<0)

# Obtain information on the x,y grid
geotransform = ds.GetGeoTransform()
print('geo transform: %s' % (geotransform,))
originX = geotransform[0]
originY = geotransform[3]
pixelWidth = geotransform[1]
pixelHeight = geotransform[5]

# Assume x,y units are meter and that variation in surface area per grid cell is negligible.
cell_area = abs(pixelWidth)*abs(pixelHeight)
print('Surface area of a single grid cell = %s m2' % cell_area)

# Convert elevations to 1D array with masked values removed,
# sorted from low [deep] to high [shallow].
print('Computing hypsograph...')
levels_flat = numpy.sort(levels.compressed())
area_flat = numpy.empty_like(levels_flat)
area_flat[:] = cell_area
cumarea_flat = area_flat.cumsum()

print('Saving raw hypsograph...')
with open('hypsograph_all.dat','w') as f:
    for level,cumarea in zip(levels_flat,cumarea_flat):
        f.write('%s\t%s\n' % (level,cumarea))

# Write hypsograph [cumulative distribution of surface area over bottom elevation]
print('Saving hypsograph...')
with open('hypsograph.dat','w') as f:
    maxvals = numpy.linspace(levels_flat[0],levels_flat[-1],1000)
    cumareas = numpy.interp(maxvals,levels_flat,cumarea_flat)
    for maxval,cumarea in zip(maxvals,cumareas):
        f.write('%s\t%s\n' % (maxval,cumarea))

# Plot hypsograph
print('Plotting hypsograph...')
pylab.plot(levels_flat,cumarea_flat,'-')
pylab.grid(True)
pylab.xlabel('bottom elevation (m)')
pylab.ylabel('cumulative surface area (m2)')

# Compute x,y corner coordinates of bathymetry
# NB GDAL specs state upper left corner of the upper left pixel is at originX,originY
x = originX + numpy.arange(levels.shape[1]+1)*pixelWidth
y = originY + numpy.arange(levels.shape[0]+1)*pixelHeight

# From 1D to 2D coordinate arrays
x2d,y2d = pylab.meshgrid(x,y)

# Determine map projection that can e used to convert x,y to lon,lat
s_srs_obj = osr.SpatialReference()
s_srs_obj.SetFromUserInput(ds.GetProjectionRef())
p = pyproj.Proj(s_srs_obj.ExportToProj4())

# Compute longitude, latitude
print('Computing longitude and latitude from x,y...')
lon,lat = p(x2d,y2d,inverse=True)
print('   longitude range: %s - %s' % (lon.min(),lon.max()))
print('   latitude range: %s - %s' % (lat.min(),lat.max()))

# Create KML file for Google Earth
print('Saving KML + map...')
lonrange = lon.max()-lon.min()
latrange = lat.max()-lat.min()
fig = pylab.figure(figsize=(lonrange*10,latrange*10))
fig.subplots_adjust(left=0., bottom=0., right=1., top=1.)
ax = fig.gca()
pylab.pcolormesh(lon,lat,levels)
pylab.axis('tight')
ax.patch.set_visible(False)
fig.patch.set_visible(False)
ax.axis('off')
pylab.savefig('bath.png',dpi=300)
pylab.close(fig)
with open('bath.kml','w') as f: f.write("""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
    <GroundOverlay>
        <name>%s</name>
        <Icon>
            <href>%s</href>
        </Icon>
        <LatLonBox>
            <north>%s</north>
            <south>%s</south>
            <east>%s</east>
            <west>%s</west>
            <rotation>0.</rotation>
        </LatLonBox>
    </GroundOverlay>
</kml>""" % ('Lake Urmia bathymetry','bath.png',lat.max(),lat.min(),lon.max(),lon.min()))

# Save NetCDf
print('Saving to NetCDF...')
import netCDF4
nc = netCDF4.Dataset('bath.nc','w')
nc.createDimension('xc',x.size-1)
nc.createDimension('yc',y.size-1)
nc.createDimension('x',x.size)
nc.createDimension('y',y.size)
nc.createVariable('x','f8',('x',))[:] = x
nc.createVariable('y','f8',('y',))[:] = y
nc.createVariable('xc','f8',('xc',))[:] = (x[:-1]+x[1:])/2
nc.createVariable('yc','f8',('yc',))[:] = (y[:-1]+y[1:])/2
nc.createVariable('lon','f8',('y','x'))[:,:] = lon
nc.createVariable('lat','f8',('y','x'))[:,:] = lat
nclevel = nc.createVariable('level','f4',('yc','xc'),fill_value=band.GetNoDataValue())
nclevel[:,:] = levels
nclevel.units = 'm'
nclevel.long_name = 'bottom elevation'
nc.close()

# Plot bathymetry
print('Plotting bathymetry...')
pylab.figure()
pylab.subplot(1,2,1)
pc = pylab.pcolormesh(x2d,y2d,levels)
pylab.colorbar(pc)
pylab.axis('equal')
pylab.subplot(1,2,2)
pc = pylab.pcolormesh(lon,lat,levels)
pylab.colorbar(pc)
pylab.axis('equal')
pylab.show()
