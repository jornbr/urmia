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
#print 'Computing longitude and latitude from x,y...'
#lon,lat = p(x2d,y2d,inverse=True)
#print '   longitude range: %s - %s' % (lon.min(),lon.max())
#print '   latitude range: %s - %s' % (lat.min(),lat.max())

x_2002 = numpy.array(open('./2002/xx.txt','rU').read().split(),dtype=float)
y_2002 = numpy.array(open('./2002/yy.txt','rU').read().split(),dtype=float)
z_2002 = numpy.array(open('./2002/zz.txt','rU').read().split(),dtype=float)

import scipy.interpolate

xc,yc = (x[1:]+x[:-1])/2,(y[1:]+y[:-1])/2

delta=2000.
xgrid_if = numpy.arange(x_2002.min(),x_2002.max(),delta*2)
ygrid_if = numpy.arange(y_2002.min(),y_2002.max(),delta*2)
xgrid = (xgrid_if[1:]+xgrid_if[:-1])/2
ygrid = (ygrid_if[1:]+ygrid_if[:-1])/2
z_old_ip = numpy.ma.array(numpy.empty((xgrid.size,ygrid.size)),mask=True)
z_new_ip = numpy.ma.array(numpy.empty((xgrid.size,ygrid.size)),mask=True)
for ix,xold in enumerate(xgrid):
    for iy,yold in enumerate(ygrid):
        ixmin,ixmax = numpy.searchsorted(xc,(xold-delta,xold+delta))
        iymin,iymax = -numpy.searchsorted(yc[::-1],(yold+delta,yold-delta))
        old_valid = numpy.logical_and(numpy.logical_and(x_2002>xold-delta,x_2002<xold+delta),numpy.logical_and(y_2002>yold-delta,y_2002<yold+delta))
        validlevels = levels[iymin:iymax,ixmin:ixmax]
        #if validlevels.mask.any(): continue
        if validlevels.mask.sum()/float(validlevels.size)<0.5:
            z_new_ip[ix,iy] = validlevels.mean()
        if old_valid.sum()>10:
            z_old_ip[ix,iy] = z_2002[old_valid].mean()

z_old_ip_flat = z_old_ip.flatten()
deltaz_flat = (z_new_ip-z_old_ip).flatten()
valid = numpy.logical_not(deltaz_flat.mask)
z_old_ip_flat = z_old_ip_flat[valid]
deltaz_flat = deltaz_flat[valid]

delta_z = .5
z_old_rm = numpy.linspace(z_old_ip_flat.min()+delta_z,z_old_ip_flat.max()-delta_z,100)
dz = numpy.ma.array(numpy.empty_like(z_old_rm),mask=True)
for i,z in enumerate(z_old_rm):
    valid = numpy.logical_and(z_old_ip_flat>z-delta_z,z_old_ip_flat<z+delta_z)
    #if valid.sum()<10: continue
    dz[i] = deltaz_flat[valid].mean()

# offset old bathymetry such that the minimum running-mean-of-differences is 0.
datum_old = dz.min()
print(datum_old)

with open('bathymetry_difference.dat', 'w') as f:
    f.write('bottom depth in 2002 (m)\tdifference 2011-2002 (m)\n')
    for values in sorted(zip(z_old_ip_flat+datum_old, deltaz_flat-datum_old), key=lambda x: x[0]):
        f.write('%s\t%s\n' % values)

pylab.figure(figsize=(5,5))
pylab.plot(z_old_ip_flat+datum_old,deltaz_flat-datum_old,'o')
pylab.plot(z_old_rm+datum_old,dz-datum_old,'-r',linewidth=2.)
pylab.grid(True)
pylab.gca().ticklabel_format(useOffset=False)
pylab.title('bottom increase with %.1f m running mean' % (2*delta_z))
pylab.xlabel('bottom elevation in 2002 (m)')
pylab.ylabel('2002-2011 increase in bottom elevation (m)')
pylab.savefig('difference_trend.png',dpi=150)

z_min = min(z_old_ip.min()+datum_old,z_new_ip.min())
z_max = max(z_old_ip.max()+datum_old,z_new_ip.max())
pylab.figure(figsize=(10,8))
ax = pylab.subplot2grid((2,4),(0, 0),colspan=2)
pylab.title('2002')
pylab.pcolormesh(xgrid_if,ygrid_if,z_old_ip.T+datum_old,vmin=z_min,vmax=z_max)
ax = pylab.subplot2grid((2,4),(0, 2),colspan=2)
pylab.title('2011')
pylab.pcolormesh(xgrid_if,ygrid_if,z_new_ip.T,vmin=z_min,vmax=z_max)
ax = pylab.subplot2grid((2,4),(1, 1),colspan=2)
pylab.title('difference: 2011-2002')
pc=pylab.pcolormesh(xgrid_if,ygrid_if,(z_new_ip-z_old_ip-datum_old).T)
cb = pylab.colorbar(pc)
cb.set_label('difference (m)')
pylab.savefig('bathymetries.png',dpi=150)

pylab.show()
sdds

# Plot bathymetry
print('Plotting bathymetry...')
pylab.figure()
pylab.subplot(1,2,1)
pc = pylab.pcolormesh(x2d,y2d,levels)
pylab.colorbar(pc)
pylab.axis('equal')
pylab.subplot(1,2,2)
pylab.scatter(x_2002,y_2002,c=z_2002,linewidths=0)

pylab.show()
