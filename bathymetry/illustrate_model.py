print('Reading hypsograph...')
levels = []
areas = []
with open('hypsograph_all.dat','rU') as f:
    for l in f:
        level,area = map(float,l.strip('\n').split('\t'))
        levels.append(level)
        areas.append(area)

import numpy
levels = numpy.array(levels)
areas = numpy.array(areas)

level_grid = numpy.linspace(levels.min(),levels.max(),50)
level_cent = (level_grid[1:]+level_grid[:-1])/2
area_grid = numpy.interp(level_grid,levels,areas)

area_grid /= 1e6

rock = levels.min()-1
water = levels.max()

import pylab
for i,z in enumerate(level_cent):
    xs = (area_grid[i],area_grid[i+1],area_grid[i+1],area_grid[i],area_grid[i])
    ys = z,z,z-1,z-1,z
    pylab.fill(xs,ys,'w',edgecolor='k')
    ys = rock,rock,z-1,z-1,rock
    pylab.fill(xs,ys,'gray',edgecolor='none')
    ys = z,z,water,water,z
    pylab.fill(xs,ys,'b',edgecolor='none')
pylab.plot((0.,area_grid[-1]),(water,water),'-k')
pylab.xlim(0.,area_grid[-1])
pylab.ylim(rock,None)
pylab.xlabel('surface area (km2)')
pylab.ylabel('elevation (m)')
pylab.grid()
pylab.savefig('boxmodel.png',dpi=300)
