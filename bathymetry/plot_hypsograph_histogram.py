print('Reading hypsograph...')
levels = []
with open('hypsograph_all.dat','rU') as f:
    for l in f:
        level,area = map(float,l.strip('\n').split('\t'))
        levels.append(level)

import numpy
levels = numpy.array(levels)

print('Plotting histogram...')
import pylab
pylab.hist(levels,200,edgecolor='w',color='k',linewidth=0.3)
pylab.xlabel('bottom elevation (m)')
pylab.ylabel('frequency (#)')
pylab.xlim(levels.min(),levels.max())
pylab.grid(True)
pylab.savefig('bathymetric_depth_histogram.png',dpi=300)
pylab.show()

