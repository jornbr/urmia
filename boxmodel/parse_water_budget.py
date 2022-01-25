import numpy
import pylab
import datetime

data = []
dts = []

with open('UrLa.dat','rU') as f:
    # Skip hypsograph data.
    while 1:
        l = f.readline()
        if not l.strip(): break

    # Skip column labels.
    f.readline()
    f.readline()

    # Read data
    for iline,l in enumerate(f):
        y,m,d,level,runoff,rain,evap = l.split()[3:]
        dts.append(datetime.datetime(*map(int,(y,m,d))))
        data.append((float(runoff),float(rain),float(evap),float(level)))
        if (iline+1)%1000==0: print('Line %i processed.' % (iline+1))

dts = pylab.date2num(dts)
data = numpy.array(data)

import pickle

with open('water_budget.pickle','wb') as f:
    pickle.dump(dts,f,-1)
    pickle.dump(data,f,-1)

