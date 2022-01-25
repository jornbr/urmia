# ------------------------
# Configuration settings
# ------------------------

# First year of simulation
startyear = 1994

# Salinity in kg/m3 (g/L) at start of first year of simulation
# In presentation by Fakhimi Kavossi Heydarzadeh, 14/12/2014:
# 2007: 375 g/L
# 2008: 382 g/L
# 2009: 414 g/L
# 2010: 423 g/L
# 2011: 425 g/L
# 2012: 471 g/L
# Dahesht et al. 2010 Iran J Fisheries Science:
# 1987: 152-168 ppt
# 2005-2006: >300 ppt
# Van Stappen et al 2001
# 1996: 170 ppt
salinity_ini = 130.

# Source of bathymetry/hypsograph
# 1: data provided by Karimi in Iran
# 2: data extracted from bathymetry (has max elevation < max observed surface level, causing occasional water overflow)
bathymetry_source = 2

# Hypsograph interpolation method
# 1: interpolate PDF of surface area
# 2: interpolate cumulative distribution of surface area
hypsograph_interpolation = 2

# NaCl solubility product in (mol/L)^2 - Krumgalz & Millero 1989 Mar Chem
K_NaCl = 37.19

# Molar mass of salt constituents
Na_molweight = 22.9898
Cl_molweight = 35.453

# Crude estimates for NaCl precipitation and dissolution rates (kg/m3/d).
# NB dissolution not implemented yet
NaCl_precip_rate = 20.
NaCl_dissol_rate = 20.

# Density of precipitated NaCl (kg/m3) - http://en.wikipedia.org/wiki/Sodium_chloride
dens_NaCl = 2165.

# Porosity of precipitated salt, used to calculate thickness of visible (porous) salt layer.
salt_porosity = 0.5

# Whether to generate animation stills.
animate = False

# ------------------------
# Code
# ------------------------

import datetime
import pickle
import numpy
import pylab,matplotlib.dates

# Read forcing (evaporation, precipitation, rivers) from previously parsed data.
with open('water_budget.pickle','rb') as f:
    water_budget_dt = pickle.load(f)
    water_budget = pickle.load(f)

# Read basin shape (hypsograph): cumulative surface area (m2) as function of bottom elevation (m)
if bathymetry_source==1:
    # Hypsograph data provided by Karimi in Iran
    with open('UrLa.dat','rU') as f:
        labels = f.readline()
        data = []
        for l in f:
            l = l.strip()
            if not l: break
            level,area,volume = map(float,l.split())
            data.append((level,area,volume))

    # Convert to SI units
    data = numpy.array(data)
    data[:,1]*=1e6  # surface area: from km2 to m2
    data[:,2]*=1e6  # volume: from 10^-3 km3 to m3

    # Recompute cumulative volumes (as function of bottom elevation)
    # and compare with original value to test whether we understand the hypsograph data.
    if False:
        cumv = 0.
        print('Hypsograph validation:')
        for i in range(1,data.shape[0]):
            # The volume added by the current layer consists of a part added on top of previous water layers,
            # which has a thickness equal to the difference in elevation between the left and right bin boundary,
            # and a part on top of an area that was previously dry. Assume the mean elevation of the dry area
            # is the average of the elevation between left and right boundaries, which means the water layer thickness
            # is half the difference between left and right elevation. The two contibutions equate to:
            #   data[i-1,1]*(data[i,0]-data[i-1,0]) + (data[i,1]-data[i-1,1])*(data[i,0]-data[i-1,0])/2
            # After rewriting:
            cumv += (data[i,1]+data[i-1,1])*(data[i,0]-data[i-1,0])/2
            print('  ori = %s, reconstr = %s, rel dev = %s' % (data[i,2],cumv,(cumv-data[i,2])/data[i,2]))
else:
    # Hypsograph data extracted from bathymetry
    with open('hypsograph.dat','rU') as f:
        data = numpy.array([list(map(float,l.split())) for l in f])

salt_thickness = numpy.maximum(-(data[:,0]-1271.)*.5*salt_porosity,0.)
#data[:,0] -= salt_thickness
max_z_surf = data[:,0].max()

# Plot hypsograph
pylab.figure()
pylab.plot(data[:,0],data[:,1],'o')
pylab.xlabel('level (m)')
pylab.ylabel('cumulative surface area (km2)')
pylab.grid(True)
pylab.savefig('cum_area_dist.png',dpi=150)

# Define new (finer) elevation grid to interpolate hypsograph onto.
z_grid_if = numpy.linspace(data[0,0],data[-1,0],100)
z_grid_cent = (z_grid_if[1:]+z_grid_if[:-1])/2
z_grid_delta = z_grid_if[1:]-z_grid_if[:-1]

if hypsograph_interpolation==1:
    # Interpolate PDF of surface area
    x_cent = (data[1:,0]+data[:-1,0])/2             # centre of original bins
    area_per_bin = numpy.diff(data[:,1])            # bin-specific surface area (m2)
    y_cent = area_per_bin/(data[1:,0]-data[:-1,0])  # bin-size independent PDF  (m)
    area_dist = numpy.interp(z_grid_cent,x_cent,y_cent)*(z_grid_if[1:]-z_grid_if[:-1])
else:
    # Interpolate cumulative distribution of surface area
    cumarea_fine = numpy.interp(z_grid_if,data[:,0],data[:,1])
    area_dist = numpy.diff(cumarea_fine)

def compute_cumulative_volume(z_grid_delta,cumarea):
    addedvol = numpy.empty_like(cumarea)
    addedvol[0] = 0.
    addedvol[1:] = (cumarea[1:]+cumarea[:-1])*z_grid_delta[1:]/2
    return addedvol.cumsum()

def level_and_area_from_volume(area_dist,z_grid_cent,z_grid_delta,V):
    cumvol = compute_cumulative_volume(z_grid_delta,area_dist.cumsum())
    return numpy.interp(V,cumvol,z_grid_cent),numpy.interp(V,cumvol,area_dist.cumsum())

def volume_from_level(area_dist,z_grid_cent,z_grid_delta,level):
    cumvol = compute_cumulative_volume(z_grid_delta,area_dist.cumsum())
    return numpy.interp(level,z_grid_cent,cumvol)

# Plot original and interpolated distribution of surface area (as function of bottom elevation).
pylab.figure()
pylab.plot(z_grid_cent,area_dist/(z_grid_if[1:]-z_grid_if[:-1]),'-')
pylab.plot((data[1:,0]+data[:-1,0])/2,numpy.diff(data[:,1])/(data[1:,0]-data[:-1,0]),'o')
pylab.xlabel('level (m)')
pylab.ylabel('surface area per level (m)')
pylab.savefig('interpolated_area_dist.png',dpi=150)

# Validate interpolated hypsograph
if data.shape[1]>2:
    print('Validate interpolated hypsograph against original cumulative volume data:')
    for i in range(data.shape[0]):
        print('  ',volume_from_level(area_dist,z_grid_cent,z_grid_delta,data[i,0]),data[i,2])

# Precipitated salt (kg/m2) as a function of bottom elevation
salt_per_level = numpy.zeros_like(z_grid_cent)

# Find elevation of lake surface at starting time.
startdt = pylab.date2num(datetime.datetime(startyear,1,1))
z_surf_ini = numpy.interp(startdt,water_budget_dt,water_budget[:,3])

# Infer initial volume from initial elevation (and hypsograph)
V = volume_from_level(area_dist,z_grid_cent,z_grid_delta,z_surf_ini)
S = V*salinity_ini

print('Initial level = %f m, initial volume = %.3f km3' % (z_surf_ini,V/1e9))

def printStill(z_bot,area,h_salt,z_surf,dts,z_surfs,salts,path,curdate,dpi=96):
    # Plot basin.
    cumarea = area.cumsum()
    radius = numpy.sqrt(cumarea/numpy.pi)  # distribute area over circle and compute effective radius
    x = numpy.concatenate((-radius[::-1],radius))/1000.
    y_bot = numpy.concatenate((z_bot[::-1],z_bot))
    z_salt = z_bot + h_salt
    y_salt = numpy.concatenate((z_salt[::-1],z_salt))
    y_surf = numpy.empty_like(x)
    y_surf[:] = numpy.maximum(z_surf,y_salt)
    y_figbot = numpy.empty_like(x)
    y_figbot[:] = y_bot.min()
    pylab.subplot2grid((2,2), (0, 0), rowspan=2)
    pylab.cla()
    pylab.fill_between(x,y_figbot,y_bot,facecolor='gray')
    pylab.fill_between(x,y_bot,y_salt,facecolor='w')
    pylab.fill_between(x,y_salt,y_surf,facecolor='b')
    pylab.plot(x,y_bot,'-k')
    pylab.plot(x,y_salt,'-k')
    pylab.plot(x,y_surf,'-k')
    pylab.xlabel('distance from centre (km)')
    pylab.ylabel('elevation (m)')
    pylab.axis('tight')
    strdate = curdate.strftime('%d-%m-%Y')
    pylab.title(strdate)

    # Plot time series for water level.
    locator = matplotlib.dates.YearLocator(4)
    ax = pylab.subplot2grid((2,2), (0, 1))
    pylab.cla()
    istart = water_budget_dt.searchsorted(pylab.date2num(startdate))-1
    pylab.plot_date(water_budget_dt[istart:],water_budget[istart:,3],'-r',label='observed')
    pylab.plot_date(dts,z_surfs,'-b',label='modelled')
    pylab.plot_date(dts[-1:],z_surfs[-1:],'ob')
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().set_major_locator(locator)
    pylab.legend()
    pylab.grid(True)
    pylab.ylabel('surface level (m)')

    # Plot time series for salinity.
    ax = pylab.subplot2grid((2,2), (1, 1))
    pylab.cla()
    pylab.plot_date(dts,salts,'-b',label='modelled')
    pylab.plot_date(dts[-1:],salts[-1:],'ob')
    pylab.grid(True)
    pylab.xlim(water_budget_dt[istart],water_budget_dt[-1])
    ax.get_xaxis().set_major_locator(locator)
    pylab.ylim(100.,450.)
    pylab.ylabel('salinity (g/L)')

    pylab.savefig(path,dpi=dpi)
    print('Saved frame %s for %s.' % (path,strdate))

dts = []
Vs = []
z_surfs = []
areas = []
sources = []
salts = []
h_salts = []
iframe = 0
if animate:
    fig = pylab.figure(figsize=(10,5))
    pylab.subplots_adjust(bottom=0.12,left=0.15,wspace=0.25)
years = numpy.arange(startyear,2012,1)
startdate = datetime.datetime(startyear,1,1)
for iyear in years:
    for iday in range(365):
        # Compute salt layer thickness (m) per bottom elevation bin from salt density (kg/m2)
        h_salt = salt_per_level/dens_NaCl

        # Compute current bottom levels by taking original ("rock bottom") and adding salt layer.
        current_bottom_level = z_grid_cent+h_salt

        # compute updated cumulative area distribution as function of bottom level
        isort = numpy.argsort(current_bottom_level)      # sort updated bottom levels from low (deepest point) to high
        new_z_grid_cent = current_bottom_level[isort]    # new bottom levels in correct order
        new_area_dist = area_dist[isort]                 # new area per bottom level bin

        # Compute updated layer thicknesses (estimate from elevation at centre of bins)
        new_z_grid_delta = numpy.zeros_like(new_z_grid_cent)
        halfdh = (new_z_grid_cent[1:]-new_z_grid_cent[:-1])/2
        new_z_grid_delta[:-1] += halfdh
        new_z_grid_delta[1:] += halfdh
        new_z_grid_delta[0] += halfdh[0]
        new_z_grid_delta[-1] += halfdh[-1]

        # Find current surface level and lake surface area
        z_surf,wet_area_tot = level_and_area_from_volume(new_area_dist,new_z_grid_cent,new_z_grid_delta,V)

        # Find forcing index for current time
        curdate = datetime.datetime(iyear,1,1)+datetime.timedelta(days=iday)
        dt = pylab.date2num(curdate)
        idat = water_budget_dt.searchsorted(dt)

        # Water sinks/sources
        runoff = numpy.interp(dt,water_budget_dt,water_budget[:,0])*86400        # from m3/s to m3/d
        p = numpy.interp(dt,water_budget_dt,water_budget[:,1])/1000*wet_area_tot # from mm/m2/d to m3/d
        e = numpy.interp(dt,water_budget_dt,water_budget[:,2])/1000*wet_area_tot # from mm/m2/d to m3/d
        dV = runoff + p - e   # chnage in volume: sources minus sinks

        # If we reach max surface level (minus safety), overflow: any additional water is lost.
        # This happens occasionally with the (strange!) provded bathymetry, which has a max elevation
        # less than the max observed surface level.
        if z_surf>=max_z_surf-0.1: dV = min(0.,dV)

        # Compute salt precipitation
        salt = S/V  # salinity (kg/m3 = g/L) from total salt content (kg) and current volume (m3)
        Na_conc = salt/(Na_molweight+Cl_molweight)
        NaCl_sat = Na_conc**2/K_NaCl   # note this should be corrected for ion activity of Na+ and Cl-, e.g., Krumgalz & Millero 1989 Mar Chem
        salt_precip = NaCl_precip_rate*max(0.,NaCl_sat-1.)
        spec_salt_dissol = numpy.zeros_like(area_dist)
        spec_salt_dissol[salt_per_level>0] = NaCl_dissol_rate*max(0.,1.-NaCl_sat)
        dS = -salt_precip*V + (spec_salt_dissol*area_dist).sum()

        # Compute prepitated salt per bottom depth bin (multiply by total volume in bin)
        # Salt is assumed to precipitate throughout the water column, and to sink instantaneously to the bottom.
        salt_per_level += salt_precip*numpy.maximum(0.,z_surf-current_bottom_level) - spec_salt_dissol

        # Save frame for animation if desired.
        if animate and ((curdate-startdate).days%10==0):
            iframe += 1
            printStill(z_grid_cent,area_dist,h_salt,z_surf,dts,z_surfs,salts,'animation/%05i.png' % iframe,curdate)

        # Time integrate (Forward Euler, delta_t=1 day)
        V += dV
        S += dS

        # Store several variables for current time step for later plotting.
        dts.append(dt)
        Vs.append(V)
        areas.append(wet_area_tot)
        z_surfs.append(z_surf)
        sources.append((p,e,runoff))
        salts.append(salt)

    # Store salt layer thickness at end of current year.
    h_salts.append(h_salt)

    print('31 Dec %i, volume = %f km3, wet area = %f km2, level = %f m, salt = %i PPT' % (iyear,V/1e9,wet_area_tot/1e6,z_surf,S/V))

sources = numpy.array(sources)/86400
dts = numpy.array(dts)
areas = numpy.array(areas)
Vs = numpy.array(Vs)

pylab.figure()
pylab.plot_date(dts,z_surfs,'-',label='modelled')
pylab.plot_date(water_budget_dt,water_budget[:,3],'-',label='observed')
pylab.ylabel('surface elevation (m)')
pylab.xlim(dts[0],None)
pylab.ylim(1267.,None)
pylab.legend()
pylab.grid(True)
pylab.savefig('level.png',dpi=150)

# Plot time series of water sources and sinks
pylab.figure()
pylab.plot_date(dts,sources[:,0],'-',label='precipitation')
pylab.plot_date(dts,sources[:,1],'-',label='evaporation')
pylab.plot_date(dts,sources[:,2],'-',label='runoff')
pylab.ylabel('water flux (m3/s)')
pylab.grid(True)
pylab.legend()
pylab.savefig('water_sources.png',dpi=150)

# Plot time series of salinity
pylab.figure()
pylab.plot_date(dts,salts,'-')
pylab.ylabel('salinity (g/L)')
pylab.grid(True)
pylab.savefig('salinity.png',dpi=150)

with open('salinity.dat','w') as f:
    for i in range(len(dts)):
        f.write('%s\t%s\n' % (pylab.num2date(dts[i]).strftime('%y-%m-%d'),salts[i]))

# Plot precipitated salt as function of bottom elevation, for different years.
pylab.figure()
for iyear,h_salt in zip(years[::-1],h_salts[::-1],):
    pylab.plot(z_grid_cent,h_salt/salt_porosity,'-',label='%i' % iyear)
pylab.xlabel('bottom elevation (m)')
pylab.ylabel('precipitated salt (m)')
pylab.legend()
pylab.grid(True)
pylab.savefig('salt_layer_thickness.png',dpi=150)

# Plot lake surface area as function of lake level (for validation of hypsograph logic)
pylab.figure()
pylab.plot(z_surfs,areas/1e6,'o',label='model')
pylab.plot(data[:,0],data[:,1]/1e6,'-',label='input hypsograph')
pylab.grid(True)
pylab.xlabel('level (m)')
pylab.ylabel('lake surface area (km2)')
pylab.savefig('level_area.png',dpi=150)

if data.shape[1]>2:
    # Plot lake volume as function of lake level (for validation of hypsograph logic)
    pylab.figure()
    pylab.plot(z_surfs,Vs/1e9,'o',label='model')
    pylab.plot(data[:,0],data[:,2]/1e9,'-',label='input hypsograph')
    pylab.grid(True)
    pylab.xlabel('level (m)')
    pylab.ylabel('lake volume (km3)')
    pylab.savefig('level_volume.png',dpi=150)

