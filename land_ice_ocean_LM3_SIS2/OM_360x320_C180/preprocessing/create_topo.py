#!python

#============================================================
# Generate tiles for the northern/southern caps
# and central mercator grid.
#
# python create_topo.py
# Output: mercator_supergrid.nc, ncap_supergrid.nc, scap_supergrid.nc
# These are supergrids (2x grid tracer refinement) containing positions
# cell lengths, areas and angles
#
# Generate topography for grid tiles using BEDMAP for the Antarctic cap
# GEBCO 2 minute data for the Mercator grid and either
# IBCAO or GEBCO for the Northern cap (these files need to be linked to the
# current directory prior to running this command)

# python create_topo.py --tile ncap|scap|mercator 
#
#============================================================



from midas.rectgrid import *
from midas.rectgrid_gen import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.basemap import interp 

sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)

ingrid=quadmesh('GEBCO_08_v2.nc',var='depth',simple_grid=True,cyclic=True)
TOPO=state('GEBCO_08_v2.nc',grid=ingrid,fields=['depth'])
TOPO.rename_field('depth','topo')
TOPO.var_dict['topo']['Ztype']='Fixed'

fnam='interpolated_topog.nc'   
R=TOPO.subtile('topo',target=grid)
R.rename_field('mean','depth')
R.write_nc(fnam,['depth','max','min','std','count'])


wd=6667000.0
ht=6667000.0
m = Basemap(projection='stere',width=wd,height=ht,lon_0=0.0,lat_ts=-71.,lat_0=-90.,resolution='l')
f=netCDF4.Dataset('bedmap2.nc')
x1=sq(f.variables['x'][:])*1000 + 3333000.0
y1=x1
nx1=len(x1)
ny1=len(y1)
x1,y1=np.meshgrid(x1,y1)
grid_bedmap = quadmesh(lon=x1,lat=y1,is_latlon=False,is_cartesian=True,simple_grid=True)
TOPO=state('bedmap2.nc',grid=grid_bedmap,fields=['elev_bed','height_gl04c_wgs84'])
TOPO.elev_bed = TOPO.elev_bed - TOPO.height_gl04c_wgs84
TOPO.elev_bed=np.ma.masked_where(np.isnan(TOPO.elev_bed),TOPO.elev_bed)
TOPO.rename_field('elev_bed','topo')
fnam='scap_topog_bedmap2.nc'
scap=supergrid(file='sgrid0.nc')
scap_grid=quadmesh(supergrid=scap,is_latlon=True,cyclic=True)

xx=scap_grid.x_T_bounds.copy()
yy=scap_grid.y_T_bounds.copy()
xx[xx>180.]=xx[xx>180.]-360.
xx[xx<-180.]=xx[xx<-180.]+360.
x2,y2 = m(xx,yy,inverse=False)
cart_grid_so = supergrid(config='cartesian',axis_units='none',xdat=x2,ydat=y2)
R=TOPO.subtile('topo',target=cart_grid_so)
R.write_nc(fnam,['mean','max','min','std','count'])
