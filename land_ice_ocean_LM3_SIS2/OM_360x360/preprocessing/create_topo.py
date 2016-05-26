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
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.basemap import interp 

def blend12(d1,d2, f1, f2):
  nj,ni=d1.shape
  x=np.arange(0,nj,dtype=float)/(nj-1)
  x=(x-f1)/(f2-f1)
  x=np.maximum(0.,x)
  x=np.minimum(1.,x)
  weight1d=1. - x
  print weight1d
  weight1d=weight1d[:,np.newaxis]
  weight=np.tile( weight1d, (1,ni) )
  return d1*weight+d2*(1-weight)



# Read model grid information from the supergrid file
sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
grid.lath=grid.y_T[:,grid.im/4]
grid.latq=grid.y_T_bounds[:,grid.im/4]
grid.lonh=grid.x_T[grid.jm/2,:]
grid.lonq=grid.x_T_bounds[grid.jm/2,:]
# Define a transition region between Bedmap2 and GEBCO
blend_region = grid.geo_region(y=(grid.latq[0],-70))
grid_blend_region =grid.extract(blend_region)
# Bedmap2 is on a polar stereographic grid. Use Basemap to handle the projection.
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
# Read grid information for displaced pole Southern Cap
scap=supergrid(file='sgrid0.nc')
scap_grid=quadmesh(supergrid=scap,is_latlon=True,cyclic=True)
xx=scap_grid.x_T_bounds.copy()
yy=scap_grid.y_T_bounds.copy()
xx[xx>180.]=xx[xx>180.]-360.
xx[xx<-180.]=xx[xx<-180.]+360.
x2,y2 = m(xx,yy,inverse=False)
cart_grid_so = supergrid(config='cartesian',axis_units='none',xdat=x2,ydat=y2)
# Remapping GEBCO to global grid
ingrid=quadmesh('GEBCO_08_v2.nc',var='depth',simple_grid=True,cyclic=True)
TOPO=state('GEBCO_08_v2.nc',grid=ingrid,fields=['depth'])
TOPO.rename_field('depth','topo')
TOPO.var_dict['topo']['Ztype']='Fixed'
fnam='interpolated_topog.nc'   
R=TOPO.subtile('topo',target=grid)
R.rename_field('mean','depth')
R.write_nc(fnam,['depth','max','min','std','count'])
# Remap GEBCO in Southern cap region
ingrid=quadmesh('GEBCO_08_v2.nc',var='depth',simple_grid=True,cyclic=True)
sp_reg=ingrid.geo_region(y=(-90.0,scap.y.max()+1.0))
TOPO=state('GEBCO_08_v2.nc',grid=ingrid,geo_region=sp_reg,fields=['depth'])
TOPO.rename_field('depth','topo')
TOPO.var_dict['topo']['Ztype']='Fixed'
fnam = 'scap_topog_gebco.nc'
R=TOPO.subtile('topo',target=scap_grid)
R.write_nc(fnam,['mean','max','min','std','count'])
# Read Bedmap2 data
TOPO=state('bedmap2.nc',grid=grid_bedmap,fields=['elev_bed','height_gl04c_wgs84','mask_ice','thick','elev_surf'])
TOPO.elev_bed = TOPO.elev_bed - TOPO.height_gl04c_wgs84
TOPO.elev_bed=np.ma.masked_where(np.isnan(TOPO.elev_bed),TOPO.elev_bed)
TOPO.elev_surf = TOPO.elev_surf - TOPO.height_gl04c_wgs84
TOPO.elev_surf=np.ma.masked_where(np.isnan(TOPO.elev_surf),TOPO.elev_surf)
TOPO.thick[TOPO.mask_ice==0.0]=0.0
TOPO.elev_bed[np.logical_and(TOPO.elev_surf>=0.0,TOPO.thick==0.)]=0.0
TOPO.rename_field('elev_bed','topo')
# Remap Bedmap2 topography to blended region 
xx=grid_blend_region.x_T_bounds.copy()
yy=grid_blend_region.y_T_bounds.copy()
xx[xx>180.]=xx[xx>180.]-360.
xx[xx<-180.]=xx[xx<-180.]+360.
x2,y2 = m(xx,yy,inverse=False)
cart_grid_blendregion = supergrid(config='cartesian',axis_units='none',xdat=x2,ydat=y2)
# Remap Bedmap2 ice sheet thickness in blending region
R=TOPO.subtile('topo',target=cart_grid_blendregion)
fnam='blend_region_topog_bedmap2.nc'
R.write_nc(fnam,['mean','max','min','std','count'])
# Remap Bedmap2 ice sheet thickness
R=TOPO.subtile('thick',target=cart_grid_blendregion)
R.rename_field('mean','thick')
R.thick[R.thick<1.0]=0.0
ice_area = np.zeros(R.thick.shape)
ice_area=grid_blend_region.Ah[np.newaxis,np.newaxis,:,:]
ice_area[R.thick==0.0]=0.0
vdict=R.var_dict['thick'].copy()
vdict['units']='m2'
R.add_field_from_array(ice_area,'area',var_dict=vdict)
fnam='blend_region_ice_shelf.nc'
R.write_nc(fnam,['thick','area'])   
