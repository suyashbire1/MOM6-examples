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
import numpy as np
import netCDF4 as nc


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

def ice9it(i, j, depth, minD=0.):
  """
  Recursive implementation of "ice 9".
  Returns 1 where depth>minD and is connected to depth[j,i], 0 otherwise.
  """
  wetMask = 0*depth      

  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or depth[j,i] <= minD: continue
    wetMask[j,i] = 1
  
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) ) # Periodic beyond i=0

    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (j,0) ) # Periodic beyond i=ni-1
    
    if j>0: stack.add((j-1,i))
        
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) ) # Tri-polar fold beyond j=nj-1

  return wetMask

def ice9(x, y, depth, xy0):
  ji = nearestJI(x, y, xy0)
  return ice9it(ji[1], ji[0], depth)

def nearestJI(x, y, (x0, y0)):
  """
  Find (j,i) of cell with center nearest to (x0,y0).
  """
  return numpy.unravel_index( ((x-x0)**2 + (y-y0)**2).argmin() , x.shape)


# Read model grid information from the supergrid file
sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
grid.lath=grid.y_T[:,grid.im/4]
grid.latq=grid.y_T_bounds[:,grid.im/4]
grid.lonh=grid.x_T[grid.jm/2,:]
grid.lonq=grid.x_T_bounds[grid.jm/2,:]


# Now blend these datasets to create a single gridded field
bedmap_blend_region=nc.Dataset('blend_region_topog_bedmap2.nc')
gebco=nc.Dataset('interpolated_topog.nc','a')
jblend=bedmap_blend_region.variables['mean'].shape[0]
depth_new=blend12(bedmap_blend_region.variables['mean'][:],gebco.variables['depth'][0:jblend,:],0.2,1)
std_new=blend12(bedmap_blend_region.variables['std'][:],gebco.variables['std'][0:jblend,:],0.2,1)
depth=gebco.variables['depth'][:]
depth[0:jblend,:]=depth_new
std=gebco.variables['std'][:]
std[0:jblend,:]=std_new
blend_region_ish=nc.Dataset('blend_region_ice_shelf.nc')
thick=np.zeros(depth.shape)
area_ice=np.zeros(depth.shape)
thick[0:jblend,:]=blend_region_ish.variables['thick'][:]
area_ice[0:jblend,:]=blend_region_ish.variables['area'][:]
dens_ice=917.
dens_ocn=1035.
min_ice_cavity_clearance=50.
draft_ice = thick*dens_ice/dens_ocn
tmp=depth.copy()
tmp[draft_ice>-depth]=-(draft_ice[draft_ice>-depth]+min_ice_cavity_clearance)
depth[draft_ice>0.]=tmp[draft_ice>0.]
TOPO=state('interpolated_topog.nc',grid=grid,fields=['depth'])
vdict=TOPO.var_dict['depth'].copy()
TOPO.depth = -depth[np.newaxis,np.newaxis,:]
TOPO.add_field_from_array(thick[np.newaxis,np.newaxis,:],'ice_shelf_thickness',var_dict=vdict)
TOPO.add_field_from_array(area_ice[np.newaxis,np.newaxis,:],'ice_shelf_area',var_dict=vdict)
fnam='output_topog.nc'
#Remove head of Ross Shelf
#TOPO.depth[0,0,2,115:160]=0.0
#George V
#TOPO.depth[0,0,40:46,232]=200.0
#TOPO.ice_shelf_thickness[0,0,40:46,232]=min_ice_cavity_clearance
#Cleanup E Antarctica
#TOPO.depth[0,0,38:41,296]=0.0
#Remove Amery
#TOPO.depth[0,0,30:41,6:13]=0.0
#Open Gibraltar
#TOPO.depth[0,0,221,294]=300.
#Open Baltic
#TOPO.depth[0,0,249,311]=100.
#Open Black Sea
#TOPO.depth[0,0,226:228,326]=50.
#Open Red Sea
#TOPO.depth[0,0,195,342:344]=50.
#Indonesia
#TOPO.depth[0,0,166:170,41:43]=0.0
# ICE-9
#tmp = ice9(grid.x_T,grid.y_T,sq(TOPO.depth),(-140.,0.))
#TOPO.depth = tmp[np.newaxis,np.newaxis,:]
#Make sure ice shelf mask is consistent
TOPO.ice_shelf_thickness[TOPO.depth<=0.0]=0.0
TOPO.ice_shelf_area[TOPO.depth<=0.0]=0.0
TOPO.write_nc(fnam,['depth','ice_shelf_thickness','ice_shelf_area'])

