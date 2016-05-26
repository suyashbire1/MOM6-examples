#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

def blend12(d1,d2, f1, f2):
  nj,ni=d1.shape
  print nj,ni
  x=np.arange(0,nj,dtype=float)/(nj-1)
  x=(x-f1)/(f2-f1)
  x=np.maximum(0.,x)
  x=np.minimum(1.,x)
  weight1d=1. - x
  print weight1d
  weight1d=weight1d[:,np.newaxis]
  weight=np.tile( weight1d, (1,ni) )
  return d1*weight+d2*(1-weight)

bedmap=nc.Dataset('scap_topog_bedmap2.nc')
gebco_scap=nc.Dataset('scap_topog_gebco.nc')
gebco=nc.Dataset('interpolated_topog.nc','a')

depth_scap=blend12(bedmap.variables['mean'][:],gebco_scap.variables['mean'][:],0,2)

