import midas
import numpy as np
import netCDF4 as nc

path_seawifs='SeaWiFS_Non_Elnino_clim.nc'
grid=midas.rectgrid.quadmesh(path_seawifs,cyclic=True,var='chl')
sgrid = midas.rectgrid.supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid2=midas.rectgrid.quadmesh(supergrid=sgrid)
C=midas.rectgrid.state(path=path_seawifs,grid=grid,fields=['chl'])
grid2.D=nc.Dataset('topog.nc').variables['depth'][:]
grid2.wet[grid2.D==0.]=0
C2=C.horiz_interp('chl',target=grid2)
chl_tav=np.mean(C2.chl,axis=0)
chl_tav=chl_tav[np.newaxis,:,:,:]
C2.chl=np.ma.filled(C2.chl,chl_tav)
C2.chl=np.ma.filled(C2.chl,1.e-3)
C2.var_dict['chl']['xax_data']=grid2.x_T[grid2.jm/2,:]
C2.var_dict['chl']['yax_data']=grid2.y_T[:,grid2.im/4]
C2.var_dict['chl']['calendar']='julian'
path_out='seawifs.nc'
C2.write_nc(path_out,['chl'])
