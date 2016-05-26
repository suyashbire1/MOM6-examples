from midas.rectgrid import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc

sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
grid.x_T_bounds[-1,0:grid.im/2]=grid.x_T_bounds[-1,grid.im/2+1:]
grid.y_T_bounds[-1,0:grid.im/2]=grid.y_T_bounds[-1,grid.im/2+1:]
grid.lath=grid.y_T[:,grid.im/4]
grid.latq=grid.y_T_bounds[:,grid.im/4]
grid.lonh=grid.x_T[grid.jm/2,:]
grid.lonq=grid.x_T_bounds[grid.jm/2,:]




def show_regional_grid(x_bnds,y_bnds,fig=1,sec=None,vmin=None,vmax=None,cmap=plt.cm.spectral):
        fig=plt.figure(fig)
        ax1=fig.add_subplot(111)
#        ax2=fig.add_subplot(122)
        xind,yind,xindc,yindc,XAX,YAX,XAXc,YAXc,D,WET=extract_domain(x_bnds,y_bnds,grid)
        cf=ax1.pcolormesh(XAX,YAX,np.ma.masked_where(D<=0.0,D),vmin=vmin,vmax=vmax,cmap=cmap)
        ax1.contour(XAXc,YAXc,WET,[0.5,0.5],colors='k')
        if sec is not None:
            print sec
            xpos=[grid.x_T[sec[1],sec[0]],grid.x_T[sec[3],sec[2]]]
            ypos=[grid.y_T[sec[1],sec[0]],grid.y_T[sec[3],sec[2]]]
            print xpos
            print ypos
#        ax1.plot([xpos[0],ypos[0]],[xpos[1],ypos[1]],linestyle='dashed',linewidth=2,color='w')
            
        ax1.set_title('Depth')
        ax1.set_xlabel('Longitude')
        ax1.set_ylabel('Latitude')
        plt.colorbar(cf,ax=ax1)
#        cf=ax2.pcolormesh(xind,yind,np.ma.masked_where(D<0.0,D),vmin=vmin,vmax=vmax,cmap=cmap)
#        ax2.contour(xindc,yindc,WET,[0.5,0.5],colors='k')
#        plt.colorbar(cf,ax=ax2)
#        ax2.set_title('Depth')
#        if sec is not None:
#            ax2.plot([sec[0],sec[1]],[sec[2],sec[3]],linestyle='dashed',linewidth=2,color='w')
#        ax2.set_xlabel('I-index')
#        ax2.set_ylabel('J-index')


def extract_domain(x_bnds,y_bnds,grid):
    xind=sq(np.where(np.logical_and(grid.lonq>=x_bnds[0],grid.lonq<x_bnds[1])))
    yind=sq(np.where(np.logical_and(grid.latq>=y_bnds[0],grid.latq<y_bnds[1])))
    xindc=xind[:-1]
    yindc=yind[:-1]
    XAX=np.take(np.take(grid.x_T_bounds,xind,axis=1),yind,axis=0)
    YAX=np.take(np.take(grid.y_T_bounds,xind,axis=1),yind,axis=0)
    XAXc=np.take(np.take(grid.x_T,xindc,axis=1),yindc,axis=0)
    YAXc=np.take(np.take(grid.y_T,xindc,axis=1),yindc,axis=0)
    D=np.take(np.take(grid.D,xindc,axis=1),yindc,axis=0)
    WET=np.take(np.take(grid.wet,xindc,axis=1),yindc,axis=0)
    xindc=0.5*(xind[0:-1]+xind[1:])
    yindc=0.5*(yind[0:-1]+yind[1:])
    return xind,yind,xindc,yindc,XAX,YAX,XAXc,YAXc,D,WET

        
def shade_stereo(var,grid,hem=1,fig=1,field=None,tlev=0,klev=0,vmin=None,vmax=None,cmap=plt.cm.spectral,scale=1.0,nlevs=0):
    wd=6667000.
    ht=6667000.
    fig=plt.figure(fig)
    m =  Basemap(projection='stere',width=wd,height=ht,lon_0=0.0,lat_ts=hem*71.,lat_0=hem*90.,resolution='l')
    xx=grid.x_T_bounds.copy()
    yy=grid.y_T_bounds.copy()
    print xx.min(),xx.max()
    print yy.min(),yy.max()
    xx[xx>180.]=xx[xx>180.]-360.
    xx[xx<-180.]=xx[xx<-180.]+360.
    x2,y2 = m(xx,yy,inverse=False)
    zout=sq(var)*scale
    print zout.min(),zout.max()
    ax1=fig.add_subplot(111)
    cf=ax1.pcolormesh(x2,y2,sq(zout),vmin=vmin,vmax=vmax,cmap=cmap)
    plt.colorbar(cf)
    xx=grid.x_T.copy()
    yy=grid.y_T.copy()
    xx[xx>180.]=xx[xx>180.]-360.
    xx[xx<-180.]=xx[xx<-180.]+360.
    x2,y2 = m(xx,yy,inverse=False)
#    ax1.contour(x2,y2,grid.D,[0.,500.,1000,1500,2000],colors='w')
    m.drawcoastlines()
    if nlevs>0:
        ax1.contour(x2,y2,sq(zout),np.linspace(vmin,vmax,nlevs),colors='k')

south=grid.geo_region(y=(grid.latq[0],-40))
S=state(grid=grid,geo_region=south)
D=S.grid.dxh/S.grid.dyh
#shade_stereo(D,S.grid,fig=1,hem=-1,scale=1,vmin=0,vmax=6000,cmap=plt.cm.gist_stern)
shade_stereo(D,S.grid,fig=1,hem=-1,scale=1,vmin=0.02,vmax=1.02,cmap=plt.cm.spectral)
plt.show()
plt.show()

