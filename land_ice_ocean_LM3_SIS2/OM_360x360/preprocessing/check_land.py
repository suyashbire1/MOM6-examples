import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='mill')
di=[1,1,0,-1,-1,-1,0,1]
dj=[0,-1,-1,-1,0,1,1,1]
 
def plot_river(bmap,P):
  XX=[];YY=[]
  for p in P[1:-1]:
    xx=x[p[0],p[1]]
    yy=y[p[0],p[1]]
    xx,yy=m(xx,yy)
    XX.append(xx);YY.append(yy)
  bmap.plot(XX,YY,color='b',alpha=0.2,linewidth=0.5)
  p=P[0]
  xx=x[p[0],p[1]]
  yy=y[p[0],p[1]]
  xx,yy=m(xx,yy)
  bmap.scatter(xx,yy,color='g',alpha=0.2,marker=',',s=1.0)
  p=P[-1]
  xx=x[p[0],p[1]]
  yy=y[p[0],p[1]]
  xx,yy=m(xx,yy)
  bmap.scatter(xx,yy,color='r',alpha=0.2,marker=',',s=1.0)        

for tile in np.arange(1,7):

  path='river_data.tile'+str(tile)+'.nc'
  tocell=nc.Dataset(path).variables['tocell'][:]
  x=nc.Dataset(path).variables['x'][:]
  y=nc.Dataset(path).variables['y'][:]
  land_frac=nc.Dataset(path).variables['land_frac'][:]
  lake_frac=nc.Dataset(path).variables['lake_frac'][:]
  ylen=x.shape[0]
  xlen=x.shape[1]

  LSINK=[]
  OSINK=[]
  
  for j in np.arange(1,ylen-1):
    for i in np.arange(1,xlen-1):
      start=(j,i)
      cur=start
      P=[]    
      if tocell[cur[0],cur[1]] > 0:
        next=np.int(np.log2(tocell[cur[0],cur[1]]))
        P.append(start)      
      else:
        next=-1

      while next>=0:
        if tocell[cur[0],cur[1]]>0:
          next=np.int(np.log2(tocell[cur[0],cur[1]]))
          lnext=(np.int(cur[0]+dj[next]),np.int(cur[1]+di[next]))
          if lnext[0]<0 or lnext[1]<0:# If the river extends into the next tile, then stop here. 
            next=-1
          elif lnext[0]>=ylen or lnext[1]>=xlen:# If the river extends into the next tile, then stop here. 
            next=-1
          else:
            cur=lnext
            P.append(cur)
        else:
          next=-1

      if len(P)>1:

        cur=P[-1]
        if cur[0]<ylen-1 and cur[1]<xlen-1 and cur[0]>0 and cur[1]>0:
          if land_frac[cur[0],cur[1]]>=1.0 :
            if lake_frac[cur[0],cur[1]]==0.0:
#              print 'River Ends at Full Land Cell which is not a lake: ',cur, ' tile= ',tile
#              print P
              plot_river(m,P)
              if cur not in LSINK:
                LSINK.append(cur)

            if land_frac[cur[0],cur[1]]==0.0:
#              print 'River Ends at Full Ocean Cell: ',cur , ' tile= ',tile
              plot_river(m,P)
              if cur not in OSINK:
                OSINK.append(cur)

  if len(LSINK)>0:
    print '********Found ',len(LSINK),' Land sink points on tile# ',tile,'********'
    print LSINK
  if len(OSINK)>0:
    print '********Found ',len(OSINK),' Ocean sink points on tile# ',tile,'********'    
    print OSINK
m.drawcoastlines()
plt.savefig('river_network.png')

