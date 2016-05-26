#!/usr/bin/env python

import pandas
import netCDF4 as nc
import argparse

"""
Read a list of edits which were provided by and apply them to a topo file.
"""

parser = argparse.ArgumentParser()
parser.add_argument('--orig',type=str,help='path to original netCDF file',default=None)
parser.add_argument('--edits',type=str,help='path to edits file',default=None)
parser.add_argument('--scale',type=float,help='scaling for input depths',default=1.0)
parser.add_argument('--output',type=str,help='path to output netCDF file',default='topog_edits.nc')


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


args=parser.parse_args()


if args.orig is None:

    print "Need to supply the original topography netCDF file"
    raise

if args.edits is None:

    print "Need to supply the path to edits "
    raise

f=nc.Dataset(args.orig)
depth=args.scale*f.variables['depth'][:]


try:
    Edits=pandas.read_csv(args.edits,header=None,names=['jEdit','iEdit','New'])

    jEdits=[]
    for a in Edits.jEdit:
        jEdits.append(int(a.replace('(','')))

    iEdits=[]
    for a in Edits.iEdit:
        iEdits.append(int(a))

    New=[]
    for a in Edits.New:
        New.append(float(a.replace(')','')))
        
    
except:
    try:
        fnam=args.edits+'.nc'
        i=nc.Dataset(fnam)
        iEdits=i.variables['iEdit'][:]
        jEdits=i.variables['jEdit'][:]
        New=i.variables['zEdit'][:]
        i.close()
    except:
        pass
    

Orig=[]
for i,j,new in zip(iEdits,jEdits,New):
    Orig.append(depth[j,i])
    print i,j,depth[j,i],new
    
g=nc.Dataset(args.output,'w',format='NETCDF3_CLASSIC')

dims=[]
for d in f.dimensions:
    print d
    dimvals=f.variables[d][:]
    dims.append((d,len(f.dimensions[d]),dimvals))


for d in dims:
    g.createDimension(d[0],d[1])
    
g.createDimension('nEdits',None)

dimv=[]
for d in dims:
    dimv.append(g.createVariable(d[0],'f8',(d[0])))
    dimv[-1].units='degrees'

for v,d in zip(dimv,dims):
    v[:]=d[2]
    
for v in f.variables:
    if v.find('longitude')==-1 and v.find('latitude')==-1 and v.find('count')==-1 and v.find('min')==-1 and v.find('max')==-1:
        varout=g.createVariable(v,'f4',('latitude','longitude'))
        try:
            units=f.variables[v].units
            varout.units=units
        except:
            pass
        try:
            standard_name=f.variables[v].standard_name
            varout.standard_name=standard_name
        except:
            pass
        try:
            description=f.variables[v].description
            varout.description=description
        except:
            pass
        try:
            long_name=f.variables[v].long_name
            varout.long_name=long_name
        except:
            pass

        dat=f.variables[v][:]

        if v.find('depth')>-1:
            for i,j,d in zip(iEdits,jEdits,New):
                print 'Modifying Elevation at j,i= ',j,i,' Old= ',dat[j,i],' New= ',d
                dat[j,i]=d

        varout[:]=dat

        
ivar=g.createVariable('iEdit','i4',('nEdits'))
jvar=g.createVariable('jEdit','i4',('nEdits'))
kvar=g.createVariable('zEdit','f4',('nEdits'))

n=0
for i,j,d in zip(iEdits,jEdits,Orig):
    ivar[n]=i
    jvar[n]=j
    kvar[n]=d
    n=n+1



g.sync()
g.close()
