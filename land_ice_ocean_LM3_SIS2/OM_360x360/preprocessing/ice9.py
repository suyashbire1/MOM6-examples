#!/usr/bin/env python

def error(msg,code=9):
  print 'Error: ' + msg
  exit(code)


# Imports
try: import argparse
except: error('This version of python is not new enough. python 2.7 or newer is required.')
try: from netCDF4 import Dataset
except: error('Unable to import netCDF4 module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_netcdf4')
try: import numpy as np
except: error('Unable to import numpy module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_numpy')
import shutil as sh

#import matplotlib.pyplot as plt

def main():

  # Command line arguments
  parser = argparse.ArgumentParser(description=
       'Applies an "ice 9" algorithm to remove detached water from the topography. Also sets land elevation to 0.',
       epilog='Written by A.Adcroft, 2013.')
  parser.add_argument('filename', type=str,
                      help='netcdf file to read.')
  parser.add_argument('variable', type=str,
                      nargs='?', default='depth',
                      help='Name of variable to plot.')
  parser.add_argument('--output', type=str,
                      nargs='?', default=' ',
                      help='name of the output file. If not specified, "iced_" is prepended to the name of the input file.')
  parser.add_argument('--shallow', type=float,
                      help='The "shallow" value (+ve, default 1.) to use when calculating the modified_mask. Wet points shallower than this are indicated with mask value of 2.',\
                      default=1.0)
  parser.add_argument('--analyze', action='store_true',
                      help='Report on impact of round shallow values to zero')

  optCmdLineArgs = parser.parse_args()

  nFileName = optCmdLineArgs.output
  if nFileName == ' ': nFileName = 'iced_'+optCmdLineArgs.filename
  shallow = optCmdLineArgs.shallow
  applyIce9(optCmdLineArgs.filename, nFileName, optCmdLineArgs.variable,
            0., -40., shallow, optCmdLineArgs.analyze)

def applyIce9(fileName, nFileName, variable, x0, y0, shallow, analyze):

  iRg = Dataset( fileName, 'r' );
  iDepth = iRg.variables[variable] # handle to the variable
  depth = iDepth[:] # Read the data
  print 'Range of input depths: min=',np.amin(depth),'max=',np.amax(depth)

    
  # Open new netcdf file
  if fileName==nFileName: error('Output file must be different from the input file')
  try: rg=Dataset( nFileName, 'w', format='NETCDF3_CLASSIC' );
  except: error('There was a problem opening "'+nFileName+'".')

  (ny, nx) = depth.shape
  rg.createDimension('nx',nx)
  rg.createDimension('ny',ny)
  rgDepth = rg.createVariable('depth','f4',('ny','nx'))
  rgDeptho = rg.createVariable('deptho','f4',('ny','nx'))  
  rgDepth.units = iDepth.units
  rgDepth.standard_name = 'depth below geoid'
  rgDepth.description = 'Non-negative nominal thickness of the ocean at cell centers (excluding ice shelf cavities)'
  rgDeptho.units = iDepth.units
  rgDeptho.standard_name = 'depth below geoid'
  rgDeptho.description = 'Non-negative nominal thickness of the ocean at cell centers (including ice shelf cavities)'
  rg.createDimension('ntiles',1)

  # A mask based solely on value of depth
  notLand = np.where( depth>0, 1, 0)

  lons=iRg.variables['longitude'][:]
  lats=iRg.variables['latitude'][:]
  
  iind=np.where(lons>-40.)[0][0]
  jind=np.where(lats>30.)[0][0]
  print jind,iind
  print lons[iind],lats[jind]
  print depth[jind,iind]
  notLand = ice9it(iind,jind,depth)

  rgWet = rg.createVariable('wet','f4',('ny','nx'))
  rgWet.long_name = 'Wet/dry mask'
  rgWet.description = 'Values: 1=Ocean, 0=Land'
  rgWet[:] = notLand # + (1-notLand)*0.3*np.where( depth<0, 1, 0)



  if 'std' in iRg.variables: # 
    rgH2 = rg.createVariable('h2','f4',('ny','nx'))
    rgH2.units = iDepth.units+'^2'
    rgH2.standard_name = 'Variance of sub-grid scale topography'
    rgH2[:] = iRg.variables['std'][:]**2

  ice_shelf_mask = np.zeros(depth.shape)
    
  if 'ice_shelf_thickness' in iRg.variables:
    rgIST = rg.createVariable('ice_shelf_thickness','f4',('ny','nx'))
    rgIST.units = 'm'
    rgIST.standard_name = 'Floating ice shelf thickness'
    rgIST[:] = iRg.variables['ice_shelf_thickness'][:]
    ice_shelf_mask[rgIST[:]>0.]=1.0
    
  if 'ice_shelf_area' in iRg.variables:
    rgISA = rg.createVariable('ice_shelf_area','f4',('ny','nx'))
    rgISA.units = 'm^2'
    rgISA.standard_name = 'Floating ice shelf area'
    rgISA[:] = iRg.variables['ice_shelf_area'][:]
    
  rgDepth[:] = depth*(notLand)*(1.0-ice_shelf_mask)
  rgDeptho[:] = depth*(notLand)
  
  
  if 'zEdit' in iRg.variables: # Need to copy over list of edits
    rgMod = rg.createVariable('modified_mask','f4',('ny','nx'))
    rgMod.long_name = 'Modified mask'
    rgMod.description = 'Values: 1=Ocean, 0=Land, -1 indicates water points removed by "Ice 9" algorithm. 2 indicates wet points that are shallower than 1m deep.'
    rgMod[:] = notLand - (1-notLand)*np.where( depth<0, 1, 0) + np.where( (notLand>0) & (depth>-shallow), 1, 0)
    n = len(iRg.variables['zEdit'])
    nEd = rg.createDimension('nEdits',n)
    iEd = rg.createVariable('iEdit','i4',('nEdits',))
    iEd.long_name = 'i-index of edited data'
    jEd = rg.createVariable('jEdit','i4',('nEdits',))
    jEd.long_name = 'j-index of edited data'
    zEd = rg.createVariable('zEdit','f4',('nEdits',))
    zEd.long_name = 'Original value of height data'
    zEd.units = iDepth.units
    iEd[:] = iRg.variables['iEdit'][:]
    jEd[:] = iRg.variables['jEdit'][:]
    zEd[:] = iRg.variables['zEdit'][:]

  rg.close()
  print 'File "%s" written.'%(nFileName)

  # Analyze the shallow points
  if analyze:
    print 'Analyzing...'
    numNotLand = np.count_nonzero(notLand)
    print '# of wet points after Ice 9 = %i'%(numNotLand)
    newDepth = depth*np.where(depth*notLand <= shallow, 1, 0)
    numNewWet = np.count_nonzero(newDepth)
    print '# of wet points deeper than %f = %i'%(shallow,numNewWet)
    print '%i - %i = %i fewer points left'%(numNotLand,numNewWet,numNotLand-numNewWet)
    newWet = ice9it(600,270,newDepth)
    numNewDeep = np.count_nonzero(newWet)
    print '# of wet deep points after Ice 9 = %i'%(numNewDeep)
    print '%i - %i = %i fewer points left'%(numNewWet,numNewDeep,numNewWet-numNewDeep)
  

def ice9it(i,j,depth):
  # Iterative implementation of "ice 9"
  wetMask = 0*depth
  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or depth[j,i] <= 0: continue
    wetMask[j,i] = 1
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) )
    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (0,j) )
    if j>0: stack.add( (j-1,i) )
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) )
  return wetMask

# Invoke main()
if __name__ == '__main__': main()

