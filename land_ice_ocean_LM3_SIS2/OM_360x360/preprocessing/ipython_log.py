# IPython log file

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
f=nc.Dataset('edit_topog.nc')
depth=f.variables['depth'][:]
depth=np.ma.masked_where(depth<=0.,depth)
plt.pcolormesh(depth);plt.show()
