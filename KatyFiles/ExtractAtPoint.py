import yt
from yt import derived_field
import numpy as np
from numpy.linalg import inv
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
yt.enable_parallelism()

# ==============================================================================
#      Spherical Horizon finder
#
#  Assumes (Obviously) spherical symmetry !! Gives bad results if not !
#  Input = Black hole center
#
#  For valid calculation of Mass one has to slice along x-Axis !!
#       Output :-Black Hole mass (with Error estimate)
#               -Black Hole radius (with Error estimate)
# ==============================================================================

#Loading dataset
loc = '/scratch2/kclough/ASBH/M0.0a0.0f0.06plt0*'
ds = yt.load(loc)

#initial_location
location =  ds[0].domain_right_edge/2
variable = 'phi'
time_data = []
var_data = []
CycleData = []

for i in ds:
  start = time.time()
  print(i.current_time)
  # find the center of the BH
  #value, location = i.find_min("chi")
  #center = [float(location[0]), float(location[1]), float(location[2])]
  #print ("New center ", center)
  point = i.point(location)
  var_data.append(point[variable]) 
  time_data.append(i.current_time)
  CycleData.append(time.time()-start)

np.savetxt('Cycletime.out', CycleData)
np.savetxt('var.out', var_data)
np.savetxt('time.out', time_data)
