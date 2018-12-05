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

#Loading dataset
loc = '/scratch2/kclough/Inflation/Inf9HROldchk0*'
ds = yt.load(loc)

#initial_location
location =  np.array([0,0,0])
#location =  ds[0].domain_right_edge/2
variable = 'phi'
time_data = []
var_data = []
CycleData = []

for i in ds:
  print("Now at ", i.current_time)
  # find the center of the BH
  #value, location = i.find_min("chi")
  #center = [float(location[0]), float(location[1]), float(location[2])]
  #print ("New center ", center)
  point = i.point(location)
  var_out = point[variable]
  if(comm.rank==0) :
    datafile=open(filename,'a')
    datafile.write("%f    %.8f \n" % (i.current_time, var_out))
    datafile.close()
