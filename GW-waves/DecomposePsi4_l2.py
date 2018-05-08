import yt
import numpy as np
from numpy import pi
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time 
yt.enable_parallelism()


# ===============================================
#  Decomposition of Psi4 into spinweighted spherical harmonics  
#               l = 2,3,4 
#   Input: Any dataset with Weyl1, Weyl2 and Radius of sphere r
#       Output: \psi_4 *r 
# ===============================================

# ===============================================
# Radius of sphere (center is set automatically in as domain_right_edge/2)
Rad = 120
# ===============================================

# Data file
loc = '/scratch2/kclough/ASBH/M1.0a0.0NRplt000*'

#Loading dataset

ds = yt.load(loc)

# =================
# Init Outputarrays
# =================
# l = 2  
Weyl4_l2_m0_data  = []
# positive m  
Weyl4_l2_m1_data = []
Weyl4_l2_m2_data = []
# negative m  
Weyl4_l2_m1n_data = []
Weyl4_l2_m2n_data = []

Sr_data = []
trianglearea = []
timedata = []
Cycletimedata = []

# Definitions for Quadrature sceme 

N = 131

coord = np.loadtxt('PointDistFiles/lebedev/lebedev_%03d.txt'% N)
theta = coord[:,1]*pi/180; phi = coord[:,0]*pi/180 + pi
w = coord[:,2]

# Spinweighted spherical harmonics s = -2 
# =========================================
# L = 2
# =========================================
# m = 0 zero 
sY_l2_m0 = np.sqrt(15./(32.*np.pi))*np.sin(theta)**2
# positive m values : 1,2
sY_l2_m1 = np.exp(1j*phi)*np.sqrt(5./(np.pi))*np.cos(theta/2)**3*np.sin(theta/2)
sY_l2_m2 = 1./2.*np.exp(2*1j*phi)*np.sqrt(5./(np.pi))*np.cos(theta/2)**4
# negative m values :-1,-2
sY_l2_m1n = 1./2.*np.exp(-1j*phi)*np.sqrt(5./(np.pi))*np.sin(theta/2)**2*np.sin(theta)
sY_l2_m2n = 1./2.*np.exp(-2*1j*phi)*np.sqrt(5./(np.pi))*np.sin(theta/2)**4
  
# ==============================
#         Loop over all frames 
# ==============================
for i in ds:

  startCycle = time.time()

  i.print_stats() 
  center =  (i.domain_right_edge/2)

  start = time.time()

# ==================================================
  # Initalising  

  print ("Integrating ... ")  
  
  Sr = 0 + 1j*0  
  # l = 2  
  Weyl4_l2_m0 = 0 + 1j*0
  # positive m  
  Weyl4_l2_m1 = 0 + 1j*0
  Weyl4_l2_m2 = 0 + 1j*0
  # negative m  
  Weyl4_l2_m1n = 0 + 1j*0 
  Weyl4_l2_m2n = 0 + 1j*0

  index = 0 
  for (k,x) in enumerate(phi):
    phi_var = phi[k]
    theta_var = theta[k]
    x1 = Rad*np.cos(phi_var)*np.sin(theta_var)+float(center[0])
    y1 = Rad*np.sin(phi_var)*np.sin(theta_var)+float(center[1])
    z1 = Rad*np.cos(theta_var)+float(center[2])

    c = [x1,y1,z1]
    ptn = i.point(c)
    ReWeyl = float(ptn["Weyl4_Re"][0])
    ImWeyl = float(ptn["Weyl4_Im"][0])
    Weyl4 = ReWeyl + 1j*ImWeyl

    Weyl4_l2_m0 += 4*pi*w[k]*np.conjugate(sY_l2_m0[k])*Weyl4*Rad
    # positive m  
    Weyl4_l2_m1 += 4*pi*w[k]*np.conjugate(sY_l2_m1[k])*Weyl4*Rad
    Weyl4_l2_m2 += 4*pi*w[k]*np.conjugate(sY_l2_m2[k])*Weyl4*Rad
    # negative m  
    Weyl4_l2_m1n += 4*pi*w[k]*np.conjugate(sY_l2_m1n[k])*Weyl4*Rad
    Weyl4_l2_m2n += 4*pi*w[k]*np.conjugate(sY_l2_m2n[k])*Weyl4*Rad

# ==================================================
# DATA WRITEOUT
# ==================================================

  # l = 2  
  Weyl4_l2_m0_data.append(Weyl4_l2_m0)
  # positive m  
  Weyl4_l2_m1_data.append(Weyl4_l2_m1)  
  Weyl4_l2_m2_data.append(Weyl4_l2_m2)  
  # negative m  
  Weyl4_l2_m1n_data.append(Weyl4_l2_m1n)  
  Weyl4_l2_m2n_data.append(Weyl4_l2_m2n)  

  # time
  timedata.append(i.current_time)

  np.savetxt('time.out',timedata)
  # l = 2 
  np.savetxt('Weyl4_l2_m0_data.out',Weyl4_l2_m0_data)
  np.savetxt('Weyl4_l2_m1_data.out',Weyl4_l2_m1_data)
  np.savetxt('Weyl4_l2_m2_data.out',Weyl4_l2_m2_data)
  np.savetxt('Weyl4_l2_m1n_data.out',Weyl4_l2_m1n_data)
  np.savetxt('Weyl4_l2_m2n_data.out',Weyl4_l2_m2n_data)
  mid1 = time.time()
  print("time to Integrate= ", mid1-start)

# ==================================================
# Performance 
# ==================================================
  Cycletime = abs(time.time()-startCycle)
  Cycletimedata.append(Cycletime)
  np.savetxt('Cycletime.out',Cycletimedata)
# ==================================================



print("FINISHED !!! ")


