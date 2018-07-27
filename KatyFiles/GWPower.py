import yt
import numpy as np
from numpy import pi
yt.enable_parallelism()

filename = 'GWEnergy.txt'
loc = '/gss/scratch/pr48pu/di36hup2/phi0.07/BHAS_ASNovaPlt_001*'
ds = yt.load(loc)
#write data
if __name__ == "__main__" :
  datafile=open(filename,'a')
  datafile.write("Time    Power     Energy  \n")

N=89
coord = np.loadtxt('/home/hpc/pr48pu/di36hup2/YTAnalysis/GW-waves/PointDistFiles/lebedev/lebedev_%03d.txt'% N)
theta = coord[:,1]*pi/180;
phi = coord[:,0]*pi/180 + pi
w = coord[:,2]

time0 = ds[0].current_time
time1 = ds[1].current_time
DeltaT = float(time1-time0)
#extraction radius
Rad = 300
Starttime = 0
Startindex = (Starttime)/DeltaT
counter = 1 
Int = np.zeros(len(coord),dtype = "complex128")

Energy=0
for i in ds :
  if ((counter > Startindex) and ((int(float(i.current_time)/1.6)-999)%4 == 0)) :

    center =  (i.domain_right_edge/2)

    #integrate over sphere
    Spher = 0 
    for (k,x) in enumerate(phi):

      phi_var = phi[k]
      theta_var = theta[k]
      x1 = Rad*np.cos(phi_var)*np.sin(theta_var)+float(center[0])
      y1 = Rad*np.sin(phi_var)*np.sin(theta_var)+float(center[1])
      z1 = Rad*np.cos(theta_var)+float(center[2])
      c = [x1,y1,z1]
      ptn = i.point(c)
      ReWeyl = float(ptn["Weyl4_Re"])
      ImWeyl = float(ptn["Weyl4_Im"])
      Weyl4 = ReWeyl + 1j*ImWeyl

      Int[k] += Weyl4*DeltaT
      Spher += 4*pi*w[k]*np.absolute(Int[k])**2

    #normalise value of power
    Power = Spher*(Rad)**2/(16*np.pi)
    Energy += Power*DeltaT

    #write data
    if __name__ == "__main__" :
      datafile=open(filename,'a')
      datafile.write("%f    %E     %E  \n" % (i.current_time, Power, Energy))

  if ((int(float(i.current_time)/1.6)-999) % 4) == 0 :
    counter += 1 

print("FINISHED !!! ")
