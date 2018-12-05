import yt
import numpy as np
from yt import derived_field

yt.enable_parallelism()

@derived_field(name="HamSquared", units="")
def _HamSquared(field, data) :
    value = data["Ham"] * data["Ham"]
    return value

@derived_field(name="RhoChi", units="")
def _RhoChi(field, data) :
    value = data["rho"] * ( data["chi"]**(-1.5) )
    return value
filename = 'RhoChiOut_r60HR_1.dat'
loc = '/scratch2/kclough/ASBH/phi0.08M1.0r60HRplt*'
ds = yt.load(loc)

from mpi4py import MPI
comm = MPI.COMM_WORLD
		
#integral_data = []

volume = ( float(ds[0].domain_right_edge[0]) )**3.0

for i in ds :
    #field = 'HamSquared'
    field = 'RhoChi'
    weight = 'cell_volume'
    ad = i.all_data()
    #ad2 = ad.cut_region('obj["chi"] > 0.1')

    integral = ad.quantities.weighted_average_quantity(field, weight)
    #integral = ad.quantities.total_quantity(field)
    integral = integral * volume

    if(comm.rank==0) :
        datafile=open(filename,'a')
        datafile.write("%f    %.8f \n" % (i.current_time, integral))
        datafile.close()

#np.savetxt('RhoChi.out', integral_data)
