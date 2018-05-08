import yt
import numpy as np
from yt import derived_field

yt.enable_parallelism()

@derived_field(name="RhoChi", units="")
def _RhoChi(field, data) :
    value = data["rho"] * (data["chi"])**(-1.5)
    return value

loc = '/scratch2/kclough/ASBH/phi0.08M1.0r80plt000000*'
ds = yt.load(loc)

time_data = []
integral_data = []

for i in ds :
    field = 'RhoChi'
    weight = 'cell_volume'
    ad = i.all_data()
    ad2 = ad.cut_region('obj["chi"] > 0.1')

    #integral = ad2.quantities.weighted_average_quantity(field, weight)
    integral = ad2.quantities.total_quantity(field)

    time_data.append(i.current_time)
    integral_data.append(integral)

np.savetxt('time.out', time_data)
np.savetxt('integral.out', integral_data)
