import yt
import numpy as np
from yt import derived_field

yt.enable_parallelism()

@derived_field(name="HamSquared", units="")
def _HamSquared(field, data) :
    value = data["Ham"] * data["Ham"]
    return value

loc = '/scratch2/kclough/Inflation/Inf9HROldchk0*'
ds = yt.load(loc)

time_data = []
integral_data = []

for i in ds :
    field = 'HamSquared'
    weight = 'cell_volume'
    ad = i.all_data()
    #ad2 = ad.cut_region('obj["chi"] > 0.1')

    integral = ad.quantities.weighted_average_quantity(field, weight)
    #integral = ad2.quantities.total_quantity(field)

    integral = np.sqrt(integral)

    time_data.append(i.current_time)
    integral_data.append(integral)

np.savetxt('HamTime_old.out', time_data)
np.savetxt('AvHam_old.out', integral_data)
