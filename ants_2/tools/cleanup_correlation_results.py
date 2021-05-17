from glob import glob
import h5py
import os
import numpy as np

if not os.path.exists("data/correlations/cleaned"):
    os.mkdir("data/correlations/cleaned")

corrs_in = glob("data/correlations/*.h5")

for cin in corrs_in:
    print("cleaning up: ", cin)
    f = h5py.File(cin, "r")
    fnew = os.path.join("data/correlations/cleaned", os.path.basename(cin))
    fnew = h5py.File(fnew, "a")

    sumcorrs = np.sum(f["corr_windows"]["data"][:], axis=-1)
    ixgood = np.where(sumcorrs != 0.0)



    # Save some basic information
    stats = fnew.create_dataset('stats', data=(0,))
    stats.attrs['sampling_rate'] = f["stats"].attrs['sampling_rate'] 
    stats.attrs['channel1'] = f["stats"].attrs['channel1']
    stats.attrs['channel2'] = f["stats"].attrs['channel2']
    stats.attrs['distance'] = f["stats"].attrs['distance']



    corr_windows = fnew.create_group("corr_windows")
    # create an array for the results
    data = corr_windows.create_dataset("data", data=f["corr_windows"]["data"][ixgood])
    data_keys = corr_windows.create_dataset("timestamps",
        data=f["corr_windows"]["timestamps"][ixgood])

    fnew.close()
    f.close()
