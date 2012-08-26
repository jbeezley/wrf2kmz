from netCDF4 import Dataset
import numpy as np

def relative_humidity(t2,q2,psfc):
    qrs=teten(t2,psfc)
    rh=relhum(q2,qrs)
    return rh


def teten(t,p):
    es=610.78*np.exp(17.269* (t-273.161)/(t-35.861))
    return 0.622*es/(p-0.378*es)

def relhum(qmixr,satmixr):
    return 100. * qmixr/satmixr
