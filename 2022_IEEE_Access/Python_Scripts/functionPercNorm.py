import numpy as np
import scipy.stats

def functionPercNorm(slice_original_quantile): # Wigert et al 2018
    vQuantiles = scipy.stats.mstats.mquantiles(slice_original_quantile, prob=[0.001, 0.999], alphap=0.5, betap=0.5)

    quantile999low = vQuantiles[0]
    print(vQuantiles)
    quantile999max = vQuantiles[1]

    slice_original_quantile = np.where(slice_original_quantile<quantile999low, quantile999low, slice_original_quantile)
    slice_original_quantile = np.where(slice_original_quantile>quantile999max, quantile999max, slice_original_quantile)

    if quantile999max > 10:
        slice_original_quantile = np.double(slice_original_quantile - quantile999low)/np.double(quantile999max - quantile999low)

    else:
        slice_original_quantile = np.zeros_like(slice_original_quantile)
    return slice_original_quantile
