import rioxarray as rxr
from sklearn.cluster import AffinityPropagation
from sklearn.preprocessing import StandardScaler
import numpy as np

# data
eofs = rxr.open_rasterio('./global_4_EOFs_coswts_RESCALED_EPSG4326.tif')

# subset to CA, to practice
CA = eofs[:3].sel(x=slice(-122, -118), y=slice(42,32))
CA_x = CA.values
CA_xx = CA_x.reshape(CA_x.shape[0], -1).T

# tracks NaNs and non-NaNs
nans = np.where(np.isnan(CA_xx).sum(axis=1)>0)[0]
notnans = np.where(np.invert(np.isnan(CA_xx)).sum(axis=1)>0)[0]
CA_xxx = CA_xx[notnans,:]

# need to standardize?
# X = StandardScaler().fit_transform(CA_x)

#af = AffinityPropagation(preference=-50, random_state=0).fit(X)
af = AffinityPropagation().fit(CA_xxx)

