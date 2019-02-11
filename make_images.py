import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve
import h5py

plt.rc("font", size=18)

basenm = "f4"
angle = "90"

fcts = fits.open("%s_cts_%sdeg.fits" % (basenm, angle))
fcube = fits.open("%s_cube_%sdeg.fits" % (basenm, angle))

cube = fcube[0].data
cts = fcts[0].data
cts = np.ma.masked_where(cts <= 0.0, cts)

fcts.close()
fcube.close()

f = h5py.File("spec_table.h5", "r")
vels = f["vel"][:]
f.flush()
f.close()

range1 = (0, np.searchsorted(vels, -70.0)-1)
range2 = (range1[1], np.searchsorted(vels, 70.0)-1)
range3 = (range2[1], vels.size)

R = cube[range3[0]:range3[1],:,:].sum(axis=0)
G = cube[range2[0]:range2[1],:,:].sum(axis=0)
B = cube[range1[0]:range1[1],:,:].sum(axis=0)
I = R+G+B
R /= I
G /= I
B /= I
R[np.isnan(R)] = 0.0
G[np.isnan(G)] = 0.0
B[np.isnan(B)] = 0.0

kernel = Gaussian2DKernel(x_stddev=1.0)

alpha = np.log10(cts)/np.log10(cts.max()*0.05)
#alpha = np.sqrt(cts)/np.sqrt(cts.max())
alpha[np.isnan(alpha)] = 0.0
np.clip(alpha, a_min=0.0, a_max=1.0, out=alpha)

Rd = convolve(R, kernel)
Bd = convolve(B, kernel)
Gd = convolve(G, kernel)
A = convolve(alpha, kernel)

rgb = np.transpose(np.array([Rd,Gd,Bd]), axes=(1,2,0)).astype("float32")
rgb *= alpha[:,:,np.newaxis]

cmap = plt.cm.viridis
cmap.set_bad(color="black")

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20,10))
ax1.imshow(cts, origin='lower', norm=LogNorm(), cmap=cmap)
ax2.imshow(rgb, origin='lower')
fig.savefig("%s_%sdeg.png" % (basenm, angle))

hdur = fits.ImageHDU(data=Rd*A)
hdur.name = "red"
hdur.writeto("%s_red_%sdeg.fits" % (basenm, angle), overwrite=True)
hdub = fits.ImageHDU(data=Bd*A)
hdub.name = "blue"
hdub.writeto("%s_blue_%sdeg.fits" % (basenm, angle), overwrite=True)
hdug = fits.ImageHDU(data=Gd*A)
hdug.name = "green"
hdug.writeto("%s_green_%sdeg.fits" % (basenm, angle), overwrite=True)
fits.ImageHDU(data=rgb).writeto("%s_rgb_%sdeg.fits" % (basenm, angle), overwrite=True)
