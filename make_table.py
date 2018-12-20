import soxs
import numpy as np
import h5py
from yt.units import clight

arf = soxs.AuxiliaryResponseFile.from_instrument("lynx_lxm_ultra")
rmf = soxs.RedistributionMatrixFile.from_instrument("lynx_lxm_ultra")

ne = arf.emid.size
agen = soxs.ApecGenerator(arf.elo[0], arf.ehi[-1], ne)

nv = 81
nT = 199
vel = np.linspace(-200.0, 200.0, nv) # km/s
kT = np.linspace(0.05, 5.0, nT) # keV

spec_table = np.zeros((nT,nv,ne))

ckms = clight.to_value("km/s")

for i, T in enumerate(kT):
    for j, v in enumerate(vel):
        spec = agen.get_spectrum(T, 1.0, v/ckms, 1.0)
        spec.apply_foreground_absorption(0.02)
        cspec = soxs.ConvolvedSpectrum(spec, arf)
        model = rmf.convolve_spectrum(cspec, (500.0, "ks"), noisy=False)
        spec_table[i,j,:] = model[:]

f = h5py.File("spec_table.h5", "w")
f.create_dataset("vel", data=vel)
f.create_dataset("kT", data=kT)
f.create_dataset("spec_table", data=spec_table)
f.flush()
f.close()

