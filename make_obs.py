import yt
from yt.units import mp, kboltz
import pyxsim
import soxs
import numpy as np

ds = yt.load("cgols_simulation.gdf")

def _temperature(field, data):
    return (2./3.)*data["gdf", "thermal_energy"]*0.592*mp/kboltz/data["gas","density"]
ds.add_field(("gas", "temperature"), _temperature, units="K", force_override=True)

tm = pyxsim.ThermalSourceModel("apec", 0.05, 3.0, 50000, Zmet=1.0, kT_min=0.05, kT_max=5.0)

c = ds.domain_center
w = ds.quan(600.0, "pc").to('code_length')
le = c - 0.5*w
re = c + 0.5*w
box = ds.box(le, re)

redshift = 0.0
area = (2.5, "m**2")
exp_time = (500.0, "ks")
dist = (1.0, "Mpc")

photons = pyxsim.PhotonList.from_data_source(box, redshift, area, exp_time, tm, dist=dist)

photons.write_h5_file("cgols_photons.h5")

theta = 60

theta_rad = np.deg2rad(theta)
events = photons.project_photons([np.sin(theta_rad), 0.0, np.cos(theta_rad)], (30.0, 45.0), absorb_model='wabs', nH=0.02, north_vector=[0.0, 0.0, 1.0])

events.write_simput_file("cgols_%sdeg" % theta, overwrite=True)

soxs.instrument_simulator("cgols_%sdeg_simput.fits" % theta, "cgols_evt_%sdeg.fits" % theta, (500.0, "ks"), "lynx_lxm_ultra", (30.0, 45.0), overwrite=True, ptsrc_bkgnd=False)

