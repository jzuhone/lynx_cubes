using FITSIO
using HDF5
using StatsBase
using LinearAlgebra

reblock = 1

emin = 0.3
emax = 1.0

evtfile = "cgols_evt_0deg.fits"

# Open event file 

evtf = FITS(evtfile)

hdu = evtf["EVENTS"]

header = read_header(hdu)

# Determine bounds of image

xmin = header["TLMIN2"]
ymin = header["TLMIN3"]
xmax = header["TLMAX2"]
ymax = header["TLMAX3"]

nx = div(Int(xmax-xmin), reblock)
ny = div(Int(ymax-ymin), reblock)

xbins = collect(range(xmin, stop=xmax, length=nx+1))
ybins = collect(range(ymin, stop=ymax, length=ny+1))

# Get correspondence between energy and channels from RMF

rmffile = header["RESPFILE"]

rmf = FITS(rmffile)
for i in 1:length(rmf)
    rmf_header = read_header(rmf[i])
    if occursin("MATRIX", rmf_header["EXTNAME"])
        mat = i 
        break
    end
end
elo = read(rmf[mat], "ENERG_LO")
ehi = read(rmf[mat], "ENERG_HI")
cbins = collect(range(rmf_header["TLMIN4"], stop=rmf_header["TLMAX4"]))
close(rmf)

cmin = cbins[elo .> emin][1]
cmax = cbins[ehi .< emax][end]

# Now filter out events that do not fit within these channels

chan = read(hdu, header["CHANTYPE"])

cidxs = (chan .> cmin) & (chan .< cmax)

x = read(hdu, "X")[cidxs]
y = read(hdu, "Y")[cidxs]
chan = chan[cidxs]

close(evtf)

# Load the model table in

th = h5open("spec_table.h5", "r")
kT = read(th, "kT")
vel = read(th, "vel")
# Must filter the spectrum part on the channels as above
m = read(th, "spec_table")[[cmin:cmax,:,:]]
close(th)

A = sum(spect, dims=1)
logm = log.(m ./ A)

nT = length(kT)
nv = length(vel)
nc = size(spect)[1]

cube = zeros((nx,ny,nv))

for i in 1:nx
    pidxs = (x >= xbins[i]) & (x <= xbins[i+1])
    for j in 1:ny
        pidxs = pidxs & (y >= ybins[j]) & (y <= ybins[j+1])
        if sum(pidxs) != 0
            n = counts(chan[pidxs], cmin:cmax)
            lnL = sum(n .* logm, dims=1)[1,:,:]
            Lv = maximum(lnL, dims=2)[:,1]
            # This next step is likely incorrect
            cube[i,j,:] = Lv
        end
    end
end
