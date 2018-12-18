using FITSIO
using HDF5
using StatsBase
using LinearAlgebra
using ProgressMeter

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

mat = -1
rmf_header = nothing
rmf = FITS(rmffile)
for i in 1:length(rmf)
    global mat
    global rmf_header
    rmf_header = read_header(rmf[i])
    if haskey(rmf_header, "EXTNAME")
        if occursin("MATRIX", rmf_header["EXTNAME"])
            mat = i 
            break
        end
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

cidxs = (chan .> cmin) .& (chan .< cmax)

x = read(hdu, "X")[cidxs]
y = read(hdu, "Y")[cidxs]
chan = chan[cidxs]

close(evtf)

# Load the model table in

th = h5open("spec_table.h5", "r")
kT = read(th, "kT")
vel = read(th, "vel")
# Must filter the spectrum part on the channels as above
m = read(th, "spec_table")[cmin:cmax,:,:]
close(th)

A = sum(m, dims=1)
logm = log.(m ./ A)

nT = length(kT)
nv = length(vel)
nc = size(m)[1]

cube = zeros((nx,ny,nv))
maxv = zeros((nx,ny))
cts = zeros((nx,ny))

for i in 1:nx
    p = Progress(nx, 1)
    xidxs = (x .>= xbins[i]) .& (x .<= xbins[i+1])
    for j in 1:ny
        pidxs = xidxs .& (y .>= ybins[j]) .& (y .<= ybins[j+1])
        if sum(pidxs) != 0
            n = counts(chan[pidxs], cmin:cmax)
            cts[i,j] = sum(n)
            lnL = sum(n .* logm, dims=1)[1,:,:]
            Lv = maximum(lnL, dims=2)[:,1]
            cube[i,j,:] = Lv
            maxv[i,j] = maximum(Lv)
        end
    end
    next!(p)
end

angle = split(evtfile, ".")[1][end-3:end]

velfile = "cgols_vel_$angle.fits"
cubefile = "cgols_cube_$angle.fits"
ctsfile =  "cgols_cts_$angle.fits"

cubef = FITS(cubefile, "w")
write(cubef, cube)
close(cubef)

velf = FITS(velfile, "w")
write(velf, maxv)
close(velf)

ctsf = FITS(ctsfile, "w")
write(ctsf, cts)
close(ctsf)

