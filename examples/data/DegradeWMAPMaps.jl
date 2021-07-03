using Healpix
using Printf
using Statistics
using CFITSIO
topdir = "/Users/eiichirokomatsu/Documents/work/myJulia/wmap/data"
# %% Target Nside
nside_out = 128
npix = 12 * nside_out^2
# %% ILC region definition map
# The original file location: https://lambda.gsfc.nasa.gov/product/map/dr5/ilc_map_get.cfm
infile = joinpath(topdir, "wmap_ilc_rgn_defn_9yr_v5.fits")
bit_in, rgn_in = (readMapFromFITS(infile, ic, Float32) for ic = 1:2)
nside_in = bit_in.resolution.nside
fact = Int((nside_in / nside_out)^2)
bit = HealpixMap{Float32,RingOrder}(nside_out)
rgn = HealpixMap{Float32,RingOrder}(nside_out)
for ip = 1:npix
    degrade = (ip-1)*fact+1:ip*fact
    jp = nest2ring(bit.resolution, ip)
    bit[jp] = maximum(bit_in[degrade])
    rgn[jp] = maximum(rgn_in[degrade])
end
outfile = "wmap_ilc_rgn_defn_r7_9yr_v5.fits"
extname = "Archive Map Table"
typechar = "E"
f = CFITSIO.fits_create_file(outfile)
CFITSIO.fits_create_binary_tbl(
    f,
    0,
    [
        ("TEMPERATURE", "1$typechar", "counts"),
        ("N_OBS", "1$typechar", "counts"),
    ],
    extname,
)
saveToFITS(bit, f, 1)
saveToFITS(rgn, f, 2)
CFITSIO.fits_close_file(f)
# %% WMAP temperature maps (smootehd to 1 degree Gaussian beam)
# The original file location: https://lambda.gsfc.nasa.gov/product/map/dr5/maps_band_smth_r9_i_9yr_get.cfm
band = ["K", "Ka", "Q", "V", "W"]
for iν = 1:5
    local infile = joinpath(
        topdir,
        @sprintf("wmap_band_smth_imap_r9_9yr_%s_v5.fits", band[iν]),
    )
    map_in = readMapFromFITS(infile, 1, Float32)
    local nside_in = map_in.resolution.nside
    local fact = Int((nside_in / nside_out)^2)
    map = HealpixMap{Float32,RingOrder}(nside_out)
    for ip = 1:npix
        degrade = (ip-1)*fact+1:ip*fact
        jp = nest2ring(map.resolution, ip)
        map[jp] = mean(map_in[degrade])
    end
    local outfile = @sprintf("wmap_band_smth_imap_r7_9yr_%s_v5.fits", band[iν])
    local extname = "Archive Map Table"
    local typechar = "E"
    local f = CFITSIO.fits_create_file(outfile)
    CFITSIO.fits_create_binary_tbl(
        f,
        0,
        [("TEMPERATURE", "1$typechar", "mK, thermodynamic")],
        extname,
    )
    saveToFITS(map, f, 1)
    CFITSIO.fits_close_file(f)
end
