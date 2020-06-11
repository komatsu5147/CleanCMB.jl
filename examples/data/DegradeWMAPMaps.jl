using Healpix
using Printf
using Statistics
using FITSIO
topdir = "/Users/eiichirokomatsu/Documents/myJulia/wmap/data"
# %% Target Nside
nside_out = 128
npix = 12 * nside_out^2
# %% ILC region definition map
# The original file location: https://lambda.gsfc.nasa.gov/product/map/dr5/ilc_map_get.cfm
infile = joinpath(topdir, "wmap_ilc_rgn_defn_9yr_v5.fits")
bit_in, rgn_in = (readMapFromFITS(infile, ic, Float32) for ic = 1:2)
nside_in = bit_in.resolution.nside
fact = Int((nside_in / nside_out)^2)
bit = Map{Float32,RingOrder}(nside_out)
rgn = Map{Float32,RingOrder}(nside_out)
for ip = 1:npix
    degrade = (ip-1)*fact+1:ip*fact
    jp = nest2ring(bit.resolution, ip)
    bit[jp] = maximum(bit_in[degrade])
    rgn[jp] = maximum(rgn_in[degrade])
end
outfile = "wmap_ilc_rgn_defn_r7_9yr_v5.fits"
extname = "Archive Map Table"
typechar = "E"
f = FITSIO.fits_create_file(outfile)
FITSIO.fits_create_binary_tbl(
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
FITSIO.fits_close_file(f)
# %% WMAP temperature maps (smootehd to 1 degree Gaussian beam)
# The original file location: https://lambda.gsfc.nasa.gov/product/map/dr5/maps_band_smth_r9_i_9yr_get.cfm
band = ["K", "Ka", "Q", "V", "W"]
for iν = 1:5
    infile = joinpath(
        topdir,
        @sprintf("wmap_band_smth_imap_r9_9yr_%s_v5.fits", band[iν]),
    )
    map_in = readMapFromFITS(infile, 1, Float32)
    nside_in = map_in.resolution.nside
    fact = Int((nside_in / nside_out)^2)
    map = Map{Float32,RingOrder}(nside_out)
    for ip = 1:npix
        degrade = (ip-1)*fact+1:ip*fact
        jp = nest2ring(map.resolution, ip)
        map[jp] = mean(map_in[degrade])
    end
    outfile = @sprintf("wmap_band_smth_imap_r7_9yr_%s_v5.fits", band[iν])
    extname = "Archive Map Table"
    typechar = "E"
    f = FITSIO.fits_create_file(outfile)
    FITSIO.fits_create_binary_tbl(
        f,
        0,
        [("TEMPERATURE", "1$typechar", "mK, thermodynamic")],
        extname,
    )
    saveToFITS(map, f, 1)
    FITSIO.fits_close_file(f)
end
