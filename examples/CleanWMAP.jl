using CleanCMB
using Healpix
using Printf
using Statistics
using Plots
# %% Number of frequency bands
nν = 5
# %% Read in maps at Nside=128. The original map resolution is Nside=512, but we use degraded maps for faster demo.
# ILC region definition map
bit, rgn = (
    readMapFromFITS("data/wmap_ilc_rgn_defn_r7_9yr_v5.fits", ic, Float64) for
    ic = 1:2
)
nside = bit.resolution.nside
# WMAP temperature maps (smootehd to 1 degree Gaussian beam)
map = []
band = ["K", "Ka", "Q", "V", "W"]
for iν = 1:nν
    filename = @sprintf("data/wmap_band_smth_imap_r7_9yr_%s_v5.fits", band[iν])
    push!(map, readMapFromFITS(filename, 1, Float64))
end
# %% Perform ILC to obtain a clean CMB from WMAP temperature maps
Tbar = zeros(nν)
cij = zeros(nν, nν)
clean_map = Map{Float64,RingOrder}(nside)
for N = 0:11 # loop over ILC regions 0-11
    println("Region: ", N)
    set = (Int.(round.(bit)) .& 2^N) / 2^N # see https://lambda.gsfc.nasa.gov/product/map/dr5/ilc_map_get.cfm
    ip = findall(x -> x == 1, set)
    # Compute nν-by-nν covariance marix that will be used by `ilc_weights(cij)`
    for iν = 1:nν
        Tbar[iν] = mean(map[iν][ip])
    end
    for iν = 1:nν, jν = iν:nν
        cij[iν, jν] =
            mean((map[iν][ip] .- Tbar[iν]) .* (map[jν][ip] .- Tbar[jν]))
    end
    weights = ilc_weights(cij) # Compute ILC weights
    @show weights
    # Obtain a clean CMB map
    ip = findall(x -> x == N, rgn)
    for iν = 1:nν
        clean_map[ip] += weights[iν] * map[iν][ip]
    end
end
# %% Plot a clean CMB map and save to "ilc.png"
p = plot(clean_map, clims = (-0.3, 0.3))
savefig("ilc.png")
display(p)
