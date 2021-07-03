import pysm3
import pysm3.units as u
import healpy as hp
import numpy as np
sky = pysm3.Sky(nside=256, preset_strings=["a1", "d1", "f1", "s1"], output_unit="uK_CMB")
nu=[27, 39, 93, 145, 220, 225, 280, 350, 410, 850]
for i in range (0,10):
    map = sky.get_emission(nu[i] * u.GHz)
    map_equ = pysm3.apply_smoothing_and_coord_transform(map, rot=hp.Rotator(coord="GC"))
    # RING format
    hp.write_map("map_equ_{:03d}ghz_r8_uKcmb.fits".format(nu[i]), map_equ, coord="C", dtype=np.float32)
    # NESTED format
    map_nested = hp.reorder(map_equ, r2n=True)
    hp.write_map("map_equ_{:03d}ghz_r8_nested_uKcmb.fits".format(nu[i]), map_nested, coord="C", nest=True, dtype=np.float32)
