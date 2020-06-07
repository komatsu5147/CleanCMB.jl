"""
    tsz(νGHz; units = "cmb", Tcmb = 2.725)

Spectrum of the thermal Sunyaev-Zeldovich effect.

*Reference*: Equation (V) in Appendix of Zeldovich, Sunyaev, Astrophys. Space Sci. 4, 301 (1969)

# Arguments
- `νGHz::Real`: frequency in units of GHz.

# Optional keyword arguments
- `units::String="cmb"`: units of the spectrum. For Rayleigh-Jeans temperature (brightness temperature) units, `units = "rj"`. The default is the CMB units.
- `Tcmb::Real=2.725`: present-day temperature of the CMB in units of Kelvin.
"""
function tsz(νGHz::Real; units = "cmb", Tcmb = 2.725)
    hp = 6.62607015e-34
    kB = 1.380649e-23
    x = hp * νGHz * 1e9 / kB / Tcmb
    g = (exp(x) - 1)^2 / x^2 / exp(x)
    tsz = x * coth(x / 2) - 4
    if units == "rj"
        return tsz / g
    else
        return tsz
    end
end

"""
    dust1(νGHz; Td = 19.6, βd = 1.6, νd = 353, units = "cmb", Tcmb = 2.725)

Spectrum of 1-component modified black-body thermal dust emission. The output is normalized to unity at `νGHz = νd`.

# Arguments
- `νGHz::Real`: frequency in units of GHz.

# Optional keyword arguments
- `Td::Real=19.6`: dust temperature  in units of Kelvin.
- `βd::Real=1.6`: dust emissivity index.
- `νd::Real=353`: frequency at which the output is normalized to unity.
- `units::String="cmb"`: units of the spectrum. For Rayleigh-Jeans temperature (brightness temperature) units, `units = "rj"`. The default is the CMB units.
- `Tcmb::Real=2.725`: present-day temperature of the CMB in units of Kelvin.
"""
function dust1(
    νGHz::Real;
    Td = 19.6,
    βd = 1.6,
    νd = 353,
    units = "cmb",
    Tcmb = 2.725,
)
    hp = 6.62607015e-34
    kB = 1.380649e-23
    x(ν) = hp * ν * 1e9 / kB / Tcmb
    g(ν) = (exp(x(ν)) - 1)^2 / x(ν)^2 / exp(x(ν))
    xd(ν) = hp * ν * 1e9 / kB / Td
    mbb_rj = (νGHz / νd)^(βd + 1) * (exp(xd(νd)) - 1) / (exp(xd(νGHz)) - 1)
    if units == "rj"
        return mbb_rj
    else
        return g(νGHz) / g(νd) * mbb_rj
    end
end

"""
    synch(νGHz; βs = -3, νs = 23, Cs = 0, νC = 40, units = "cmb", Tcmb = 2.725)

Spectrum of synchrotron emission. The output is normalized to unity at `νGHz = νs`.

# Arguments
- `νGHz::Real`: frequency in units of GHz.

# Optional keyword arguments
- `βs::Real=-3`: synchrotron power-law index.
- `νs::Real=23`: frequency at which the output is normalized to unity.
- `Cs::Real=0`: curvature of the synchrotron spectrum.
- `νC::Real=40`: pivot frequency for curvature of the synchrotron spectrum.
- `units::String="cmb"`: units of the spectrum. For Rayleigh-Jeans temperature (brightness temperature) units, `units = "rj"`. The default is the CMB units.
- `Tcmb::Real=2.725`: present-day temperature of the CMB in units of Kelvin.
"""
function synch(
    νGHz::Real;
    βs = -3,
    νs = 23,
    Cs = 0,
    νC = 40,
    units = "cmb",
    Tcmb = 2.725,
)
    hp = 6.62607015e-34
    kB = 1.380649e-23
    x(ν) = hp * ν * 1e9 / kB / Tcmb
    g(ν) = (exp(x(ν)) - 1)^2 / x(ν)^2 / exp(x(ν))
    synch_rj = (νGHz / νs)^(βs + Cs * log(νGHz / νC))
    if units == "rj"
        return synch_rj
    else
        return g(νGHz) / g(νs) * synch_rj
    end
end
