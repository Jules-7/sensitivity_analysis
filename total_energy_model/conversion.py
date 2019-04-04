from math import sqrt
from general.constants import g, k, R, beta_trop, T0, p0, rho0, a0, ft, knot, fpm, mps_kts, mps_fpm, meter, nm

from general.atmosphere import Atmosphere
# ------ SPEED CONVERSIONS -----------------


def tas_to_cas(tas, altitude, delta_T=0.0, lapse_rate=beta_trop):
    """ BADA User Manual rev 3-12. (3.1-24), p.11
        BADA formulas differ from formulas in BlueSky

        v_tas in m/s
        height in m
    """
    local = Atmosphere(altitude, delta_T, lapse_rate)  # atmosphere properties for given height
    p = local.pressure
    rho = round(local.density, 3)
    mu = get_mu()

    v_cas = (((2./mu) * (p0/rho0) *
            (((1 + (p/p0) *
            (((1 + (mu/2.) * (rho/p) * (tas**2) )**(1./mu)) - 1)) ** mu) - 1)) ** (1./2.))  # [m/s]
    return v_cas


def cas_to_tas(cas, altitude, delta_T=0.0, lapse_rate=beta_trop):
    """ BADA User Manual rev 3-12. (3.1-23), p.11

        BADA formulas differ from formulas in BlueSky

        v_cas in m/s
        height in m
    """
    local = Atmosphere(altitude, delta_T, lapse_rate)  # atmosphere properties for given height
    p = local.pressure
    rho = local.density
    mu = get_mu()

    v_tas = (((2./mu) * (p/rho) *
            (((1 + (p0/p) *
            (((1 + (mu/2.) * (rho0/p0) * (cas**2) )**(1./mu)) - 1)) ** mu) - 1)) ** (1./2.))  # [m/s]
    return v_tas


def get_mu():
    """ BADA USer MAnual rev 3-12. (3.1-25), p.11 """
    mu = (k - 1) / k
    return mu


def mach_to_tas(mach, altitude, delta_T=0.0, lapse_rate=beta_trop):
    """ BADA User Manual rev 3-12. (3.1-26), p.11"""
    atm = Atmosphere(altitude, delta_T, lapse_rate)  # atmosphere properties for given height
    temp = atm.temp
    v_tas = mach * sqrt(k * R * temp)  # [m/s]
    return v_tas


def mach_to_cas(mach, altitude, delta_T=0.0, lapse_rate=beta_trop):
    tas = mach_to_tas(mach, altitude, delta_T, lapse_rate)
    cas = tas_to_cas(tas, altitude, delta_T, lapse_rate)
    return cas


def tas_to_mach(tas, altitude, delta_T=0.0, lapse_rate=beta_trop):
    local = Atmosphere(altitude, delta_T, lapse_rate)  # atmosphere properties for given height
    a = local.speed_of_sound
    mach = float(tas) / float(a)
    return mach

# --------- END SPEED CONVERSIONS -------------------


def stall_velocity(mass_act, phase, apm, configuration):
    """ aircraft mass changes the speed where stall will occur"""

    if phase == "CLIMB" and configuration == "TAKE-OFF":
        velocity_ref = apm.V_stall_to_ref

    elif phase == "CLIMB" and configuration == "INITIAL CLIMB":
        velocity_ref = apm.V_stall_ic_ref

    elif phase == "CLIMB" and configuration == "CLEAN":
        velocity_ref = apm.V_stall_cr_ref

    elif phase == "CRUISE":
        velocity_ref = apm.V_stall_cr_ref

    elif phase == "DESCENT" and configuration == "CLEAN":
        velocity_ref = apm.V_stall_cr_ref

    elif phase == "DESCENT" and configuration == "APPROACH":
        velocity_ref = apm.V_stall_ap_ref

    elif phase == "DESCENT" and configuration == "LANDING":
        velocity_ref = apm.V_stall_ld_ref

    velocity_act = velocity_ref * sqrt(mass_act / apm.mass_ref)

    return velocity_act

def crossover_altitude(cas, mach):
    """ Transition altitude (crossover altitude) in [m] between V_cas [m/s] and M,
        is the geopotential pressure altitude, at which V_cas and M represent the same V_tas [m/s]

        BADA User Manual rev 3-12. (3.1-27), (3.1-28), (3.1-29), p.12
    """
    sigma_trans = (
                  (((1 + ((k - 1)/2.) * ((cas/a0)**2))**(k/(k-1))) - 1)
                  /
                  (((1 + ((k-1)/2.)*(mach**2))**(k/(k-1))) - 1)
                  )

    teta_trans = sigma_trans**(-((beta_trop * R) / g))
    H_trans = (1000/6.5) * (T0*(1 - teta_trans))  # [m]
    return H_trans


def feet_to_meters(h):
    """ convert from feet to meters"""
    return h * ft


def meters_to_feet(h):
    """ convert from meters to feet """
    return h * meter


def knot_to_mps(speed_in_kts):
    """ convert speed from knots to m/s"""
    return speed_in_kts * knot


def fpm_to_mps(speed_in_fpm):
    """ convert speed from fpm to mps """
    return speed_in_fpm * fpm


def mps_to_kts(speed_in_mps):
    """ convert speed from mps to knots"""
    return speed_in_mps * mps_kts


def mps_to_fpm(speed_in_mps):
    """ convert speed from mps to fpm"""
    return speed_in_mps * mps_fpm


def nm_to_meters(dist_in_nm):
    """ convert from nautical miles to meters"""
    return dist_in_nm * nm



