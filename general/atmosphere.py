"""
   BADA  user manual rev 3-12
   EARTH ATMOSPHERE MODEL

   Air temperature, density and pressure is calculated for 2 layers:
   - Troposphere
   - Lower Stratosphere
   since aircraft don`t fly in upper Stratosphere

   Standard atmosphere at MSL
       occur at conditions of ISA - T0, p0, rho0, a0 and Hp (geopotential pressure altitude) == 0

   Non-Isa atmospheres
       provides pressure, temperature and density as a function of geopotential altitude and 2 differentials,
       follow same hypothesis as ISA atmosphere, but either delta_T and/or delta_p != 0
       delta_T - temperature differential at MSL
       delta_p - pressure differential at MSL
"""

from math import exp, sqrt
from general.constants import g, k, R, beta_trop, K, TROPOPAUSE, T0, p0, rho0, a0


class Atmosphere(object):

    """ ISA and non-ISA
        based on BADA User Manual rev 3-12.

        :output values for specified altitude:
        temp           - Air temperature [K]
        pressure       - Air pressure [Pa]
        density        - Air density [kg/m**3]
        speed of sound - Speed of sound [m/s]
    """

    def __init__(self, h, delta_t=0, lapse_rate=beta_trop):
        self.height = h  # [m]
        self.delta_t = delta_t  # [K] difference in temperature between ISA and actual temp -> non-ISA, i.e. ISA+20
        self.lapse_rate = lapse_rate
        self.temp = 0  # [K]
        self.pressure = 0
        self.density = 0
        self.speed_of_sound = 0

        self.trop_temp, self.trop_pressure, self.trop_density, self.trop_speed_of_sound = self.get_tropopause_values()

        self.get_temperature()
        self.get_pressure()
        self.get_density()
        self.get_speed_of_sound()

    def get_temperature(self):
        """ BADA User Manual rev 3-12. (3.1-13), (3.1-14), (3.1-15), (3.1-16), p.10"""
        if self.height <= TROPOPAUSE:
            self.temp = T0 + self.delta_t + self.lapse_rate * self.height
        elif self.height > TROPOPAUSE:  # if above tropopause - temp == temp at tropopause
            self.temp = self.trop_temp

    def get_pressure(self):
        """ BADA User Manual rev 3-12. (3.1-18), (3.1-19), (3.1-20), p.10"""
        if self.height <= TROPOPAUSE:
            self.pressure = p0 * ((float(self.temp - self.delta_t)/T0)**(-g/(self.lapse_rate*R)))
        elif self.height > TROPOPAUSE:
            # print "in here"
            pr = p0 * ((float(self.temp - self.delta_t)/T0)**(-g/(self.lapse_rate*R)))
            self.pressure = self.trop_pressure * exp(- (g/(R * self.temp)) * (self.height - TROPOPAUSE))

    def get_density(self):
        """ BASA User Manual rev 3-12. (3.1-21), p.10"""
        self.density = float(self.pressure) / (R * self.temp)
        density_v2 = rho0 * ((float(self.temp)/T0)**(-g/(self.lapse_rate*R) - 1))
        density_above_tropopause = self.trop_density * exp(- (g/(R * self.temp)) * (self.height - TROPOPAUSE))
        a = 1

    def get_speed_of_sound(self):
        """ BADA User Manual rev 3-12. (3.1-22), p.11"""
        self.speed_of_sound = sqrt(k * R * self.temp)

    def get_tropopause_values(self):
        """ BADA User Manual rev 3-12. (3.1-15), (3.1.-19), (3.1-21) p.10,
           (3.1-22) p.11 """
        trop_temp = T0 + self.delta_t + self.lapse_rate * TROPOPAUSE
        trop_pressure = p0 * ((float(trop_temp - self.delta_t)/T0)**(-g/(self.lapse_rate*R)))
        trop_density = float(trop_pressure) / (R * trop_temp)
        trop_speed_of_sound = sqrt(k * R * trop_temp)
        return trop_temp, trop_pressure, trop_density, trop_speed_of_sound


if __name__ == "__main__":
    atm = Atmosphere(12000., 0)
    print("temp in K", atm.temp, "in C", atm.temp - 273.15)
    print("pressure in Pa", atm.pressure)
    print("density", atm.density)
    print("speed of sound", atm.speed_of_sound)
    print("---------------------------------------------")
    # print "tropopause"
    # print "temp", atm.trop_temp, "in C", atm.trop_temp - 273.15
    # print "pressure", atm.trop_pressure
    # print "density", atm.trop_density
    # print "speed of sound", atm.trop_speed_of_sound

