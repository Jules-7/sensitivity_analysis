from math import atan2, cos, sin, degrees, radians, asin, sqrt, tan

from general.constants import R_e, g, ft, meter
from total_energy_model.conversion import cas_to_tas, tas_to_mach, crossover_altitude, mps_to_kts, mps_to_fpm, \
     tas_to_cas, knot_to_mps, stall_velocity, fpm_to_mps, mach_to_tas, mach_to_cas
from general.main_functions import Flight


class DescentPrediction(Flight):

    def __init__(self, name, FPL, path_to_save_xls):

        super(DescentPrediction, self).__init__(name, FPL)

        self.route = self.fpl["route"][::-1]
        self.phase = "DESCENT"
        self.path_to_save_xls = path_to_save_xls

        # set positive VS in descent for backwards computation
        self.VS_descent = abs(self.VS_descent)

        self.descent_crossover_altitude = crossover_altitude(self.cas_descent, self.mach_descent)

        self.distance_flown = 0
        self.descent_deccel_alt = []

        # ------------------- FLAGS ------------------------------------------------------------------------------------
        self.TARGET_CRUISE_ALTITUDE_IS_REACHED = False

        self.CLIMB = False
        self.CRUISE = False
        self.DESCENT = True
        self.LANDED = False

        self.ROCD = self.VS_descent

        self.time_in_turn = 0

        self.set_speeds(speed=self.cas_descent)

        self.init_route_position_track()

        self.main()

        if self.show_map:
            self.display_map()

    def descent_phase(self):

        # ---------------- CHECK SPEED SETTINGS ----------------------------------------------------------------
        if self.alt < self.descent_crossover_altitude:

            if self.cas_descent > knot_to_mps(250):
                if self.FL < 100:
                    self.cas_target = knot_to_mps(250)
                else:
                    self.cas_target = self.cas_descent
            else:
                self.cas_target = self.cas_descent

            # check if target CAS has changed
            self.check_target_CAS()

            # CHECK FOR CAS ACCELERATION/DECELERATION
            self.check_TAS_for_acc_dec(below=True, descent_calc=True)

            if self.ACCELERATE_TAS or self.DECELERATE_TAS:
                # if aircraft needs to accelerate - accelerate TAS
                # and check if target CAS has been reached
                # update acceleration/deceleration
                self.update_long_acc_dec()
                self.Mach = tas_to_mach(self.TAS, self.alt, self.delta_T, self.lapse_rate)
                self.CAS = tas_to_cas(self.TAS, self.alt, self.delta_T, self.lapse_rate)

            else:
                # if no acceleration is required - maintain constant CAS
                # convert value of constant CAS to TAS and to Mach
                self.TAS = cas_to_tas(self.CAS, self.alt, self.delta_T, self.lapse_rate)
                self.Mach = tas_to_mach(self.TAS, self.alt, self.delta_T, self.lapse_rate)
            # --------------------------------------------------------------------------------------------------

        else:  # above crossover altitude
            self.Mach = self.mach_descent

            self.TAS = mach_to_tas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.CAS = mach_to_cas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            # -------------------------------------------------------------------------------------------------

        self.set_rocd()

    def main(self):
        while not self.TARGET_CRUISE_ALTITUDE_IS_REACHED:

            # ------------ UPDATE ATMOSPHERIC CONDITIONS ---------------------------------------------------------------
            self.update_atmospheric_conditions()

            # if not self.V_stall_ld:
            #     self.V_stall_ld = stall_velocity(mass_act=self.mass, phase="DESCENT", apm=self.apm,
            #                                      configuration="LANDING")
            self.descent_phase()

            # -------------- UPDATE ALTITUDE AND FL --------------------------------------------------------------------
            self.alt += self.ROCD * self.delta_time
            self.FL = int(round(self.alt * meter)) / 100  # convert altitude to FL

            if self.alt > self.target_cr_alt:
                self.TARGET_CRUISE_ALTITUDE_IS_REACHED = True

            # --------------- UPDATE FLIGHT PHASE ----------------------------------------------------------------------
            if self.show_map:
                self.update_phase()
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------- UPDATE BEARING -------------------------------------------------------------------
            if self.TURN:
                """ integrate track """
                self.track += self.turn_rate * self.delta_time
                self.time_in_turn += self.delta_time

            else:
                """ before calculating next lat/lon get update of a current track
                update is done based on previous coordinate and target coordinate"""
                self.track = self.bearing(self.cur_lat, self.cur_lon,
                                          self.route[self.next_wpt][0],
                                          self.route[self.next_wpt][1])

            # ---------- UPDATE HEADING AND GROUND SPEED ---------------------------------------------------------------
            self.get_gamma()
            """ since horizontal speed component is not the same as true air speed - use horizontal speed component for
            ground speed calculations
            if dont want to use - just comment it
            """
            self.HAS = self.TAS * cos(self.gamma)  # same as self.HAS = sqrt(self.TAS**2 - self.ROCD**2)
            self.wind_components()
            self.wind_vector_from_components()
            """ !!! flip wind direction to ensure the same influence of wind as in a flight"""
            self.wind_direction_deg = (self.wind_direction_deg + 180) % 360
            self.wind_direction = radians(self.wind_direction_deg)
            self.get_wind_correction_angle()
            self.update_heading()
            self.update_ground_speed()
            # ----------------------------------------------------------------------------------------------------------

            # ------------- PERFORM FLIGHT -----------------------------------------------------------------------------
            self.update_position()
            self.distance_flown += self.d_s

            self.x, self.y = self.convert_flight_to_cartesian()

            if self.show_map:
                self.path_points.append((self.deg_lat_cur, self.deg_lon_cur))

            # if self.flight_time_counter == 1 or self.flight_time_counter % 600 == 0:
            #     self.write_states_to_txt()
            # ----------------------------------------------------------------------------------------------------------

            # ---------------- PERFORM CHECKS --------------------------------------------------------------------------
            # ToDo:checks
            if self.FLIGHT_END:
                # if we are approaching the last wp in our route - check only for reaching required coordinates
                self.check_coordinates()
                self.check_altitude()
            elif not self.TURN:
                # if we are not performing the turn - check distance to turn start
                # on every step check if we approached the turning distance and need to start turning
                self.compute_turn_distance_and_rate()
                self.get_distance_to_wpt()
                self.check_turn_init_distance()
            elif self.TURN:
                # while we are performing turn - check bearing to decide if we have already reached required new heading
                # or we still need to keep turning
                self.check_track()
            # ----------------------------------------------------------------------------------------------------------

        return
