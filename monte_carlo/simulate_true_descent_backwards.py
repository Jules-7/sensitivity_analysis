from math import atan2, cos, sin, degrees, radians, asin, sqrt, tan, pi

from general.constants import R_e, g, ft, meter
from total_energy_model.conversion import cas_to_tas, tas_to_mach, crossover_altitude, mps_to_kts, mps_to_fpm, \
     tas_to_cas, knot_to_mps, stall_velocity, fpm_to_mps, mach_to_tas, mach_to_cas
from general.main_functions import Flight


class TrueDescentSimulation(Flight):

    def __init__(self, name, FPL, path_to_save_xls,
                 descent_distance, alt_start_decel):

        super(TrueDescentSimulation, self).__init__(name, FPL,
                                                    path_to_save_xls=path_to_save_xls,
                                                    descent_distance=descent_distance,
                                                    alt_start_decel=alt_start_decel)

        self.route = self.fpl["route"][::-1]
        self.phase = "DESCENT"

        self.path_to_save_xls = path_to_save_xls

        if self.fpl["instruction"]["type"] == "vs":
            # add 10 FLs to cruise altitude
            self.target_cr_alt = self.target_cr_alt + 10 * 100 * ft

        # set positive VS in descent for backwards computation
        self.VS_descent = abs(self.VS_descent)

        self.leveloffs = self.fpl["leveloff"]
        self.leveloff_descent_fl = self.leveloffs["DESCENT"][0][::-1]
        self.leveloff_descent_dur = self.leveloffs["DESCENT"][1][::-1]
        self.LEVELOFF = False
        self.leveloff_time = 0

        self.descent_crossover_altitude = crossover_altitude(self.cas_descent, self.mach_descent)
        self.descent_deccel_alt = []

        # distance flown in descent for TOD determination in flight
        self.distance_flown = 0

        self.TARGET_CRUISE_ALTITUDE_IS_REACHED = False

        self.DESCENT, self.CLIMB, self.CRUISE = True, False, False

        self.path_to_save_xls = path_to_save_xls

        if not self.test:
            self.create_txt_file()

        self.ROCD = self.VS_descent

        self.set_speeds(speed=self.cas_descent)

        self.init_route_position_track()

        # if self.show_map:
        #     self.init_pygmap()

        self.main()

        if not self.test:
            self.save_data()

        if self.show_map:
            self.display_map()

    def create_txt_file(self):
        """writing results into txt file is much faster than writing to excel file"""
        self.results_file = open(self.path_to_save_xls+self.name+"_"+str(self.route_id)+".txt", 'w')
        header = ",".join(["TIME_s", "FL", "ALT_m", "CAS_kt",
                     "M", "TAS_kt",
                     "GS_kt", "WS_kt", "ROCD_fpm",
                     "LAT_DD", "LONG_DD", "X_nm", "Y_nm",
                     "TRK_deg", "HDG_deg",
                     "WD_deg", "PHASE",
                     "TURN", "BANK_deg",
                     "DIST_nm", "EVENT",
                     "CUR_wpt", "NEXT_wpt",
                     "DIST_to_dest", "\n"])
        self.results_file.write(header)

        return

    def main(self):

        while not self.TARGET_CRUISE_ALTITUDE_IS_REACHED:

            # ------------ UPDATE ATMOSPHERIC CONDITIONS ------------
            self.update_atmospheric_conditions()

            self.descent_phase()

            if self.LEVELOFF:  # if level-off
                self.ROCD = 0.

            # -------------- UPDATE ALTITUDE AND FL -----------------
            self.alt += self.ROCD * self.delta_time

            self.FL = int(round(self.alt * meter) / 100)  # convert altitude to FL

            if not self.LEVELOFF:
                self.check_leveloff_constraints()

            if self.alt > self.target_cr_alt:
                self.TARGET_CRUISE_ALTITUDE_IS_REACHED = True
                # print("time\t", self.flight_time)
                # print("distance\t", self.distance_flown)

            if self.show_map:
                self.update_phase()

            # ----------------------- UPDATE BEARING ----------------
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

            # ---------- UPDATE HEADING, WIND AND GROUND SPEED ------
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
            # -------------------------------------------------------

            # ------------- PERFORM FLIGHT --------------------------
            self.update_position()
            self.distance_flown += self.d_s

            self.x, self.y = self.convert_flight_to_cartesian()

            # -------------------------------------------------------
            if self.show_map:
                self.path_points.append((self.deg_lat_cur, self.deg_lon_cur))

            if not self.test:
                # self.write_states_to_xls()
                self.write_states_to_txt()

            if self.LEVELOFF:
                self.update_leveloff_time()

            # ---------------- PERFORM CHECKS -----------------------
            if self.FLIGHT_END:
                # if we are approaching the last wp in our route - check only for reaching required coordinates
                self.check_coordinates()
                self.check_altitude()
            elif not self.TURN:
                # if we are not performing the turn - check distance to turn start
                # on every step check if we approached the turning distance and need to start turning
                self.time_in_turn = 0
                self.distance_in_turn = 0
                self.radial_dist = 0
                self.compute_turn_distance_and_rate()
                self.get_distance_to_wpt()
                self.check_turn_init_distance()
            elif self.TURN:
                # while we are performing turn - check bearing to decide if we have already reached required new heading
                # or we still need to keep turning
                self.check_track()

                # self.radial_dist = abs(self.turn_rate)
                # self.distance_in_turn = (degrees(self.radial_dist) / 360) * (2 * pi * self.r)

                # self.get_distance_to_wpt()
                # self.check_waypoint_switch()
            # -------------------------------------------------------

        return

    def check_leveloff_constraints(self):

        if self.phase == "DESCENT":
            if self.leveloff_descent_fl:
                leveloff_fl = self.leveloff_descent_fl[0]
                if leveloff_fl - 0.1 <= self.FL <= leveloff_fl + 0.1:
                    self.LEVELOFF = True

    def update_leveloff_time(self):
        self.leveloff_time += self.delta_time

        if self.DESCENT and self.leveloff_descent_fl:
            if self.leveloff_time >= self.leveloff_descent_dur[0]:
                # self.temp_phase = ""
                self.LEVELOFF = False
                del self.leveloff_descent_fl[0]
                del self.leveloff_descent_dur[0]
                self.leveloff_time = 0

    def descent_phase(self):
        if self.alt < self.descent_crossover_altitude:

            if self.cas_descent > knot_to_mps(250):
                if self.FL < 100:
                    self.cas_target = knot_to_mps(250)
                else:
                    self.cas_target = self.cas_descent
            else:
                self.cas_target = self.cas_descent

            self.target_spd = mps_to_kts(self.cas_target)

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
            # -------------------------------------------------------

        else:  # above crossover altitude
            self.Mach = self.mach_descent

            self.TAS = mach_to_tas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.CAS = mach_to_cas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.target_spd = self.Mach
            # -------------------------------------

        self.set_rocd()

    def write_states_to_txt(self):

        data_to_record = [self.flight_time, self.FL, self.alt, mps_to_kts(self.CAS),
                         self.Mach, mps_to_kts(self.TAS),
                         mps_to_kts(self.GS), mps_to_kts(self.wind_speed),
                         mps_to_fpm(self.ROCD),
                         self.deg_lat_cur, self.deg_lon_cur,
                         self.x, self.y,
                         degrees(self.track)%360,
                         degrees(self.hdg)%360,
                         self.wind_direction_deg, self.phase,
                         self.TURN, degrees(self.bank_angle),
                         self.distance_traveled / 1852., self.event,
                         self.cur_wpt, self.next_wpt,
                         self.dist_to_destination / 1852.]

        data = ['%.3f'%each if not isinstance(each, str) else each for each in data_to_record]

        data.append('\n')

        data_row = ",".join(data)

        self.results_file.write(data_row)
