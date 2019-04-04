from math import cos, degrees, radians, sqrt, acos
import numpy as np
from scipy import spatial

from general.constants import meter
from total_energy_model.conversion import crossover_altitude, mps_to_kts, mps_to_fpm, mach_to_tas
from general.main_functions import Flight


class FlightPrediction(Flight):

    def __init__(self, FPL, descent_distance, alt_start_decel,
                 real_lat, real_lon, real_alt, real_fl, real_phase):

        super(FlightPrediction, self).__init__('prediction', FPL)

        self.descent_distance = descent_distance

        self.delta_time = 0.2  # choose out of [0.1, 0.2, 0.25, 0.5, 1]
        self.multp_factor = int(1/self.delta_time)

        self.positions_at_lat = {5: [],
                                 10: [],
                                 15: [],
                                 20: []}

        self.u_array, self.v_array = [], []

        if alt_start_decel:
            self.alt_start_decel = alt_start_decel[0]
        else:
            self.alt_start_decel = 0

        self.climb_crossover_altitude = crossover_altitude(self.cas_climb, self.mach_climb)
        self.descent_crossover_altitude = crossover_altitude(self.cas_descent, self.mach_descent)

        self.dist_to_reach_desc_mach = self.get_distance_to_deceler_to_descent_mach()

        self.ROCD = self.VS_climb

        # data from the real flight
        self.phase = real_phase
        self.alt = real_alt
        self.FL = real_fl

        self.real_lat = real_lat
        self.real_lon = real_lon

        self.cur_wpt, self.next_wpt = None, None

        self.set_speeds(speed=self.cas_climb)

        self.init_success = False

        self.init_route_position_track()

    def init_route_position_track(self):

        self.get_bank_angle()  # [rad]

        self.cur_lat = radians(self.real_lat)
        self.cur_lon = radians(self.real_lon)
        self.lat_ref_rad, self.lon_ref_rad = self.route[0][0], self.route[0][1]
        self.deg_lat_cur = self.real_lat
        self.deg_lon_cur = self.real_lon

        # look at the route and based on the current position decide which
        # waypoint is current and which is next

        for i, wpt in enumerate(self.route[:-1]):
            if self.cur_lat == wpt[0] and self.cur_lon == wpt[1]:
                self.cur_wpt = i
                self.next_wpt = self.cur_wpt + 1
                self.init_success = True
                break

            # convert waypoints and position to cartesian
            wpt_1_x, wpt_1_y = self.convert_positions_to_cartesian(wpt[0], wpt[1])  # A
            wpt_2_x, wpt_2_y = self.convert_positions_to_cartesian(self.route[i + 1][0], self.route[i + 1][1])  # B
            position_x, position_y = self.convert_positions_to_cartesian(self.cur_lat, self.cur_lon)  # C

            # compute distances between points
            a = sqrt((position_x - wpt_2_x) ** 2 + (position_y - wpt_2_y) ** 2)
            b = sqrt((position_x - wpt_1_x) ** 2 + (position_y - wpt_1_y) ** 2)
            c = sqrt((wpt_2_x - wpt_1_x) ** 2 + (wpt_2_y - wpt_1_y) ** 2)

            # check if position is on the same line with wpt_1 and wpt_2
            # within 2 % margin
            margin = 2
            if c - ((c * margin) / 100) <= b + a <= c + ((c * margin) / 100):
                self.cur_wpt = i
                self.next_wpt = self.cur_wpt + 1
                self.init_success = True
                break

            # compute A and B angles of a triangle
            # which is another way of deciding
            # where aircraft is located with respect to the
            # waypoints of the route
            # the decision is made if aircraft's position is located
            # at certain A and B angles of a triangle
            val_for_A_acos = (b ** 2 + c ** 2 - a ** 2) / (2 * b * c)
            val_for_B_acos = (c ** 2 + a ** 2 - b ** 2) / (2 * c * a)

            # sometimes due to floating number precision
            # the values above can result in a
            # value larger than 1.0 which throws
            # a ValueError exception --> math domain error
            # since acos (below in A and B)
            # takes values in range [-1, 1]
            # to avoid this exception check if the values above are larger
            # than 1.0 and force them to be 1.0 and
            # if value is smaller than -1.0 force it to be -1.0
            if val_for_A_acos > 1.0:
                val_for_A_acos = 1.0
            elif val_for_A_acos < -1.0:
                val_for_A_acos = -1.0

            if val_for_B_acos > 1.0:
                val_for_B_acos = 1.0
            elif val_for_B_acos < -1.0:
                val_for_B_acos = -1.0

            A = degrees(acos(val_for_A_acos))
            B = degrees(acos(val_for_B_acos))
            # A = degrees(acos((b ** 2 + c ** 2 - a ** 2) / (2 * b * c)))
            # B = degrees(acos((c ** 2 + a ** 2 - b ** 2) / (2 * c * a)))

            if A <= 135 and B <= 46:
                self.cur_wpt = i
                self.next_wpt = self.cur_wpt + 1
                self.init_success = True
                break

        if self.init_success:
            self.pygmap_center = self.route[0]

            self.find_distances()
            self.find_bearings()

            self.track = self.course_array[self.cur_wpt]  # current course

            self.hdg = self.track

            self.course_array_degrees = [degrees(each) for each in self.course_array]

            if self.show_map:
                self.init_pygmap()

        else:
            print('prediciton could not be initialized')

    def create_txt_file(self):
        """writing results into txt file is much faster than writing to excel file"""
        self.results_file = open(self.path_to_save_xls+self.name+".txt", 'w')
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

    def get_distance_to_deceler_to_descent_mach(self):
        """ if mach number in descent is different from the mach number in cruise
        aircraft needs to decelerate to descent mach number
        descent must be performed in cruise

        to know when to start deceleration in cruise prior to top of descent
        precalculate the distance aircraft will travel in deceleration

        when this distance is reached before top of descent - start deceleration"""
        if self.mach_cruise == self.mach_descent:
            return 0

        else:
            tas_cr = mach_to_tas(self.mach_cruise, self.target_cr_alt, self.delta_T, self.lapse_rate)
            tas_des = mach_to_tas(self.mach_descent, self.target_cr_alt, self.delta_T, self.lapse_rate)

            time_ = abs(tas_des - tas_cr) / self.apm.acceler_long_max
            dist_ = (tas_cr * time_) + (self.apm.acceler_long_max * (time_**2))/2.  # m

            return dist_

    def convert_positions_to_cartesian(self, lat, lon):
        """ after all coordinates are calculated convert them to Cartesian coordinate system """

        lat_ref_rad = self.lat_ref_rad
        lon_ref_rad = self.lon_ref_rad

        lat_ref_deg = degrees(lat_ref_rad)
        lon_ref_deg = degrees(lon_ref_rad)

        lat_deg = degrees(lat)
        lon_deg = degrees(lon)

        # important not to mix y and x

        # for [nm] all coordinates must be in degrees and cos of lat in radians
        x = 60. * (lon_deg - lon_ref_deg) * cos(lat_ref_rad)  # [nm]
        y = 60. * (lat_deg - lat_ref_deg)  # [nm]

        # for [m]all coordinates must be in radians and cos of lat in radians
        # x = R_e * (lon_rad - lon_ref_rad) * cos(lat_ref_rad)  # [m]
        # y = R_e * (lat_rad - lat_ref_rad)  # [m]
        return x, y

    def update_position(self):

        self.d_s = self.GS * self.delta_time
        d_lat = self.lat_displacement()
        d_lon = self.long_displacement()

        self.cur_lat += d_lat
        self.cur_lon += d_lon

        self.flight_time += self.delta_time
        self.flight_time_counter += 1

        self.deg_lat_cur = degrees(self.cur_lat)
        self.deg_lon_cur = degrees(self.cur_lon)

    def main(self):

        while not self.LANDED:

            # ----------- CHECK IF EVENT OCCURRED ----------------------------------------------------------------------
            if self.tod_counter or self.toc_counter:  # check if this event has already happened
                self.event = 0
            # ----------------------------------------------------------------------------------------------------------

            # ------------ UPDATE ATMOSPHERIC CONDITIONS ---------------------------------------------------------------
            self.update_atmospheric_conditions()
            # ----------------------------------------------------------------------------------------------------------

            # ------------ DETERMINE PHASE -----------------------------------------------------------------------------
            if self.tod:  # if TOD is reached - its time to start descent
                self.phase = "DESCENT"

                if not self.tod_counter:
                    self.event = "TOD"
                    self.tod_counter += 1

            # elif self.alt_is_const:  # if level-off
            #     self.phase = "CRUISE"

            elif self.alt < self.target_cr_alt and self.phase != "DESCENT":
                self.phase = "CLIMB"

            elif self.alt > self.target_cr_alt and self.phase == "CLIMB":  # if target FL is reached - stay in cruise
                self.phase = "CRUISE"

                if not self.toc_counter:
                    self.event = "TOC"
                    self.toc_counter += 1

            elif self.alt < self.destination_alt and self.phase == "DESCENT":  # landed
                break
            # ----------------------------------------------------------------------------------------------------------

            # ---------- SPEED SETTINGS ACCORDING TO A FLIGHT PHASE ----------------------------------------------

            if self.phase == "CLIMB":
                self.climb_phase()

            elif self.phase == "CRUISE":
                self.cruise_phase()

            elif self.phase == "DESCENT":
                self.descent_phase()

            # -------------- UPDATE ALTITUDE AND FL --------------------------------------------------------------------
            self.alt += self.ROCD * self.delta_time

            if self.alt < 0.0:
                self.LANDED = True

            self.FL = int(round(self.alt * meter)) / 100  # convert altitude to FL

            # --------------- UPDATE FLIGHT PHASE FOR GOOGLE MAP --------------------------------
            if self.show_map:
                self.update_phase()
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------- UPDATE BEARING -------------------------------------------------------------------
            if self.TURN:
                """ integrate track """
                self.track += self.turn_rate * self.delta_time
                # self.time_in_turn += self.delta_time

            else:
                """ before calculating next lat/lon get update of a current track
                update is done based on current coordinate and target coordinate"""
                self.track = self.bearing(self.cur_lat, self.cur_lon,
                                          self.route[self.next_wpt][0],
                                          self.route[self.next_wpt][1])

            # ---------- UPDATE WIND, HEADING AND GROUND SPEED --------------------------------------------------------
            self.get_gamma()
            self.HAS = self.TAS * cos(self.gamma)  # same as self.HAS = sqrt(self.TAS**2 - self.ROCD**2)
            self.wind_components()
            self.u_array.append(self.u)
            self.v_array.append(self.v)
            self.wind_vector_from_components()
            self.get_wind_correction_angle()
            self.update_heading()
            self.update_ground_speed()
            # ----------------------------------------------------------------------------------------------------------

            # ------------- UPDATE DISTANCE TO DESTINATION AND TOP OF DESCENT ------------------------------------------
            if not self.TURN:
                self.dist_to_destination = (self.distance(self.cur_lat, self.cur_lon,
                                                          self.route[self.next_wpt][0],
                                                          self.route[self.next_wpt][1]) +
                                            sum(self.dist_array[self.next_wpt:]))

            # check distance to destination
            if (self.phase == "CRUISE" and
                    self.dist_to_destination - self.GS < self.descent_distance < self.dist_to_destination + self.GS or
                    self.dist_to_destination < self.descent_distance):
                self.tod = True

            if self.phase == "CRUISE" and self.dist_to_destination - self.GS < 10 + self.descent_distance + self.dist_to_reach_desc_mach:
                """ decelerate in cruise before descent"""
                self.mach_target = self.mach_descent
                self.reach_mach_for_descent = True
            # ----------------------------------------------------------------------------------------------------------

            # ------------- UPDATE POSITION --------------------------------------------------------------------------
            self.update_position()
            self.distance_traveled += self.d_s
            self.x, self.y = self.convert_flight_to_cartesian()

            # extract the necessary values for the prediction
            # comparison with the true flight
            if self.flight_time_counter == 5 * 60 * self.multp_factor:
                self.positions_at_lat[5] = [self.x, self.y, self.alt, self.FL, degrees(self.track)%360, self.phase,
                                            np.average(self.u_array), np.average(self.v_array)]
            elif self.flight_time_counter == 10 * 60 * self.multp_factor:
                self.positions_at_lat[10] = [self.x, self.y, self.alt, self.FL, degrees(self.track)%360, self.phase,
                                            np.average(self.u_array), np.average(self.v_array)]
            elif self.flight_time_counter == 15 * 60 * self.multp_factor:
                self.positions_at_lat[15] = [self.x, self.y, self.alt, self.FL, degrees(self.track)%360, self.phase,
                                            np.average(self.u_array), np.average(self.v_array)]
            elif self.flight_time_counter == 20 * 60 * self.multp_factor:
                self.positions_at_lat[20] = [self.x, self.y, self.alt, self.FL, degrees(self.track)%360, self.phase,
                                            np.average(self.u_array), np.average(self.v_array)]
                break

            if self.show_map:
                self.path_points.append((self.deg_lat_cur, self.deg_lon_cur))
            # ----------------------------------------------------------------------------------------------------------

            # ---------------- PERFORM CHECKS --------------------------------------------------------------------------
            if self.FLIGHT_END:
                # if we are approaching the last wp in our route -
                # check only for reaching required coordinates
                self.check_coordinates()
                self.check_altitude()
            elif not self.TURN:
                # if we are not performing the turn -
                # check distance to turn start
                # on every step check if we approached
                # the turning distance and need to start turning
                self.compute_turn_distance_and_rate()
                self.get_distance_to_wpt()
                self.check_turn_init_distance()
            elif self.TURN:
                # while we are performing turn - check bearing
                # to decide if we have already reached
                # required new heading
                # or we still need to keep turning
                self.check_track()
            # ----------------------------------------------------------------------------------------------------------

        return self.positions_at_lat

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
