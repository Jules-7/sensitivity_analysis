import os
from general import pygmaps
from math import atan2, cos, sin, degrees, radians, asin, sqrt, tan, pi

from general.constants import R_e, knot, g, ft, meter
from total_energy_model.conversion import cas_to_tas, tas_to_mach, crossover_altitude, mps_to_kts, mps_to_fpm, \
    tas_to_cas, knot_to_mps, stall_velocity, fpm_to_mps, mach_to_tas, mach_to_cas
from general.main_functions import Flight


class TrueFlightSimulation(Flight):

    def __init__(self, name, FPL, path_to_save_xls,
                 descent_distance, alt_start_decel):

        super(TrueFlightSimulation, self).__init__(name, FPL,
                                                   path_to_save_xls=path_to_save_xls,
                                                   descent_distance=descent_distance,
                                                   alt_start_decel=alt_start_decel)

        self.path_to_save_xls = path_to_save_xls
        self.descent_distance = descent_distance

        self.instruction = self.fpl['instruction']
        self.instr_type = self.instruction['type']
        self.instr_val = self.instruction['value']
        self.instr_dur = self.instruction['dur']

        # forgot to add it in inputs creation
        if self.instr_type == 'vs':
            self.instr_dur = 60.
        # -----------------------------------------------------------

        # ---------- GLOBAL VARIABLES -------------------------------
        self.phase = 'CLIMB'

        self.time_in_turn = 0
        self.distance_in_turn = 0
        self.radial_dist = 0
        self.leveloff_time = 0

        self.leveloff_climb_fl = self.leveloffs['CLIMB'][0]
        self.leveloff_climb_dur = self.leveloffs['CLIMB'][1]
        self.leveloff_descent_fl = self.leveloffs['DESCENT'][0]
        self.leveloff_descent_dur = self.leveloffs['DESCENT'][1]

        if alt_start_decel:
            self.alt_start_decel = alt_start_decel[0]
        else:
            self.alt_start_decel = 0

        # time counter of instruction duration
        self.instr_time = 0
        # -----------------------------------------------------------

        # ------------------- FLAGS ---------------------------------
        self.LEVELOFF = False
        self.tod = False  # top of descent

        # flag to indicate that aircraft started to perform an instruction
        self.INSTR_START_VECTOR, self.INSTR_START_SPEED, self.INSTR_START_CLIMB = False, False, False

        # flag to indicate that aircraft has completed the instruction
        self.INSTR_END_VECTOR, self.INSTR_END_SPEED, self.INSTR_END_CLIMB = False, False, False

        # flag to allow or cancel ATC climb instruction
        self.CANCEL_INSTR = False

        self.instr_dists_to_wpt = []
        self.instr_climb_dist = 0

        # ----- WARNINGS --------------------------------------------
        self.warn_about_distance = False
        # -----------------------------------------------------------
        if not self.test:
            # self.create_xls_file()
            self.create_txt_file()
        # -----------------------------------------------------------

        self.dist_to_reach_desc_mach = self.get_distance_to_reach_desc_mach()
        # print("distance to reach descent mach before TOD", self.dist_to_reach_desc_mach/1852)

        self.ROCD = self.VS_climb  # self.ROCD = self.VS_climb[0]

        # ------------ crossover altitudes --------------------------
        self.climb_crossover_altitude = crossover_altitude(self.cas_climb, self.mach_climb)
        self.descent_crossover_altitude = crossover_altitude(self.cas_descent, self.mach_descent)

        self.set_speeds(speed=self.cas_climb)

        self.init_route_position_track()

    def check_speeds(self):

        # check CAS
        if self.cas_climb < self.cas_min:
            print("CAS climb is less than min speed. Replacing with min speed\n")
            self.cas_climb = knot_to_mps(260)

        if self.cas_descent < self.cas_min:
            print("CAS descent is less than min speed. Replacing with min speed\n")
            self.cas_descent = knot_to_mps(260)

        if self.cas_climb > self.apm.V_mo:
            print("CAS climb is larger than max speed. Replacing with max speed\n")
            self.cas_climb = self.apm.V_mo

        if self.cas_descent > self.apm.V_mo:
            print("CAS descent is larger than max speed. Replacing with max speed\n")
            self.cas_descent = self.apm.V_mo

        # check Mach

        if self.mach_climb > self.apm.M_mo:
            print("Mach climb is larger than max %s. Replacing with max %s\n"%(self.wtc, self.wtc))
            self.mach_climb = self.apm.M_mo

        if self.mach_cruise > self.apm.M_mo:
            print("Mach cruise is larger than max M. Replacing with max M\n")
            self.mach_cruise = self.apm.M_mo

        if self.mach_descent > self.apm.M_mo:
            print("Mach descent is larger than max M. Replacing with max M\n")
            self.mach_descent = self.apm.M_mo

        # check ATC instruction if it is a speed instruction
        if self.instr_type == "speed":
            if self.instr_val > self.apm.M_mo:
                print("ATC instruction Mach is larger than max M. Replacing with max M\n")
                self.instr_val = self.apm.M_mo

        # ----- check if Mach is not smaller than minimum CAS -----------------

        cas_from_mach_climb = mach_to_cas(self.mach_climb, self.target_cr_alt,
                                              self.delta_T, self.lapse_rate)

        cas_from_mach_cruise = mach_to_cas(self.mach_cruise, self.target_cr_alt,
                                              self.delta_T, self.lapse_rate)

        cas_from_mach_descent = mach_to_cas(self.mach_descent, self.target_cr_alt,
                                              self.delta_T, self.lapse_rate)

        # check ATC instruction if it is a speed instruction
        if self.instr_type == "speed":
            cas_from_mach_atc = mach_to_cas(self.instr_val, self.target_cr_alt,
                                            self.delta_T, self.lapse_rate)

        if cas_from_mach_climb < self.cas_min:
            print("Mach climb is smaller than min speed. Replacing with min speed M\n")
            tas = cas_to_tas(knot_to_mps(260), self.target_cr_alt, self.delta_T,
                             self.lapse_rate)
            self.mach_climb = tas_to_mach(tas, self.target_cr_alt, self.delta_T,
                                          self.lapse_rate)

        if cas_from_mach_cruise < self.cas_min:
            print("Mach cruise is smaller than min speed. Replacing with min speed M\n")
            tas = cas_to_tas(knot_to_mps(260), self.target_cr_alt, self.delta_T,
                             self.lapse_rate)
            self.mach_cruise = tas_to_mach(tas, self.target_cr_alt, self.delta_T,
                                          self.lapse_rate)

        if cas_from_mach_descent < self.cas_min:
            print("Mach descent is smaller than min speed. Replacing with min speed M\n")
            tas = cas_to_tas(knot_to_mps(260), self.target_cr_alt, self.delta_T,
                             self.lapse_rate)
            self.mach_descent = tas_to_mach(tas, self.target_cr_alt, self.delta_T,
                                            self.lapse_rate)

        # check ATC instruction if it is a speed instruction
        if self.instr_type == "speed":
            if cas_from_mach_atc < self.cas_min:
                print("ATC instruction Mach is smaller than min speed. Replacing with min speed M\n")
                tas = cas_to_tas(knot_to_mps(260), self.target_cr_alt, self.delta_T,
                                 self.lapse_rate)
                self.instr_val = tas_to_mach(tas, self.target_cr_alt, self.delta_T,
                                             self.lapse_rate)

    def get_distance_to_reach_desc_mach(self):
        """ If Mach number in descent is different from the Mach number in cruise
        aircraft needs to decelerate or accelerate to descent mach number,
        acceleration/deceleration must be performed in cruise.

        To know when to start acceleration/deceleration in cruise prior to top of descent,
        pre-calculate the distance aircraft will travel while accelerating/decelerating.

        When this distance is reached before top of descent - start acceleration/deceleration"""

        if self.mach_cruise == self.mach_descent:
            return 0

        else:
            # check if this works when there is wind (GS != TAS)
            tas_cr = mach_to_tas(self.mach_cruise, self.target_cr_alt, self.delta_T, self.lapse_rate)
            tas_des = mach_to_tas(self.mach_descent, self.target_cr_alt, self.delta_T, self.lapse_rate)

            time_ = abs(tas_des - tas_cr) / self.apm.acceler_long_max
            dist_ = (tas_cr * time_) + (self.apm.acceler_long_max * (time_**2))/2.  # m

            return dist_

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

        # x, y = self.convert_flight_to_cartesian()

    def write_states_to_txt(self):

        data_to_record = [self.flight_time, self.FL, self.alt, mps_to_kts(self.CAS),
                         self.Mach, mps_to_kts(self.TAS),
                         mps_to_kts(self.GS), mps_to_kts(self.wind_speed),
                         mps_to_fpm(self.ROCD),
                         self.deg_lat_cur, self.deg_lon_cur,
                         self.x, self.y,
                         degrees(self.track)%360,
                         degrees(self.hdg)%360,
                         self.wind_direction_deg,
                         mps_to_kts(self.u), mps_to_kts(self.v),
                         self.phase,
                         self.TURN, degrees(self.bank_angle),
                         self.distance_traveled / 1852., self.event,
                         self.cur_wpt, self.next_wpt,
                         self.dist_to_destination / 1852.,
                         # mps_to_kts(self.u), mps_to_kts(self.v),
                         self.instr_type]

        data = ['%.3f'%each if not isinstance(each, str) else each for each in data_to_record]

        data.append('\n')

        data_row = ",".join(data)

        self.results_file.write(data_row)

    def cruise_phase(self):

        # ------------ CHECK SPEED SETTINGS -------------------------
        # above crossover altitude fly at constant Mach
        # since mach number for cruise may be different from the climb mach number
        # it means that when aircraft is levelled it needs
        # to accelerate/decelerate to reach cruise mach number
        # therefore we need to check what is our current mach number
        # and if we need - initiate acceleration/deceleration

        # this check ensures that correct mach number value is used
        # as a target mach value in acceleration/deceleration phase
        # in cruise prior to descent

        if self.INSTR_START_SPEED:
            self.mach_target = self.instr_val

            if self.instr_time >= self.instr_dur:
                self.INSTR_START_SPEED = False
                self.event = 'speed end'
                self.INSTR_END_SPEED = True

        elif not self.reach_mach_for_descent:
            # self.mach_target = get_mach(self.phase, self.apm, self.speed_profile)
            self.mach_target = self.mach_cruise

        self.target_spd = self.mach_target

        self.check_target_Mach()

        self.check_TAS_for_acc_dec(below=False)

        if self.ACCELERATE_TAS or self.DECELERATE_TAS:
            # if aircraft needs to accelerate - accelerate TAS
            # and check if target CAS has been reached
            # update acceleration/deceleration
            self.update_long_acc_dec()
            self.Mach = tas_to_mach(self.TAS, self.alt, self.delta_T, self.lapse_rate)
            self.CAS = tas_to_cas(self.TAS, self.alt, self.delta_T, self.lapse_rate)

        else:
            # if no acceleration is required - maintain constant Mach
            # convert value of constant Mach to TAS and CAS
            self.TAS = mach_to_tas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.CAS = mach_to_cas(self.Mach, self.alt, self.delta_T, self.lapse_rate)

        self.ROCD = 0.

    def check_leveloff_constraints(self):

        if self.phase == "CLIMB":
            if self.leveloff_climb_fl:
                leveloff_fl = self.leveloff_climb_fl[0]
                if leveloff_fl - 0.1 <= self.FL <= leveloff_fl + 0.1:
                    self.LEVELOFF = True

        elif self.phase == "DESCENT":
            if self.leveloff_descent_fl:
                leveloff_fl = self.leveloff_descent_fl[0]
                if leveloff_fl - 0.1 <= self.FL <= leveloff_fl + 0.1:
                    self.LEVELOFF = True

    def update_leveloff_time(self):
        self.leveloff_time += self.delta_time

        if self.CLIMB and self.leveloff_climb_fl:
            if self.leveloff_time >= self.leveloff_climb_dur[0]:
                # self.temp_phase = ""
                self.LEVELOFF = False
                del self.leveloff_climb_fl[0]
                del self.leveloff_climb_dur[0]
                self.leveloff_time = 0

        if self.DESCENT and self.leveloff_descent_fl:
            if self.leveloff_time >= self.leveloff_descent_dur[0]:
                # self.temp_phase = ""
                self.LEVELOFF = False
                del self.leveloff_descent_fl[0]
                del self.leveloff_descent_dur[0]
                self.leveloff_time = 0

    def update_instruction_time(self):
        self.instr_time += self.delta_time

    def main(self):

        while not self.LANDED:

            # ----------- CHECK IF EVENT OCCURRED -------------------
            # check if this event has already happened
            if self.tod_counter or self.toc_counter or self.instr_counter or self.decel_counter:
                self.event = ""
            # -------------------------------------------------------

            # ------------ UPDATE ATMOSPHERIC CONDITIONS ------------
            self.update_atmospheric_conditions()
            # -------------------------------------------------------

            # ------------ DETERMINE PHASE --------------------------
            if self.tod:  # if TOD is reached - its time to start descent
                self.phase = "DESCENT"
                self.CLIMB = False
                self.CRUISE = False
                self.DESCENT = True

                if not self.tod_counter:
                    self.event = "TOD"
                    self.tod_counter += 1

            elif self.alt < self.target_cr_alt and self.phase != "DESCENT":
                self.phase = "CLIMB"
                self.CLIMB = True
                self.CRUISE = False
                self.DESCENT = False

            elif self.alt > self.target_cr_alt and self.phase == "CLIMB":  # if target FL is reached - stay in cruise
                self.phase = "CRUISE"
                self.CLIMB = False
                self.CRUISE = True
                self.DESCENT = False

                if not self.toc_counter:
                    self.event = "TOC"
                    self.toc_counter += 1

            elif self.alt < self.destination_alt and self.phase == "DESCENT":  # landed
                break

            if self.flight_time > 65000.:  # if there is somewhere an infinite loop
                break
            # -------------------------------------------------------

            # ----- SPEED SETTINGS ACCORDING A TO FLIGHT PHASE ------

            if self.phase == "CLIMB":
                self.climb_phase()

            elif self.phase == "CRUISE":
                self.cruise_phase()

            elif self.phase == "DESCENT":
                self.descent_phase()

            # set speed according to a flight phase
            # and set ROCD to 0 if performing
            # a temporary level-off
            if self.LEVELOFF:
                self.ROCD = 0.

            # check for ATC speed and climb instruction in cruise phase
            if self.toc_counter and self.instr_type == "speed" and not self.INSTR_START_SPEED and not self.INSTR_END_SPEED:
                """ speed instruction """
                self.INSTR_START_SPEED = True
                self.instr_counter += 1
                self.event = "speed start"

            if (self.toc_counter and self.instr_type == "vs" and not self.instr_climb_dist and
                    not self.CANCEL_INSTR and not self.INSTR_END_CLIMB):
                """ determine distance and time we have in cruise"""
                cruise_dist = self.dist_to_destination - self.descent_distance
                time_in_cruise = cruise_dist / self.GS
                # if we have more than 2.5 minutes in cruise perform ATC climb instruction
                if time_in_cruise > 2.5*60.:
                    self.instr_climb_dist = cruise_dist / 3
                else:
                    # if there is not enough time for instruction, cancel it and dont check for it anymore
                    self.CANCEL_INSTR = True

            if self.INSTR_START_CLIMB:
                self.ROCD = fpm_to_mps(1000)

                if self.instr_time >= self.instr_dur:
                    self.INSTR_START_CLIMB = False
                    self.instr_counter += 1
                    self.event = "climb end"
                    self.INSTR_END_CLIMB = True

            # ----- UPDATE ALTITUDE, FL AND CHECK FOR TEMPORARY LEVEL-OFFS ----
            self.alt += self.ROCD * self.delta_time

            if self.alt < 0.0:
                self.LANDED = True

            # convert altitude to FL
            self.FL = int(round(self.alt * meter) / 100)

            if not self.LEVELOFF:
                self.check_leveloff_constraints()
            # -------------------------------------------------------

            # ----- UPDATE FLIGHT PHASE FOR GOOGLE MAP --------------
            if self.show_map:
                self.update_phase()

            if self.INSTR_START_VECTOR or self.INSTR_START_SPEED or self.INSTR_START_CLIMB:
                self.update_instruction_time()

            # ----------------------- UPDATE BEARING ----------------
            if self.TURN:
                """ integrate track """
                self.track += self.turn_rate * self.delta_time
                self.time_in_turn += self.delta_time

            elif (self.toc_counter and self.instr_type == "heading"
                  and not self.INSTR_START_VECTOR and not self.INSTR_END_VECTOR):
                """ if aircraft is turning at TOC do not start vector.
                Let aircraft finish the turn and then start vectoring instruction"""
                self.track += radians(self.instr_val)
                self.INSTR_START_VECTOR = True
                self.instr_counter += 1
                self.event = "vector start"

            elif self.INSTR_START_VECTOR and self.instr_type == "heading":

                if self.instr_time >= self.instr_dur:
                    self.INSTR_START_VECTOR = False
                    self.event = 'vector end'
                    self.INSTR_END_VECTOR = True

                    if self.show_map:
                        self.flight_phases.append((self.phase, self.cur_lat, self.cur_lon))

                else:
                    # check and update next waypoint if necessary
                    # if the waypoint aircraft was flying towards is behind, there is no point
                    # in flying back to it. Instead, switch to the next waypoint
                    self.dist_to_wp = self.distance(self.cur_lat, self.cur_lon,
                                                    self.route[self.next_wpt][0], self.route[self.next_wpt][1])

                    self.instr_dists_to_wpt.append([self.dist_to_wp, self.next_wpt])
                    if len(self.instr_dists_to_wpt) >= 2:
                        if (self.instr_dists_to_wpt[-1][1] == self.instr_dists_to_wpt[-2][1] and
                                self.instr_dists_to_wpt[-1][0] > self.instr_dists_to_wpt[-2][0]):
                            self.cur_wpt = self.next_wpt
                            self.next_wpt += 1

                """ It is better not to return aircraft back after vectoring is completed.
                The best option is to let it fly towards the next waypoint downstream.
                This is the most flexible way. When returning aircraft back (e.g. X degrees)
                sometimes this results in no flight completion. Keep it as is. """

            # elif self.INSTR_END and not self.BACK_ON_ROUTE:
            #
            #     brng_to_next_wpt = degrees(self.bearing(self.cur_lat, self.cur_lon,
            #                                             self.route[self.next_wpt][0], self.route[self.next_wpt][1]))
            #     if self.counter == 0:
            #         track_deg = degrees(self.track)
            #         delta_angle = self.instr_val * 2
            #         if self.instr_val < 0:
            #             self.track += radians(abs(delta_angle))
            #         else:
            #             self.track -= radians(delta_angle)
            #         self.counter += 1
            #
            #     track_deg = degrees(self.track)
            #
            #     next_wp_course = degrees(self.course_array[self.cur_wpt])
            #
            #     brng_to_next_wpt_ = degrees(self.bearing(self.cur_lat, self.cur_lon,
            #                                             self.route[self.next_wpt][0], self.route[self.next_wpt][1]))
            #
            #     low_lim = next_wp_course - 0.5
            #     up_lim = next_wp_course + 0.5
            #
            #     if low_lim < brng_to_next_wpt < up_lim:
            #         self.BACK_ON_ROUTE = True

            else:
                """ before calculating next lat/lon get update of a current track
                update is done based on current coordinate and target coordinate"""
                self.track = self.bearing(self.cur_lat, self.cur_lon,
                                          self.route[self.next_wpt][0],
                                          self.route[self.next_wpt][1])

            # ----- UPDATE WIND, HEADING AND GROUND SPEED -----------------------
            self.get_gamma()
            self.HAS = self.TAS * cos(self.gamma)  # same as self.HAS = sqrt(self.TAS**2 - self.ROCD**2)
            self.wind_components()
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

            """ TOD location. 
            When aircraft reaches the distance to destination it takes to descent (descent distance),
            initiate descent"""

            # this is a warning that the present distance to destination
            # is smaller then the distance necessary for descent
            if not self.warn_about_distance:
                if self.dist_to_destination < self.descent_distance:
                    # print("distance to destination is smaller than descent distance")
                    self.warn_about_distance = True

            if (self.phase == "CRUISE" and
                    self.dist_to_destination - self.GS < self.descent_distance < self.dist_to_destination + self.GS or
                    self.dist_to_destination < self.descent_distance or
                    self.dist_to_destination - self.distance_in_turn < self.descent_distance):
                # self.TOD = self.top_of_descent()
                self.tod = True

            """ Deceleration segment
            If aircraft needs to decelerate/accelerate to descent mach, it is performed while in cruise.
            Start deceleration/acceleration prior to TOD when the deceleration/acceleration distance
            is reached (prior to TOD location)"""

            if (self.phase == "CRUISE" and
                    self.dist_to_destination - self.GS < 10 + self.descent_distance + self.dist_to_reach_desc_mach):
                self.mach_target = self.mach_descent
                self.reach_mach_for_descent = True
                self.decel_counter += 1
                self.event = "decel for descent"

            if (self.phase == "CRUISE" and self.instr_type == "vs" and
                    self.instr_climb_dist and not self.INSTR_START_CLIMB and not self.INSTR_END_CLIMB and
                    not self.CANCEL_INSTR):
                if self.dist_to_destination <= self.descent_distance + self.instr_climb_dist * 2:
                    self.INSTR_START_CLIMB = True
                    self.instr_counter += 1
                    self.event = "climb start"

            # ------------- UPDATE POSITION --------------------------------------------------------
            self.update_position()
            self.distance_traveled += self.d_s

            self.x, self.y = self.convert_flight_to_cartesian()
            # ----------------------------------------------------------------------------------------------------------

            if self.show_map:
                self.path_points.append((self.deg_lat_cur, self.deg_lon_cur))

            # write states of the flight if the time step is reached
            if not self.test:
                # writing to txt file is much faster than writing to xls file
                # self.write_states_to_xls()
                if self.flight_time_counter == 1 or self.flight_time_counter % 600 == 0:
                    self.write_states_to_txt()

            if self.LEVELOFF:
                self.update_leveloff_time()

            # ---------------- PERFORM CHECKS --------------------------------------------------------------------------
            if self.FLIGHT_END:
                # if we are approaching the last wp in our route -
                # check only for reaching required coordinates
                self.check_coordinates()
                self.check_altitude()
            # check if the second 'not' works
            elif not self.TURN and not self.INSTR_START_VECTOR:
                # if we are not performing the turn -
                # check distance to turn start
                # on every step check if we approached
                # the turning distance and need to start turning
                self.time_in_turn = 0
                self.distance_in_turn = 0
                self.radial_dist = 0
                self.compute_turn_distance_and_rate()
                self.get_distance_to_wpt()
                self.check_turn_init_distance()
            elif self.TURN:
                # ToDo make working code on computing distance to destination while turning
                # otherwise such events as TOD or deceleration for descent do not
                # happen on time when aircraft is turning
                # since distance computation is not performed in turn
                # and the last computed distance is used

                # while we are performing turn - check bearing
                # to decide if we have already reached
                # required new heading
                # or we still need to keep turning
                self.check_track()

                self.distance_in_turn += self.d_s

                # self.radial_dist = abs(self.turn_rate)
                # self.distance_in_turn = (degrees(self.radial_dist)/360) * (2 * pi * self.r)

                # self.get_distance_to_wpt()
                # self.check_waypoint_switch()
            # ----------------------------------------------------------------------------------------------------------

        # print altitude to see if aircraft actually landed
        # print("flight ends at altitude of %s m" % self.alt)
        # print("travelled distance %s nm" % (self.distance_traveled / 1852.))

        # final check of altitude, position and distance travelled to mark
        # flight as OK or NOT OK

        cur_lat = round(degrees(self.cur_lat), 4)
        next_lat = round(degrees(self.route[self.next_wpt][0]), 4)
        cur_lon = round(degrees(self.cur_lon), 4)
        next_lon = round(degrees(self.route[self.next_wpt][1]), 4)

        close_enough_lat = next_lat - 0.25 < cur_lat <= next_lat + 0.25
        close_enough_lon = next_lon - 0.25 < cur_lon <= next_lon + 0.25

        close_enough_alt = -10 <= self.alt <= 500

        # if not self.instr_type == 'heading':
        #     close_enough_dist = self.route_dist - 3 <= (self.distance_traveled / 1852.) <= self.route_dist + 3
        # elif self.instr_type == 'heading':
        #     close_enough_dist = self.route_dist - 30 <= (self.distance_traveled / 1852.) <= self.route_dist + 30

        dist_to_dest_wpt = (self.distance(self.cur_lat, self.cur_lon,
                                         radians(next_lat), radians(next_lon)))/ 1852.
        if not self.instr_type == "heading":
            close_enough_dist = dist_to_dest_wpt <= 3
        elif self.instr_type == "heading":
            close_enough_dist = dist_to_dest_wpt <= 30

        time_is_not_too_large = self.flight_time < 64999.

        return close_enough_lat and close_enough_lon and close_enough_alt and close_enough_dist and time_is_not_too_large
