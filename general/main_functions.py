from general import pygmaps
from math import atan2, cos, sin, degrees, radians, asin, sqrt, tan, pi

from general.atmosphere import Atmosphere
from general.constants import R_e, knot, g, ft, meter
from total_energy_model.conversion import cas_to_tas, tas_to_mach, crossover_altitude, mps_to_kts, mps_to_fpm, \
    tas_to_cas, knot_to_mps, stall_velocity, fpm_to_mps, mach_to_tas, mach_to_cas


class Flight(object):

    def __init__(self, name, FPL, **kwargs):
        self.name = name
        self.fpl = FPL
        self.apm = self.fpl['apm']
        self.target_cr_FL = self.fpl['FL']
        self.target_cr_alt = self.fpl['FL'] * 100 * ft  # [m]
        self.route = self.fpl['route']
        self.lapse_rate = self.fpl['lapse_rate']
        self.wtc = self.fpl['wtc']

        # self.path_to_save_xls = kwargs['path_to_save_xls']
        # self.descent_distance = kwargs['descent_distance']

        self.route_dist = self.fpl['route_dist']
        self.route_id = self.fpl['route_id']

        self.VS = self.fpl['VS']
        self.VS_climb = fpm_to_mps(self.VS['Climb'])
        self.VS_descent = fpm_to_mps(self.VS['DescentPrediction'])

        self.bank_angles = self.fpl['bank']
        self.speeds = self.fpl['speeds']
        self.leveloffs = self.fpl['leveloff']

        self.delta_T = float(self.fpl['temp_diff'])  # temperature difference from ISA [K]

        self.show_map = self.fpl['show_map']
        self.test = self.fpl['test']
        self.wind_comp = self.fpl['wind_comp']
        self.a_u, self.b_u, self.c_u = self.wind_comp['u']
        self.a_v, self.b_v, self.c_v = self.wind_comp['v']
        self.u, self.v = 0, 0
        # ---------------------------------------------------------------------------------

        # ---------- GLOBAL VARIABLES -----------------------------------------------------
        self.phase = None
        # self.temp_phase = ''
        self.destination_alt = 0  # [m] altitude at landing for now assume 0 - no elevation at aerodrome
        self.event = ""
        self.toc_counter, self.tod_counter, self.turn_counter, self.instr_counter, self.decel_counter = 0, 0, 0, 0, 0

        self.delta_time = 0.1  # [sec]
        self.flight_time = 0  # [sec]
        self.flight_time_counter = 0
        self.cur_wpt = 0  # [index]
        self.next_wpt = 1  # [index]

        self.distance_traveled = 0  # [m]
        self.alt, self.target_alt, self.FL = 0, 0, 0

        self.x, self.y = None, None

        self.d_s = 0
        self.dist_to_destination = 0  # [m]

        self.cas_target, self.mach_target, self.rocd_target, self.v_stall, self.target_spd = 0, 0, 0, 0, 0
        self.alt_max = 0

        self.cas_climb = knot_to_mps(self.speeds['CLIMB'][0])
        self.mach_climb = self.speeds['CLIMB'][1]

        self.mach_cruise = self.speeds['CRUISE'][0]

        self.cas_descent = knot_to_mps(self.speeds['DESCENT'][0])
        self.mach_descent = self.speeds['DESCENT'][1]

        self.acc_long = self.apm.acceler_long_max
        self.acc_norm = self.apm.acceler_norm_max

        self.dist_to_wp = 0
        self.bank_angle = 0
        # -----------------------------------------------------------

        # ------------------- FLAGS ---------------------------------
        self.tod = False  # top of descent

        # while phase can be Climb, aircraft can be in level (cruising) in climb
        self.CLIMB, self.CRUISE, self.DESCENT, self.LANDED = True, False, False, False

        self.ACCELERATE_TAS = False
        self.DECELERATE_TAS = False
        self.TARGET_CAS_REACHED = False  # check while changing CAS if target CAS is reached
        self.TARGET_Mach_REACHED = False

        self.reach_mach_for_descent = False

        self.CAS_changed, self.Mach_changed = False, False

        self.TURN, self.TURN_LEFT, self.TURN_RIGHT, self.FLIGHT_END = False, False, False, False

        # ----- WARNINGS --------------------------------------------
        self.warn_about_distance = False
        # -----------------------------------------------------------

        # -----------------------------------------------------------
        local_atm = Atmosphere(self.alt, self.delta_T, self.lapse_rate)  # air temperature, pressure, density
        self.T, self.press, self.density = local_atm.temp, local_atm.pressure, local_atm.density
        # -----------------------------------------------------------

        self.wind_speed = 0
        self.wind_direction = 0
        self.wind_direction_deg = 0
        # -------------------------------------------------

        # ------ ANGLES -------------------------------
        self.gamma = 0  # climb/descent angle
        self.fpa = 0  # flight path angle
        self.WCA = 0  # wind correction angle
        self.hdg = 0  # heading is a direction in which aircraft nose is pointing - True Heading
        self.track = 0  # is a path aircraft is following over the ground - True Track
        # # track and course are the same
        self.trk_diff = 0
        # -----------------------------------------------

        # ------------- SPEEDS -----------------------------
        # speeds values are set in the child class
        self.CAS, self.TAS, self.Mach, self.GS, self.ROCD, self.HAS = None, None, None, None, None, None

        # ------------ crossover altitudes -----------------------------------------------------------------------------
        # crossover altitudes are computed in the child class
        self.climb_crossover_altitude, self.descent_crossover_altitude = None, None
        # --------------------------------------------------------------------------------------------------------------

        self.course_array = []  # initial courses/track at each waypoint
        self.dist_array = []  # distance between waypoints list
        self.target_CAS_array = [0]  # collect all target CAS to check if acceleration/deceleration is needed
        self.target_Mach_array = [0]
        self.target_ROCD_array = [0]  # collect all target ROCD to check if acceleration/deceleration is needed
        self.turn_array = []

        # self.current_heading = 0
        # self.next_heading = 0
        self.current_track, self.next_wp_course = 0, 0
        self.r = 0  # turn radius
        self.turn_rate, self.turn_initiation_distance = 0, 0
        self.course_change_angle = 0
        # self.heading_change_angle = 0
        self.cur_lat, self.cur_lon = None, None
        self.lat_ref_rad, self.lon_ref_rad = None, None
        self.deg_lat_cur, self.deg_lon_cur = None, None
        self.pygmap_center = None

    def set_speeds(self, speed):

        if speed > knot_to_mps(250):
            self.CAS = knot_to_mps(250)
        else:
            self.CAS = speed

        self.TAS = cas_to_tas(self.CAS, self.alt, self.delta_T, self.lapse_rate)
        self.Mach = tas_to_mach(self.TAS, self.alt, self.delta_T, self.lapse_rate)
        self.GS = self.TAS

    def init_route_position_track(self):

        self.get_bank_angle()  # [rad]

        # ------------------- INITIAL POSITION -----------------------------------------------------------------------
        self.cur_lat = self.route[self.cur_wpt][0]
        self.cur_lon = self.route[self.cur_wpt][1]

        self.lat_ref_rad, self.lon_ref_rad = self.cur_lat, self.cur_lon
        self.deg_lat_cur = degrees(self.cur_lat)
        self.deg_lon_cur = degrees(self.cur_lon)

        self.pygmap_center = self.route[0]

        self.find_distances()
        self.find_bearings()

        # -------------- READABLE VALUES ---------------------------
        self.course_array_degrees = [degrees(each) for each in self.course_array]
        # ------------------------------------------------------------------------

        self.track = self.course_array[self.cur_wpt]  # current course
        self.hdg = self.track

        if self.show_map:
            self.init_pygmap()

    def create_xls_file(self):

        return

    def create_txt_file(self):
        """ Writing results into txt file is much faster than writing to excel file"""
        self.results_file = open(self.path_to_save_xls + self.name + "_" + "rt_" + str(self.route_id) + ".txt", 'w')
        header = ",".join(["TIME_s", "FL", "ALT_m", "CAS_kt",
                           "M", "TAS_kt",
                           "GS_kt", "WS_kt", "ROCD_fpm",
                           "LAT_DD", "LONG_DD", "X_nm", "Y_nm",
                           "TRK_deg", "HDG_deg",
                           "WD_deg",
                           "U_kt", "V_kt", "PHASE",
                           "TURN", "BANK_deg",
                           "DIST_nm", "EVENT",
                           "CUR_wpt", "NEXT_wpt",
                           "DIST_to_dest", "ATC", "\n"])
        self.results_file.write(header)

        return

    def init_pygmap(self):
        # --------------- PYGMAPS -----------------
        self.gmap = pygmaps.pygmaps(degrees(self.pygmap_center[0]), degrees(self.pygmap_center[1]), 5)
        self.path_points = []
        self.path_waypoints = []
        self.flight_phases = [(self.phase, self.cur_lat, self.cur_lon)]
        for each in self.route:
            # show waypoints
            self.gmap.addpoint(degrees(each[0]), degrees(each[1]))
            # connect waypoints with lines
            self.path_waypoints.append((degrees(each[0]), degrees(each[1])))
        # display connection between waypoints
        self.gmap.addpath(self.path_waypoints, "00FF00")

    def wind_components(self):
        self.u = knot_to_mps(self.a_u * self.FL ** 2 + self.b_u * self.FL + self.c_u)
        self.v = knot_to_mps(self.a_v * self.FL ** 2 + self.b_v * self.FL + self.c_v)

    def speed_from_components(self):
        """ Compute wind speed from wind components
        equations are from Earth Observing Laboratory website
        https://www.eol.ucar.edu/content/wind-direction-quick-reference

        :param u: component represents wind blowing to the East (west-east direction)
        :param v: component represents wind blowing to the North (north-south direction)
        :return: wind speed in units of u and v
        """
        self.wind_speed = sqrt(self.u ** 2 + self.v ** 2)

    def direction_from_components(self):
        """ Compute wind direction from wind components
        equations are from Earth Observing Laboratory website
        https://www.eol.ucar.edu/content/wind-direction-quick-reference

        :param u: component represents wind blowing to the East (west-east direction)
        :param v: component represents wind blowing to the North (north-south direction)
        :return: wind direction in degrees
        """
        self.wind_direction = atan2(-self.u, -self.v)

        # to get direction in degrees from 0 to 359 use the next return
        self.wind_direction_deg = (degrees(self.wind_direction) + 360) % 360

    def wind_vector_from_components(self):
        self.speed_from_components()
        self.direction_from_components()

    def get_max_altitude(self):
        self.alt_max = min(self.apm.h_mo,
                           (self.apm.h_max + self.apm.G_t * (self.delta_T - self.apm.C_tc_4) +
                            self.apm.G_w * (self.apm.mass_max - self.mass)))

    def distance(self, lat_1, long_1, lat_2, long_2):
        """ Input arguments are already in radians,
        returned distance is in meters """
        # haversine
        delta_lat = float(lat_2) - float(lat_1)
        delta_lon = float(long_2) - float(long_1)

        root = (sin(delta_lat / 2.) * sin(delta_lat / 2.) + cos(lat_1) * cos(lat_2) *
                sin(delta_lon / 2.) * sin(delta_lon / 2.))

        sigma = 2. * atan2(sqrt(root), sqrt(1. - root))
        dist = R_e * sigma
        return dist

    def bearing(self, lat_1, lon_1, lat_2, lon_2):
        """ Input arguments are already in radians,
         after bearing calculation - it is checked on the + or - sign.
        If it is negative, i.e. -90 deg, it is converted into positive angle, i.e. 270 deg """

        delta_lon = float(lon_2) - float(lon_1)

        x = sin(delta_lon) * cos(lat_2)

        y = (cos(lat_1) * sin(lat_2) -
             sin(lat_1) * cos(lat_2) * cos(delta_lon))

        brng = atan2(x, y)

        if degrees(brng) < 0:
            brng = 360 + degrees(brng)
            brng = radians(brng)

        return brng

    def find_distances(self):
        """ Calculate distances between each pair of 2 waypoints"""
        i = 0
        while i < len(self.route) - 1:
            dist = self.distance(self.route[i][0],
                                 self.route[i][1],
                                 self.route[i + 1][0],
                                 self.route[i + 1][1])
            self.dist_array.append(dist)
            i += 1

    def find_bearings(self):
        """ Calculate bearings (courses) between each pair of 2 waypoints"""
        i = 0
        while i < len(self.route) - 1:
            brg = self.bearing(self.route[i][0],
                               self.route[i][1],
                               self.route[i + 1][0],
                               self.route[i + 1][1])
            self.course_array.append(brg)

            i += 1

    def update_long_acc_dec(self):

        if self.ACCELERATE_TAS:
            self.TAS += self.acc_long * self.delta_time

        if self.DECELERATE_TAS:
            self.TAS -= self.acc_long * self.delta_time

    def update_atmospheric_conditions(self):
        local_atm = Atmosphere(self.alt, self.delta_T, self.lapse_rate)  # air temperature, pressure, density
        self.T, self.press, self.density = local_atm.temp, local_atm.pressure, local_atm.density

    def check_TAS_for_acc_dec(self, below, descent_calc=False):
        """ Raise flag to accelerate/decelerate if necessary based on CAS """
        if below:
            if self.CAS > self.cas_target and self.ACCELERATE_TAS and not self.TARGET_CAS_REACHED:
                self.TARGET_CAS_REACHED = True
                self.ACCELERATE_TAS = False
                if descent_calc:
                    self.descent_deccel_alt.append(self.alt)

            elif self.CAS < self.cas_target and self.DECELERATE_TAS and not self.TARGET_CAS_REACHED:
                self.TARGET_CAS_REACHED = True
                self.DECELERATE_TAS = False
                # if descent_calc:
                #     self.descent_deccel_alt.append(self.alt)

            elif self.CAS < self.cas_target and not self.TARGET_CAS_REACHED:
                self.ACCELERATE_TAS = True

            elif self.CAS > self.cas_target and not self.TARGET_CAS_REACHED:
                self.DECELERATE_TAS = True

            else:
                self.ACCELERATE_TAS = False
                self.DECELERATE_TAS = False

        else:
            if self.Mach > self.mach_target and self.ACCELERATE_TAS and not self.TARGET_Mach_REACHED:
                self.TARGET_Mach_REACHED = True
                self.ACCELERATE_TAS = False

            elif self.Mach < self.mach_target and self.DECELERATE_TAS and not self.TARGET_Mach_REACHED:
                self.TARGET_Mach_REACHED = True
                self.DECELERATE_TAS = False

            elif self.Mach < self.mach_target and not self.TARGET_Mach_REACHED:
                self.ACCELERATE_TAS = True

            elif self.Mach > self.mach_target and not self.TARGET_Mach_REACHED:
                self.DECELERATE_TAS = True

            else:
                self.ACCELERATE_TAS = False
                self.DECELERATE_TAS = False

    def get_gamma(self):
        """ Calculate current flight path angle.

        Angle of climb/descent is the angle in the air frame of reference
        and the flight path angle is the effect of the wind on the aircraft
        and the actual angle it flies relative to the ground.

        Headwind increases the flight path angle, makes it steeper,
        whereas tailwind decreases the flight path angle, makes it shallower.

        see http://www.atpforum.eu/forum/technical-subjects/-032-performance-a/110992-angle-of-descent-versus-flight-path-angle-320065
        """
        try:
            self.fpa = asin(self.ROCD / self.GS)  # [rad]  ground frame of reference
            self.gamma = asin(self.ROCD / self.TAS)  # [rad]  climb/descent angle
        except ValueError:
            print("error in fpa/gamma computation at %s m at %.2f seconds" %(self.alt, self.flight_time))

    def vertical_displacement(self):
        d_h = self.ROCD * self.delta_time
        return d_h

    def horiz_lat_displacement(self):
        d_s_lat = self.d_s * cos(self.track)
        return d_s_lat

    def horiz_long_displacement(self):
        d_s_long = self.d_s * sin(self.track)
        return d_s_long

    def lat_displacement(self):
        d_s_lat = self.horiz_lat_displacement()
        d_lat = d_s_lat / R_e
        return d_lat

    def long_displacement(self):
        d_s_long = self.horiz_long_displacement()
        d_long = (d_s_long / R_e) * 1 / cos(self.cur_lat)
        return d_long

    def update_phase(self):
        if self.flight_phases[-1][0] != self.phase:
            self.flight_phases.append((self.phase, self.cur_lat, self.cur_lon))

    def convert_flight_to_cartesian(self, units='nm'):
        """ Convert lat/long to XY (Cartesian coordinate system) """

        lat_ref_rad = self.lat_ref_rad
        lon_ref_rad = self.lon_ref_rad

        lat_ref_deg = degrees(lat_ref_rad)
        lon_ref_deg = degrees(lon_ref_rad)

        lat_deg = self.deg_lat_cur
        lon_deg = self.deg_lon_cur

        # important not to mix y and x
        if units == 'nm':
            # for [nm] all coordinates must be in degrees and cos of lat in radians
            x = 60. * (lon_deg - lon_ref_deg) * cos(lat_ref_rad)  # [nm]
            y = 60. * (lat_deg - lat_ref_deg)  # [nm]

        elif units == 'm':
            lat_rad = self.cur_lat
            lon_rad = self.cur_lon

            # for [m]all coordinates must be in radians and cos of lat in radians
            x = R_e * (lon_rad - lon_ref_rad) * cos(lat_ref_rad)  # [m]
            y = R_e * (lat_rad - lat_ref_rad)  # [m]
        return x, y

    def check_target_CAS(self):
        """ If next (target) CAS is different from the current one append it in a list and raise flag that target CAS
           is not reached, which will result in acceleration/deceleration of TAS """
        if self.target_CAS_array[-1] != self.cas_target:
            if self.target_CAS_array[-1] - 0.1 <= self.cas_target <= self.target_CAS_array[-1] + 0.1:
                pass
            else:
                self.target_CAS_array.append(self.cas_target)
                self.TARGET_CAS_REACHED = False
                self.CAS_changed = True

    def check_target_Mach(self):
        """ If next (target) Mach is different from the current one append it in a list and raise flag that target Mach
           is not reached, which will result in acceleration/deceleration of TAS """
        if self.target_Mach_array[-1] != self.mach_target:
            self.target_Mach_array.append(self.mach_target)
            self.TARGET_Mach_REACHED = False
            self.Mach_changed = True

    def update_norm_acc_dec(self):

        if self.ACCELERATE_ROCD:
            self.ROCD += self.acc_norm * self.delta_time

        if self.DECELERATE_ROCD:
            self.ROCD -= self.acc_norm * self.delta_time

    def check_target_ROCD(self):
        if self.target_ROCD_array[-1] != self.rocd_target:
            self.target_ROCD_array.append(self.rocd_target)
            self.TARGET_ROCD_REACHED = False

    def update_position(self):

        self.d_s = self.GS * self.delta_time
        d_lat = self.lat_displacement()
        d_lon = self.long_displacement()

        self.cur_lat += d_lat
        self.cur_lon += d_lon

        self.flight_time += self.delta_time

        self.deg_lat_cur = degrees(self.cur_lat)
        self.deg_lon_cur = degrees(self.cur_lon)

        # x, y = self.convert_flight_to_cartesian()

    def write_states_to_xls(self):
        pass

    def get_distance_to_wpt(self):
        # get distance to wp
        self.dist_to_wp = self.distance(self.cur_lat, self.cur_lon,
                                        self.route[self.next_wpt][0],
                                        self.route[self.next_wpt][1])

    def check_turn_init_distance(self):
        """  Check our current distance to waypoint we are heading to """
        # check if its time to start the turn
        # sometimes, if different bank angles are used in different phases
        # and a turn calculation is made in one phase with one bank angle
        # while the turn will be performed in another phase with
        # a different bank angle, as soon as aircraft is in another phase
        # the turn initiation distance can be very different (smaller/larger)
        # due to radius due to bank angle
        # and it can happen that with the new bank angle the distance to WPT
        # left is smaller then the required turn initiation distance
        # in that case start turning right away
        if (self.turn_initiation_distance - self.GS * 0.2 <= self.dist_to_wp <=
                self.turn_initiation_distance + self.GS * 0.2 or
                self.dist_to_wp < self.turn_initiation_distance):
            self.TURN = True
            self.turn_array.append((self.phase, self.cur_lat, self.cur_lon))

    def check_track(self):
        """ While turning we need to make a decision
        if we have reached a new heading or we need to keep turning,
        the new heading is corrected for the wind if its known """

        # self.current_heading = (degrees(self.hdg))%360
        self.current_track = degrees(self.track) % 360

        self.next_wp_course = degrees(self.course_array[self.next_wpt])

        low_lim = self.next_wp_course - 1.5
        up_lim = self.next_wp_course + 1.5

        if low_lim < self.current_track < up_lim:

            if self.TURN_RIGHT:
                self.TURN = False
                self.TURN_RIGHT = False
                # print("turn right")

            elif self.TURN_LEFT:
                self.TURN = False
                self.TURN_LEFT = False
                # print("turn left")

            self.turn_array.append((self.phase,
                                    self.cur_lat,
                                    self.cur_lon))

            # we have passed the target waypoint,
            # current waypoint represents
            # "FROM which waypoint we are coming from"
            if self.next_wpt == len(self.route) - 1:  # no turns anymore as we are approaching end of the route
                self.FLIGHT_END = True

            else:
                self.cur_wpt = self.next_wpt
                self.next_wpt += 1

    def check_altitude(self):
        if self.alt < 0:
            self.LANDED = True

    def check_coordinates(self):
        """ This method is called if we are approaching the last waypoint on the route,
            to decide if we have already arrived and the flight is finished -
            we check the destination coordinates
        """
        cur_lat = round(degrees(self.cur_lat), 4)
        next_lat = round(degrees(self.route[self.next_wpt][0]), 4)
        cur_lon = round(degrees(self.cur_lon), 4)
        next_lon = round(degrees(self.route[self.next_wpt][1]), 4)

        if next_lat - 0.0001 < cur_lat <= next_lat + 0.0001 and next_lon - 0.0001 < cur_lon <= next_lon + 0.0001:
            self.LANDED = True

    def get_bank_angle(self):
        self.bank_angle = self.bank_angles[self.phase]

    def compute_turn_distance_and_rate(self):
        """
        Distance from turn start point to a waypoint, in relation to which turn is performed,
        d = tan (bearing_change / 2) * R;
        bearing_change - angle we are turning to (difference between brg2 and brg1)
        R - turn circle radius
        R = V**2 / (g * tan(phi)); phi - bank angle

        IMPORTANT
        In order to compensate for the wind conditions at the end of a turn
        use HEADING CHANGE ANGLE instead of course change angle when calculating turn geometry.
        This will result in a correction in a turn initiation distance, and at the end of a turn position will
        be very close to the intended one.

        """

        if self.next_wpt == len(self.route) - 1:  # no turns anymore as we are approaching end of the route
            self.FLIGHT_END = True
            return

        self.next_wp_course = degrees(self.course_array[self.next_wpt])
        self.current_track = degrees(self.track) % 360
        self.trk_diff = self.next_wp_course - self.current_track  # course change angle
        self.get_bank_angle()

        try:
            # self.r = (self.TAS ** 2) / (g * tan(self.bank_angle))  # [m]
            self.r = (self.GS ** 2) / (g * tan(self.bank_angle))  # [m]
        except:
            print("error in radius calculation")

        # calculate heading change angle to take into account the wind conditions in turn
        # self.current_heading = (degrees(self.hdg))%360
        # actual_wind_angle = self.wind_direction - radians(self.next_wp_course)  # line 5
        # WCA = asin((self.wind_speed * sin(actual_wind_angle)) / self.TAS)  # [rad]  line 6
        #
        # self.next_heading = ((self.next_wp_course + degrees(WCA))+360)%360
        # dif = self.next_heading - self.current_heading
        # self.heading_change_angle = abs(180. - abs(dif))

        self.course_change_angle = abs(180. - abs(self.trk_diff))  # is, in this case, turn inner angle

        # calculate distance from turn start to waypoint
        self.turn_initiation_distance = abs(self.r * tan(radians(self.trk_diff) / 2.))
        # second variant of turn init distance calculation
        # turn_initiation_distance_2 = abs(self.r / tan(radians(self.course_change_angle) / 2.))

        # to decide on turn direction - difference between bearings must be smaller then 180 degrees
        if abs(self.trk_diff) > 180:
            self.trk_diff = 180 - self.trk_diff

        # if difference between 2 legs is negative -> decrease the angle -> go left
        if self.trk_diff < 0:
            self.TURN_LEFT = True
            self.TURN_RIGHT = False
            # self.turn_rate = - self.TAS / self.r
            self.turn_rate = - self.GS / self.r

        # if difference between 2 legs is positive -> increase the angle -> go right
        elif self.trk_diff > 0:
            self.TURN_RIGHT = True
            self.TURN_LEFT = False
            # self.turn_rate = self.TAS / self.r
            self.turn_rate = self.GS / self.r

        # if the rate of turn in deg/sec is bigger than angle we are turning on - use following
        # if abs(degrees(self.turn_rate)) > abs(self.trk_diff):
        #     self.turn_rate = radians(self.trk_diff)

    def get_wind_correction_angle(self):
        """ Actual wind angle is the one of the angles, whichever is 90 deg or less
        formed by the wind and the track(course) line"""

        actual_wind_angle = self.wind_direction - self.track
        if not self.HAS:
            self.WCA = asin((self.wind_speed * sin(actual_wind_angle)) / self.TAS)  # [rad]
        else:
            try:
                self.WCA = asin((self.wind_speed * sin(actual_wind_angle)) / self.HAS)  # [rad]
            except:
                print('error in wind correction angle calculation')

    def update_heading(self):
        self.hdg = self.track + self.WCA  # [rad]

    def update_ground_speed(self):
        if not self.HAS:
            self.GS = sqrt(self.TAS ** 2 + self.wind_speed ** 2 -
                           2 * self.TAS * self.wind_speed * cos(self.wind_direction - self.hdg))  # [m/s]
        else:
            self.GS = sqrt(self.HAS ** 2 + self.wind_speed ** 2 -
                           2 * self.HAS * self.wind_speed * cos(self.wind_direction - self.hdg))  # [m/s]

    def set_rocd_old(self):
        if self.phase == "CLIMB":

            # if self.wtc == 'M':
            # levels = [0, 40, 105, 120, 200, 300, 400]
            # vertical_speeds_climb = [1000, 2000, 1000, 2000, 2000, 2000]
            for i in range(len(self.VS_climb_levels) - 2, -1, -1):
                if self.VS_climb_levels[i] * 100 * ft < self.alt <= self.VS_climb_levels[i + 1] * 100 * ft:
                    self.ROCD = self.VS_climb[i]

        if self.phase == "DESCENT":
            # if self.wtc == 'M':
            # levels = [0, 30, 110, 400]
            # vertical_speeds_descent = [-500, -1000, -2000]
            for i in range(len(self.VS_descent_levels) - 2, -1, -1):
                if self.VS_descent_levels[i] * 100 * ft < self.alt <= self.VS_descent_levels[i + 1] * 100 * ft:
                    self.ROCD = self.VS_descent[i]

    def set_rocd(self):
        if self.phase == "CLIMB":
            self.ROCD = self.VS_climb

        if self.phase == "DESCENT":
            self.ROCD = self.VS_descent

    def climb_phase(self):
        """ In Climb Phase Speed is controlled with FLC (Flight Level Control)mode:
               CAS is set to certain constant value within the flight levels range.
            Engines are set for the max power to ensure reaching the cruise FL as soon as possible
            and reduced power coefficient is applied to model reduced power climb.
        """

        # ------------------- CHECK SPEED SETTINGS -------------------------------------------------------------
        if self.alt < self.climb_crossover_altitude:

            if self.cas_climb > knot_to_mps(250):
                if self.FL < 100:
                    self.cas_target = knot_to_mps(250)
                else:
                    self.cas_target = self.cas_climb
            else:
                self.cas_target = self.cas_climb

            self.target_spd = mps_to_kts(self.cas_target)

            # check if target CAS has changed
            self.check_target_CAS()

            # check is acceleration/deceleration is necessary
            # self.check_CAS_for_acc_dec()
            self.check_TAS_for_acc_dec(below=True)

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

        else:
            # above crossover altitude fly at constant mach
            # ---------- UPDATE SPEEDS -------------------------------------------------------------------------
            # since transition to constant mach in climb happens at the point of crossover altitude
            # there is no need to accelerate or decelerate to required mach
            # just switch to mach is enough
            # self.Mach = get_mach(self.phase, self.apm, self.speed_profile)
            # self.Mach = self.apm.M_cl
            self.Mach = self.mach_climb
            self.TAS = mach_to_tas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.CAS = mach_to_cas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.target_spd = self.Mach
            # --------------------------------------------------------------------------------------------------

        # vertical speed settings based on the wtc and altitude range
        self.set_rocd()

    def cruise_phase(self):

        # ------------ CHECK SPEED SETTINGS --------------------------------------------------------------------
        # above crossover altitude fly at constant Mach
        # since mach number for cruise may be different from the climb mach number
        # it means that when aircraft is levelled it needs
        # to accelerate/decelerate to reach cruise mach number
        # therefore we need to check what is our current mach number
        # and if we need - initiate acceleration/deceleration

        # this check ensures that correct mach number value is used
        # as a target mach value in acceleration/deceleration phase
        # in cruise prior to descent

        if not self.reach_mach_for_descent:
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

    def descent_phase(self):

        # ---------------- CHECK SPEED SETTINGS ----------------------------------------------------------------

        if self.alt < self.descent_crossover_altitude:

            """ in descent aircraft is not allowed to exceed 250 CAS below FL100
                since aircraft needs time to decelerate and depending on mass it may take different time
                and decelerations needs to start at different altitudes
                therefore, it is necessary to check when aircraft is approaching FL100 that its target speed is
                250KCAS and it starts to decelerate on time 
            
                aircraft decelerates to 250 knots comply with restriction at and below FL100
            """
            if self.cas_descent > knot_to_mps(250):
                if self.alt < self.alt_start_decel:
                    self.cas_target = knot_to_mps(250)
                else:
                    self.cas_target = self.cas_descent
            else:
                self.cas_target = self.cas_descent

            self.target_spd = mps_to_kts(self.cas_target)

            # check if target CAS has changed
            self.check_target_CAS()

            # CHECK FOR CAS ACCELERATION/DECELERATION
            self.check_TAS_for_acc_dec(below=True)

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
            # self.Mach = get_mach(self.phase, self.apm, self.speed_profile)
            # if after deceleration in cruise this flag is not off - set it off
            # in any case deceleration from cruise mach to descent mach is done in cruise phase
            # and not in descent phase
            self.Mach = self.mach_descent
            self.DECELERATE_TAS = False

            self.TAS = mach_to_tas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.CAS = mach_to_cas(self.Mach, self.alt, self.delta_T, self.lapse_rate)
            self.target_spd = self.Mach
            # -------------------------------------------------------------------------------------------------

        self.set_rocd()

    def save_data(self):
        if self.results_file:
            self.results_file.close()
        return

    def display_map(self):

        # show flight on GMAPS
        self.gmap.addpath(self.path_points, "#FF0000")

        # show turns start and end
        for each in self.turn_array:
            self.gmap.addpoint(degrees(each[1]), degrees(each[2]), "#000000")

        # show flight phases start
        for each in self.flight_phases:
            self.gmap.addpoint(degrees(each[1]), degrees(each[2]), "#FFFFFF")

        if self.show_map:
            self.gmap.draw(self.path_to_save_xls + "flight_" + self.name + ".html")

        return

