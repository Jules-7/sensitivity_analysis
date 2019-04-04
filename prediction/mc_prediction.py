import os
from datetime import datetime
from math import radians
import numpy as np
import pandas as pd

from prediction.predict_flight import FlightPrediction
from prediction.predict_descent_backwards import DescentPrediction
from monte_carlo.dists import WIND, rocd, cas, mach, air_temp_fl50
from monte_carlo.dists import bank as bank_data
from monte_carlo.dists import lapse_rate as lapse_rate_data

from general.constants import ft
from total_energy_model.conversion import mps_to_kts
from general.apm import AircraftPerformanceModel

from prediction.errors_computation import ate_cte_error, alt_error

sim_version = 1  # 1, 2, 3, 4
wtc = 'M'
test = False
show_map = True

start_time = datetime.now()

#              name | dist in nm regular route | dist in nm direct route
routes = {0: ['EBAW-LEMG', 968.0, 963.2],
          1: ['EDDH-LEMD', 979.6, 971.5],
          2: ['EFRO-EGDD', 1245.8, 1239.4],
          3: ['EGLL-LEMD', 703.1, 697.9],
          4: ['EGNT-BIAR', 847.8, 813.9],
          5: ['EGNT-LDPL', 869.6, 867.9],
          6: ['EGNT-LERS', 854.3, 849.9],
          7: ['EHAM-LIRF', 716.0, 713.1],
          8: ['EPWR-EGNT', 731.4, 723.6],
          9: ['LIRF-EFHK', 1240.5, 1234.6]}

path_to_flights = '/media/julia/UserData/MC/true_flights/%s/' % wtc
path_to_save_errors = '/media/julia/UserData/MC/lat_errors_%s/sim)res_%s' % (wtc, sim_version)
path_to_input_data_file = os.path.dirname(os.path.realpath(__file__)) + '/MCsim_%s.xls' % wtc
path_to_save_xls = os.path.dirname(os.path.realpath(__file__)) + '/maps/'

input_sheet_0 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_0')
input_sheet_1 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_1')
input_sheet_2 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_2')
input_sheet_3 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_3')
input_sheet_4 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_4')
input_sheet_5 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_5')
input_sheet_6 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_6')
input_sheet_7 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_7')
input_sheet_8 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_8')
input_sheet_9 = pd.read_excel(path_to_input_data_file, sheet_name='MCsim_9')

sheets = {0: input_sheet_0,
          1: input_sheet_1,
          2: input_sheet_2,
          3: input_sheet_3,
          4: input_sheet_4,
          5: input_sheet_5,
          6: input_sheet_6,
          7: input_sheet_7,
          8: input_sheet_8,
          9: input_sheet_9}


def prepare_route(route_id):
    path_to_route = os.path.join(os.path.dirname(__file__), 'fpls/FPL/')
    route_name = '%s.csv' % routes[route_id][0]

    route_list = []

    with open(path_to_route+route_name, 'r') as route_file:
        route_data = route_file.readlines()
        for row in route_data:
            # do not include the last element with new line character
            row = row.rstrip().split(',')
            if row[0][0] != '#':
                route_list.append([radians(float(row[3])), radians(float(row[4]))])

    return route_list


def sample_wind_u():
    """ wind components coefficients are
    sampled from a uniform distribution"""
    wind_u = {'a': [], 'b': [], 'c': []}

    for key, value in wind_u.items():
        # 'a' and 'b' for uniform distribution
        low, high = WIND['u'][key]

        np.random.seed(111)
        # np.random.seed(7474 + wind_u_counter)
        sample = np.random.uniform(low, high, size=1)

        wind_u[key] = sample

    return wind_u


def sample_wind_v():
    """ wind components coefficients are
    sampled from a uniform distribution"""
    wind_v = {'a': [], 'b': [], 'c': []}

    for key, value in wind_v.items():
        # 'a' and 'b' for uniform distribution
        low, high = WIND['v'][key]

        np.random.seed(222)
        # np.random.seed(3232 + wind_v_counter)
        sample = np.random.uniform(low, high, size=1)

        wind_v[key] = sample

    return wind_v


def prepare_prediction_inputs(route_id, acft_name):
    # standard bank value in all phases
    if sim_version == 4:
        bank_climb_dist = bank_data[wtc]['Climb'][0]
        bank_climb_param = bank_data[wtc]['Climb'][1]

        bank_cruise_dist = bank_data[wtc]['Cruise'][0]
        bank_cruise_param = bank_data[wtc]['Cruise'][1]

        bank_descent_dist = bank_data[wtc]['DescentPrediction'][0]
        bank_descent_param = bank_data[wtc]['DescentPrediction'][1]

        bank_climb = bank_climb_dist(*bank_climb_param).mean()
        bank_cruise = bank_cruise_dist(*bank_cruise_param).mean()
        bank_descent = bank_descent_dist(*bank_descent_param).mean()

        bank = {'CLIMB': radians(bank_climb),
                'CRUISE': radians(bank_cruise),
                'DESCENT': radians(bank_descent)}
    else:
        bank = {'CLIMB': radians(25), 'CRUISE': radians(25), 'DESCENT': radians(25)}

    # no level-offs
    leveloff = {'CLIMB': [[0], [0]], 'DESCENT': [[0], [0]]}

    # actype = acft_types[wtc][0][route_id]
    apm = AircraftPerformanceModel(acft_name)
    cruising_FL = 350

    if sim_version == 4:
        air_temp_dist = air_temp_fl50[0]
        air_temp_param = air_temp_fl50[1]

        lapse_rate_dist = lapse_rate_data[0]
        lapse_rate_param = lapse_rate_data[1]

        air_temp = air_temp_dist(*air_temp_param).mean()
        lapse_rate = lapse_rate_dist(*lapse_rate_param).mean() / 1000

        # compute what would be the temperature at FL0
        # using the lapse rate
        delta_h = 50 * 100 * 0.3048  # meters
        delta_temp = delta_h * lapse_rate
        temp_fl0 = air_temp - delta_temp

        # compute what is the temperature offset from standard
        # temperature at FL0
        temp_isa = 288.15
        temp_diff = temp_fl0 - temp_isa
    else:
        # ICAO ISA
        temp_diff = 0
        lapse_rate = -6.5 / 1000

    if sim_version <= 2:
        # version 1,2
        # speed values are taken from an aircraft performance model
        # apm converts speed to mps,
        # whereas further in the Flight module
        # it is converted to mps as well
        # thus, convert apm speeds back to knots
        speeds = {'CLIMB': [mps_to_kts(apm.V_cl_2), apm.M_cl],
                  'CRUISE': [apm.M_cr],
                  'DESCENT': [mps_to_kts(apm.V_des_2), apm.M_des]}
        # speeds = {'CLIMB': [311., 0.78], 'CRUISE': [0.75], 'DESCENT': [280., 0.76]}

    else:
        # version 3
        # speed values are assigned as a mean value of a respective
        # WTC distribution
        cas_climb_dist = cas[wtc]['Climb'][0]
        cas_climb_param = cas[wtc]['Climb'][1]

        cas_descent_dist = cas[wtc]['DescentPrediction'][0]
        cas_descent_param = cas[wtc]['DescentPrediction'][1]

        mach_climb_dist = mach[wtc]['Climb'][0]
        mach_climb_param = mach[wtc]['Climb'][1]

        mach_cruise_dist = mach[wtc]['Cruise'][0]
        mach_cruise_param = mach[wtc]['Cruise'][1]

        mach_descent_dist = mach[wtc]['DescentPrediction'][0]
        mach_descent_param = mach[wtc]['DescentPrediction'][1]

        cas_climb = cas_climb_dist(*cas_climb_param).mean()
        cas_descent = cas_descent_dist(*cas_descent_param).mean()

        mach_climb = mach_climb_dist(*mach_climb_param).mean()
        mach_cruise = mach_cruise_dist(*mach_cruise_param).mean()
        mach_descent = mach_descent_dist(*mach_descent_param).mean()

        speeds = {'CLIMB': [cas_climb, mach_climb],
                  'CRUISE': [mach_cruise],
                  'DESCENT': [cas_descent, mach_descent]}

    # wind is sampled form the wind components distributions
    u = sample_wind_u()
    v = sample_wind_v()

    wind_comp = {'u': [u['a'][0], u['b'][0], u['c'][0]],
                 'v': [v['a'][0], v['b'][0], v['c'][0]]}

    # zero wind
    # wind_comp = {"u": [0, 0, 0],
    #              "v": [0, 0, 0]}

    if sim_version == 1:
        # Version 1
        # assign fixed vertical speed in climb and descent for the prediction
        vertical_speed_climb = 2500
        vertical_speeds_descent = -1500

    else:
        # Version 2
        # assign fixed vertical speed in climb and descent for the prediction
        # based on the mean vertical speed for a respective WTC
        vertical_speed_climb_dist = rocd[wtc]['Climb'][999][0]
        vertical_speed_climb_param = rocd[wtc]['Climb'][999][1]

        vertical_speeds_descent_dist = rocd[wtc]['DescentPrediction'][999][0]
        vertical_speeds_descent_param = rocd[wtc]['DescentPrediction'][999][1]

        vertical_speed_climb = vertical_speed_climb_dist(*vertical_speed_climb_param).mean()
        vertical_speeds_descent = vertical_speeds_descent_dist(*vertical_speeds_descent_param).mean()

    vs = {'Climb': vertical_speed_climb,
          'DescentPrediction': vertical_speeds_descent}

    # no ATC instructions are predicted
    instruction = {'type': 0, 'value': 0, 'dur': 0}

    flight_plan = {'actype': acft_name,
                   'apm': apm,
                   'route': prepare_route(route_id),
                   'route_dist': routes[route_id][1],
                   'route_id': route_id,
                   'FL': cruising_FL,
                   'wtc': wtc,
                   'VS': vs,
                   'bank': bank,
                   'speeds': speeds,
                   'leveloff': leveloff,
                   'temp_diff': temp_diff,
                   'lapse_rate': lapse_rate,
                   'show_map': show_map,
                   'test': test,
                   'wind_comp': wind_comp,
                   'instruction': instruction}

    descent = DescentPrediction('descent_prediction', flight_plan, path_to_save_xls)

    descent_distance = descent.distance_flown

    """ !!! add 1.5 FL to alt_start_decel in descent
    to ensure that aircraft can decelerate
    if it has a vertical constraint at FL100 """
    alt_start_decel = [each + 150 * ft for each in descent.descent_deccel_alt]

    return flight_plan, descent_distance, alt_start_decel


# loop over folders with flights
# for route_id in [9]:
for route_id in range(0, 10):

    # ################ prepare inputs for a prediction ###############################
    # if code below comes before the true flight data is read,
    # it uses the route_id value to choose the aircraft type,
    # which is different from the aircraft type of the true flight
    # this might not cause significant differences in the obtained results,
    # but must be checked.
    # flight_plan, descent_distance, alt_start_decel = prepare_prediction_inputs(route_id)

    counter = 0
    # ############## loop over each real flight #############################
    for flight in os.listdir(path_to_flights + '%s' % route_id):

        if '__' not in flight:  # disregard failed flights
            counter += 1
            print(counter, flight)

            # read input data used for the true flight generation
            # and extract the aircraft type to be used in predictions
            sheet_num = int(flight.split('_')[1])
            flight_id = int(flight.split('_')[3])

            # extract inputs for this flight
            row_inputs = sheets[sheet_num].loc[flight_id, :]
            acft_name = row_inputs['ACFT']

            # ############### prepare inputs for a prediction ###############################
            flight_plan, descent_distance, alt_start_decel = prepare_prediction_inputs(route_id, acft_name)

            # prepare files to record errors
            try:
                os.remove(path_to_save_errors + '%s/5/' % route_id + '%s_lat_5.txt' % flight[:-4])
                os.remove(path_to_save_errors + '%s/10/' % route_id + '%s_lat_10.txt' % flight[:-4])
                os.remove(path_to_save_errors + '%s/15/' % route_id + '%s_lat_15.txt' % flight[:-4])
                os.remove(path_to_save_errors + '%s/20/' % route_id + '%s_lat_20.txt' % flight[:-4])

            except OSError:
                pass

            lat_5_file = open(path_to_save_errors + '%s/5/' % route_id + '%s_lat_5.txt' % flight[:-4], 'a')
            lat_10_file = open(path_to_save_errors + '%s/10/' % route_id + '%s_lat_10.txt' % flight[:-4], 'a')
            lat_15_file = open(path_to_save_errors + '%s/15/' % route_id + '%s_lat_15.txt' % flight[:-4], 'a')
            lat_20_file = open(path_to_save_errors + '%s/20/' % route_id + '%s_lat_20.txt' % flight[:-4], 'a')

            lat_files = {5: lat_5_file,
                         10: lat_10_file,
                         15: lat_15_file,
                         20: lat_20_file}

            # extract information about the real flight
            flight_data = pd.read_csv(path_to_flights + '%s/' % route_id + flight)
            time_t = flight_data['TIME_s'].values.astype(int)
            x_t = flight_data['X_nm'].values
            y_t = flight_data['Y_nm'].values
            lat_t = flight_data['LAT_DD'].values
            lon_t = flight_data['LONG_DD'].values
            alt_t = flight_data['ALT_m'].values
            fl_t = flight_data['FL'].values.astype(int)
            u_t = flight_data['U_kt'].values
            v_t = flight_data['V_kt'].values
            phase_t = flight_data['PHASE'].values

            # reset the last time instance to be divisible by minute
            time_t[-1] = time_t[-2] + 60

            for i, time_val in enumerate(time_t.tolist()[:-5]):

                real_lat, real_lon, real_alt, real_fl = lat_t[i], lon_t[i], alt_t[i], fl_t[i]
                real_phase = phase_t[i]

                pred = FlightPrediction(flight_plan, descent_distance, alt_start_decel,
                                        real_lat, real_lon, real_alt, real_fl, real_phase)

                if pred.init_success:
                    predictions = pred.main()

                    for lat in [5, 10, 15, 20]:

                        try:
                            # extract positions for this time from real flight
                            time_id = np.where(time_t == time_val + lat * 60)[0][0]
                            x_t_lat = x_t[time_id]
                            y_t_lat = y_t[time_id]
                            alt_t_lat = alt_t[time_id]
                            phase_t_lat = phase_t[time_id]
                            u_t_lat = np.average(u_t[i:time_id])
                            v_t_lat = np.average(v_t[i:time_id])

                            # predicted values
                            x_p_lat = predictions[lat][0]
                            y_p_lat = predictions[lat][1]
                            alt_p_lat = predictions[lat][2]
                            trk_p_lat = predictions[lat][4]
                            phase_p_lat = predictions[lat][5]
                            u_p_lat = predictions[lat][6]
                            v_p_lat = predictions[lat][7]

                            # compute errors
                            ate, cte = ate_cte_error(x_p_lat, y_p_lat, trk_p_lat,
                                                     x_t_lat, y_t_lat)
                            alte = alt_error(alt_p_lat, alt_t_lat) * ft

                            if abs(alte) < 10:
                                alte = 0

                            data_to_write = ','.join(['%.2f'%ate, '%.2f'%cte, '%.2f'%alte, str(u_t_lat),
                                                      str(v_t_lat), phase_t_lat, phase_p_lat, '\n'])

                            lat_files[lat].write(data_to_write)

                        except IndexError:
                            pass

            for key, val in lat_files.items():
                val.close()

print('code ran for', datetime.now() - start_time)
print('done')
