import os
import pandas as pd
from datetime import datetime
from math import radians

from monte_carlo.simulate_true_flight import TrueFlightSimulation
from general.constants import ft
from monte_carlo.simulate_true_descent_backwards import TrueDescentSimulation
from general.apm import AircraftPerformanceModel

start_time = datetime.now()

# read a sheet from input data file
wtc = 'M'
sheet_num = 'MCsim_1'
path_to_input_data_file = os.path.dirname(os.path.realpath(__file__)) + '/inputs/MCsim_%s.xls' % wtc
input_data = pd.read_excel(path_to_input_data_file,
                           sheet_name=sheet_num)

test = False
show_map = False

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


class FlightMC(TrueFlightSimulation):

    def __init__(self, name, FPL, path_to_save_xls,
                 descent_distance, alt_start_decel,
                 log_file_path):

        super(FlightMC, self).__init__(name, FPL, path_to_save_xls,
                                       descent_distance, alt_start_decel)

        result = self.main()

        if not FPL['test']:
            self.write_states_to_txt()
            self.save_data()

        if FPL['show_map']:
            self.display_map()

        if not result:
            """ if something went wrong - mark flight in the log file"""
            log_file = open(log_file_path, 'a')
            log_file.write('%s  FAILED\n' % name)
            # print(name, " NOT OK")
            log_file.close()


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


def prepare_direct_route(route_id):
    path_to_route = os.path.join(os.path.dirname(__file__), 'fpls/FPL/')
    route_name = '%s_direct.csv' % routes[route_id][0]

    route_list = []

    with open(path_to_route+route_name, 'r') as route_file:
        route_data = route_file.readlines()
        for row in route_data:
            # do not include the last element with new line character
            row = row.rstrip().split(',')
            if row[0][0] != '#':
                route_list.append([radians(float(row[3])), radians(float(row[4]))])

    return route_list


for index, row in input_data.iterrows():
    path_to_save_flight_xls = '/media/julia/UserData/MC/MC_%s/flight/' % wtc
    path_to_save_descent_xls = os.path.dirname(os.path.realpath(__file__)) + '/maps/'

    if not test:
        print(index)
        # path_to_save_flight_xls = "/media/julia/UserData/MC/flight/"
        # path_to_save_descent_xls = "/media/julia/UserData/MC/descent/"

        file_name = sheet_num + "_" + wtc + '_' + str(index)

        bank = {'CLIMB': radians(row['BANK_Climb']),
                'CRUISE': radians(row['BANK_Cruise']),
                'DESCENT': radians(row['BANK_Descent'])}

        leveloff = {'CLIMB': [[row['LVLOFF_FL_Climb']], [row['LVLOFF_DUR_Climb']]],
                    'DESCENT': [[row['LVLOFF_FL_Descent']], [row['LVLOFF_DUR_Descent']]]}

        actype = row['ACFT']
        apm = AircraftPerformanceModel(actype)
        cruising_FL = 350  # FL

        temp_diff = row['TEMP_OFFSET']
        lapse_rate = row['LAPSE_RATE'] / 1000

        speeds = {'CLIMB': [row['CAS_Climb'], row['M_Climb']],
                  'CRUISE': [row['M_Cruise']],
                  'DESCENT': [row['CAS_Descent'], row['M_Descent']]}

        vertical_speed_climb = row['ROC']
        vertical_speeds_descent = row['ROD']

        vs = {'Climb': vertical_speed_climb,
              'DescentPrediction': vertical_speeds_descent}

        wind_comp = {'u': [row['U_a'], row['U_b'], row['U_c']],
                     'v': [row['V_a'], row['V_b'], row['V_c']]}

        # zero wins
        # wind_comp = {"u": [0, 0, 0],
        #              "v": [0, 0, 0]}

        route_id = row['N_FPL']

        instruction = {'type': row['ATC_instr'],
                       'value': row['ATC_value'],
                       'dur': row['ATC_dur']}

        flight_plan = {'actype': actype,
                       'apm': apm,
                       'route': prepare_direct_route(route_id) if instruction['type'] == 'direct' else prepare_route(route_id),
                       'route_dist': routes[route_id][2] if instruction['type'] == 'direct' else routes[route_id][1],
                       'route_id': route_id,
                       'FL': cruising_FL,
                       'wtc': wtc,
                       'VS': vs,
                       'bank': bank,
                       'speeds': speeds,
                       'leveloff': leveloff,
                       'temp_diff': temp_diff,
                       # 'wind': wind,
                       'lapse_rate': lapse_rate,
                       # 'pygmap_center': pygmap_center,
                       'show_map': show_map,
                       'test': test,
                       'wind_comp': wind_comp,
                       'instruction': instruction}

        descent_name = 'descent_backwards_mc_test'
        descent = TrueDescentSimulation(descent_name, flight_plan,
                                        path_to_save_descent_xls, 0, [0])

        descent_distance = descent.distance_flown

        """ !!! add 1.5 FL to alt_start_decel in descent 
        to ensure that aircraft can decelerate 
        if it has a vertical constraint at FL100 """
        alt_start_decel = [each + 150 * ft for each in descent.descent_deccel_alt]

        log_file_path = os.path.dirname(os.path.realpath(__file__)) + '/res_log_%s.txt' % wtc

        a = FlightMC(file_name, flight_plan, path_to_save_flight_xls,
                     descent_distance, alt_start_decel, log_file_path)

print('code ran for', datetime.now()-start_time)
print('done')
