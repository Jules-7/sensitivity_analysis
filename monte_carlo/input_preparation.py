""" GENERATE INPUTS FOR TRUE FLIGHTS SIMULATION


!!! DO NOT FORGER TO USE RANDOM STATE WHEN SAMPLING FROM DISTRIBUTIONS WITH SCIPY
AND USE NP.RANDOM.SEED(RANDOM_SEED) !!!PRIOR!!! TO EVERY NP.RANDOM FUNCTION
"""

import os
import csv
import numpy as np
import scipy.stats as st
import xlwt

from monte_carlo.dists import (rand_seed, BANK, TEMP_LEVELOFF, ACFT_TYPE, AIR_TEMP,
                               SPEED_CAS, SPEED_MACH, ROCD, WIND, ATC, random_seeds_for_resample, random_seeds)


class InputPreparation(object):

    def __init__(self, n, wtc, path):

        self.N = n  # number of flights
        self.wtc = wtc
        self.path_to_write = path
        self.rand_seed_wtc = rand_seed[self.wtc]

        # np.random.seed(self.rand_seed)

        self.sheet_counter = 0

        self.prepare_xls_file()

        # containers for values
        self.bank = {'Climb': [], 'Cruise': [], 'DescentPrediction': []}

        self.temp_leveloff = {'Climb': [], 'DescentPrediction': []}
        self.temp_leveloff_dur = {'Climb': [], 'DescentPrediction': []}
        self.temp_leveloff_fl = {'Climb': [], 'DescentPrediction': []}

        self.aircraft_types = {'ACFT': []}

        self.air_temp = {'TEMP_OFFSET': [],
                         'LAPSE_RATE': []}

        self.speed_cas = {'Climb': [], 'DescentPrediction': []}
        self.speed_mach = {'Climb': [], 'Cruise': [], 'DescentPrediction': []}

        # if self.wtc in ('H', 'M'):
        #     self.roc = {0: [], 40: [], 105: [], 120: [], 200: [], 300: []}
        #     self.rod = {0: [], 30: [], 110: []}
        #
        # elif self.wtc == 'L':
        #     self.roc = {0: [], 30: [], 100: [], 200: [], 300: []}
        #     self.rod = {0: [], 30: []}

        self.roc = {'ROC': []}
        self.rod = {'ROD': []}

        self.wind_u = {'a': [], 'b': [], 'c': []}
        self.wind_v = {'a': [], 'b': [], 'c': []}

        self.fpl = {'FPL': []}

        self.atc = {'ATC_instr': [], 'ATC_value': [], 'ATC_dur': []}

        # perform inputs sampling
        # self.sample_mass()
        self.sample_bank()
        self.sample_level()
        # self.sample_level_1()
        self.sample_acft_type()
        self.sample_temp()
        self.sample_speed()
        self.sample_roc()
        self.sample_rod()

        # separately sample wind components (u and v)
        self.sample_wind_u()
        self.sample_wind_v()

        self.sample_fpl()
        self.sample_atc()

        self.save_xls_file()

    def prepare_xls_file(self):
        self.workbook = xlwt.Workbook()
        self.sheet = self.workbook.add_sheet('MCsim_%s' % self.sheet_counter)
        self.rows_count = 0
        self.columns_count = 0
        return

    def save_xls_file(self):
        xls_file = self.path_to_write + 'MCsim_' + self.wtc + '.xls'
        self.workbook.save(xls_file)

    def sample_mass(self):
        # self.mass = st.norm(60000, 5000).rvs(self.sample_size)
        pass

    def sample_bank(self):
        random_seed_climb = random_seeds['bank_climb']
        random_seed_cruise = random_seeds['bank_cruise']
        random_seed_descent = random_seeds['bank_descent']

        for key, value in self.bank.items():

            # distribution type and parameters per phase
            dist, param = BANK['val'][self.wtc][key]

            # sample from distribution
            if key == 'DescentPrediction':
                # use a different random seed for descent,
                # otherwise values for bank in climb and descent are almost the same
                # since their distributions are very similar
                sample = dist.rvs(*param, size=self.N,
                                  random_state=self.rand_seed_wtc + random_seed_descent)
            elif key == 'Climb':
                sample = dist.rvs(*param, size=self.N,
                                  random_state=self.rand_seed_wtc + random_seed_climb)

            elif key == 'Cruise':
                sample = dist.rvs(*param, size=self.N,
                                  random_state=self.rand_seed_wtc + random_seed_cruise)

            # check that all values are within limits
            self.bank[key] = self.check_sample(sample, dist, param, 'bank')

        # write results to a file
        self.write_to_file(self.bank, 'BANK')

    def sample_level_1(self):

        """old version of code. New version is def sample_level()"""

        # percentage of flights with level-off in climb and descent
        climb_freq = TEMP_LEVELOFF['freq'][self.wtc]['Climb']
        descent_freq = TEMP_LEVELOFF['freq'][self.wtc]['DescentPrediction']

        # number of flights with level-off in climb and descent
        climb_N = int(self.N * (climb_freq / 100))
        descent_N = int(self.N * (descent_freq / 100))

        # randomly select which flights (flight ids) will have level-off in climb and descent
        # use replace=False, otherwise same ids can be sampled multiple times
        np.random.seed(self.rand_seed)
        climb_leveoff_ids = np.random.choice(range(0, self.N, 1), climb_N, replace=False)
        np.random.seed(self.rand_seed)
        descent_leveoff_ids = np.random.choice(range(0, self.N, 1), descent_N, replace=False)

        # distribution and parameters of level-offs duration in climb and descent
        climb_dist, climb_param = TEMP_LEVELOFF['dur'][self.wtc]['Climb']
        descent_dist, descent_param = TEMP_LEVELOFF['dur'][self.wtc]['DescentPrediction']

        # range and frequencies of level-offs FL in climb and descent
        climb_range, climb_freq = TEMP_LEVELOFF['fl'][self.wtc]['Climb']
        descent_range, descent_freq = TEMP_LEVELOFF['fl'][self.wtc]['DescentPrediction']

        # sample distributions
        leveloff_dur_in_climb = climb_dist.rvs(*climb_param, size=self.N, random_state=self.rand_seed)
        leveloff_dur_in_descent = descent_dist.rvs(*descent_param, size=self.N, random_state=self.rand_seed)

        # sample range and frequency
        np.random.seed(self.rand_seed)
        leveloff_fl_in_climb = np.random.choice(climb_range, p=climb_freq, size=self.N, replace=True)
        np.random.seed(self.rand_seed)
        leveloff_fl_in_descent = np.random.choice(descent_range, p=descent_freq, size=self.N, replace=True)
        print(climb_freq)

        # check that all values are within limits
        leveloff_dur_in_climb = self.check_sample(leveloff_dur_in_climb, climb_dist, climb_param, 'leveloff dur')
        leveloff_dur_in_descent = self.check_sample(leveloff_dur_in_descent, descent_dist, descent_param, 'leveloff dur')

        leveloff_fl_in_climb = self.check_sample(leveloff_fl_in_climb, climb_range, climb_freq, 'leveloff fl')
        leveloff_fl_in_descent = self.check_sample(leveloff_fl_in_descent, descent_range, descent_freq, 'leveloff fl')

        # based on randomly generated flight ids keep only those values
        # whose flight id corresponds to value id + 1 in the list
        # replace the rest with 0
        for i in range(0, self.N, 1):
            if i not in climb_leveoff_ids:
                leveloff_dur_in_climb[i] = 0
                leveloff_fl_in_climb[i] = 0

        for j in range(0, self.N, 1):
            if j not in descent_leveoff_ids:
                leveloff_dur_in_descent[j] = 0
                leveloff_fl_in_descent[j] = 0

        # add values to a dictionary
        self.temp_leveloff_dur['Climb'] = leveloff_dur_in_climb
        self.temp_leveloff_dur['DescentPrediction'] = leveloff_dur_in_descent

        self.temp_leveloff_fl['Climb'] = leveloff_fl_in_climb
        self.temp_leveloff_fl['DescentPrediction'] = leveloff_fl_in_descent

        # write values into a file
        self.write_to_file(self.temp_leveloff_fl, 'LVLOFF_FL')
        self.write_to_file(self.temp_leveloff_dur, 'LVLOFF_DUR')

    def sample_level(self):
        random_seed_freq = random_seeds['temp_leveloff_freq']
        random_seed_dur = random_seeds['temp_leveloff_dur']
        random_seed_fl = random_seeds['temp_leveloff_fl']

        freq_counter = 0
        dur_counter = 0
        fl_counter = 0

        # obtain ids of flights that will have level-offs
        for key, value in self.temp_leveloff.items():

            # percentage of flights with level-off in climb and descent
            val_freq = TEMP_LEVELOFF['freq'][self.wtc][key]

            # number of flights with level-off in climb and descent
            sample_N = int(self.N * (val_freq / 100))

            # randomly select which flights (flight ids) will have level-off in climb and descent
            # use replace=False, otherwise same ids can be sampled multiple times
            np.random.seed(self.rand_seed_wtc + random_seed_freq + freq_counter)
            self.temp_leveloff[key] = np.random.choice(range(0, self.N, 1), sample_N, replace=False)
            freq_counter += 1

        # first determine level-off durations
        for key, value in self.temp_leveloff_dur.items():

            # distribution and parameters of level-offs duration in climb and descent
            dist, param = TEMP_LEVELOFF['dur'][self.wtc][key]

            # sample distributions
            sample = dist.rvs(*param, size=self.N, random_state=self.rand_seed_wtc + random_seed_dur + dur_counter)

            # check that all values are within limits
            sample = self.check_sample(sample, dist, param, 'leveloff dur')

            # based on randomly generated flight ids keep only those values
            # whose flight id corresponds to value id + 1 in the list
            # replace the rest with 0
            for i in range(0, self.N, 1):
                if i not in self.temp_leveloff[key]:
                    sample[i] = 0

            # add values to a dictionary
            self.temp_leveloff_dur[key] = sample

            dur_counter += 1

        for key, value in self.temp_leveloff_fl.items():

            # range and frequencies of level-offs FL in climb and descent
            val_range, val_freq = TEMP_LEVELOFF['fl'][self.wtc][key]

            # sample range and frequency
            np.random.seed(self.rand_seed_wtc + random_seed_fl + fl_counter)
            sample = np.random.choice(val_range, p=val_freq, size=self.N, replace=True)

            # check that all values are within limits
            sample = self.check_sample(sample, val_range, val_freq, 'leveloff fl')

            for j in range(0, self.N, 1):
                if j not in self.temp_leveloff[key]:
                    sample[j] = 0

            self.temp_leveloff_fl[key] = sample

            fl_counter += 1

        # write values into a file
        self.write_to_file(self.temp_leveloff_fl, 'LVLOFF_FL')
        self.write_to_file(self.temp_leveloff_dur, 'LVLOFF_DUR')

    def sample_acft_type(self):
        random_seed_acft_type = random_seeds['acft_type']

        # range and frequencies of aircraft types
        types_range, types_freq = ACFT_TYPE[self.wtc]
        np.random.seed(self.rand_seed_wtc + random_seed_acft_type)
        aircraft_types = np.random.choice(types_range, p=types_freq, size=self.N, replace=True)

        self.aircraft_types['ACFT'] = aircraft_types.tolist()

        self.write_to_file(self.aircraft_types, '')

    def sample_temp(self):
        random_seed_temp = random_seeds['temp']
        random_seed_lapse = random_seeds['lapse']

        # distribution type and parameters
        temp_fl50_dist, temp_fl50_param = AIR_TEMP['temp fl50']
        lapse_rate_dist, lapse_rate_param = AIR_TEMP['lapse rate']

        # sample from distributions
        temp_fl50 = temp_fl50_dist.rvs(*temp_fl50_param, size=self.N,
                                       random_state=self.rand_seed_wtc + random_seed_temp)
        lapse_rate = lapse_rate_dist.rvs(*lapse_rate_param, size=self.N,
                                         random_state=self.rand_seed_wtc + random_seed_lapse)

        # check that all values are within limits
        temp_fl50 = np.array(self.check_sample(temp_fl50, temp_fl50_dist, temp_fl50_param, 'temp'))
        lapse_rate = np.array(self.check_sample(lapse_rate, lapse_rate_dist, lapse_rate_param, 'lapse'))

        # compute what would be the temperature at FL0
        # using the lapse rate
        delta_h = 50 * 100 * 0.3048  # meters
        lapse_rate_in_kpm = lapse_rate / 1000  # kelvin per meter
        delta_temp = delta_h * lapse_rate_in_kpm
        temp_fl0 = temp_fl50 - delta_temp

        # compute what is the temperature offset from standard
        # temperature at FL0
        temp_stndt = 288.15
        temp_offset = temp_fl0 - temp_stndt

        self.air_temp['TEMP_OFFSET'] = temp_offset.tolist()
        self.air_temp['LAPSE_RATE'] = lapse_rate.tolist()

        self.write_to_file(self.air_temp, "")

    def sample_speed(self):
        random_seed_cas_climb = random_seeds['cas_climb']
        random_seed_cas_descent = random_seeds['cas_descent']
        random_seed_mach_climb = random_seeds['mach_climb']
        random_seed_mach_cruise = random_seeds['mach_cruise']
        random_seed_mach_descent = random_seeds['mach_descent']

        # first CAS
        for key, value in self.speed_cas.items():
            # distribution and parameters per flight phase
            distr, param = SPEED_CAS['val'][self.wtc][key]

            if key == 'Climb':
                # sample from distributions
                sample = distr.rvs(*param, size=self.N,
                                   random_state=self.rand_seed_wtc + random_seed_cas_climb)

            elif key == 'DescentPrediction':
                # sample from distributions
                sample = distr.rvs(*param, size=self.N,
                                   random_state=self.rand_seed_wtc + random_seed_cas_descent)

            # check that all values are within limits
            self.speed_cas[key] = self.check_sample(sample, distr, param, 'cas')

        # then Mach number
        for key, value in self.speed_mach.items():
            # distribution and parameters
            distr, param = SPEED_MACH['val'][self.wtc][key]

            if key == 'Climb':
                # sample from distribution
                sample = distr.rvs(*param, size=self.N,
                                   random_state=self.rand_seed_wtc + random_seed_mach_climb)

            elif key == 'Cruise':
                # sample from distribution
                sample = distr.rvs(*param, size=self.N,
                                   random_state=self.rand_seed_wtc + random_seed_mach_cruise)

            elif key == 'DescentPrediction':
                # sample from distribution
                sample = distr.rvs(*param, size=self.N,
                                   random_state=self.rand_seed_wtc + random_seed_mach_descent)

            # check that all values are within limits
            self.speed_mach[key] = self.check_sample(sample, distr, param, 'mach')

        self.write_to_file(self.speed_cas, 'CAS')
        self.write_to_file(self.speed_mach, 'M')

    def sample_roc_old(self):
        random_seed_roc = random_seeds['roc']

        roc_counter = 0

        for key, value in self.roc.items():

            # distribution and parameters per altitude range
            distr, param = ROCD['val'][self.wtc]['Climb'][key]

            # sample from distributions
            sample = distr.rvs(*param, size=self.N,
                               random_state=self.rand_seed_wtc + random_seed_roc + roc_counter)

            # check that all values are within limits
            self.roc[key] = self.check_sample(sample, distr, param, 'roc')

            roc_counter += 12

        self.write_to_file(self.roc, 'ROC')

    def sample_roc(self):

        random_seed_roc = random_seeds['roc']

        # distribution and parameters per altitude range
        distr, param = ROCD['val'][self.wtc]['Climb'][999]

        # sample from distributions
        sample = distr.rvs(*param, size=self.N, random_state=self.rand_seed_wtc + random_seed_roc)

        # check that all values are within limits
        self.roc['ROC'] = self.check_sample(sample, distr, param, 'roc')

        self.write_to_file(self.roc, '')

    def sample_rod_old(self):
        random_seed_rod = random_seeds['rod']

        rod_counter = 0

        for key, value in self.rod.items():

            # distribution and parameters per altitude range
            distr, param = ROCD['val'][self.wtc]['DescentPrediction'][key]

            # sample from distributions
            sample = distr.rvs(*param, size=self.N,
                               random_state=self.rand_seed_wtc + random_seed_rod + rod_counter)

            # check that all values are within limits
            self.rod[key] = self.check_sample(sample, distr, param, 'rod')

            rod_counter += 3

        self.write_to_file(self.rod, 'ROD')

    def sample_rod(self):

        random_seed_rod = random_seeds['rod']

        # distribution and parameters per altitude range
        distr, param = ROCD['val'][self.wtc]['DescentPrediction'][999]

        # sample from distributions
        sample = distr.rvs(*param, size=self.N, random_state=self.rand_seed_wtc + random_seed_rod)

        # check that all values are within limits
        self.rod['ROD'] = self.check_sample(sample, distr, param, 'rod')

        self.write_to_file(self.rod, '')

    def sample_wind_u(self):
        """wind components coefficients are
        sampled from a uniform distribution"""

        random_seed_u = random_seeds['wind_u']

        u_counter = 0

        for key, value in self.wind_u.items():
            # 'a' and 'b' for uniform distribution
            low, high = WIND['u'][key]

            # sample from distribution
            np.random.seed(self.rand_seed_wtc + random_seed_u + u_counter)
            sample = np.random.uniform(low, high, size=self.N)

            self.wind_u[key] = sample

            u_counter += 7

        self.write_to_file(self.wind_u, "U")

    def sample_wind_v(self):
        """wind components coefficients are
        sampled from a uniform distribution"""

        random_seed_v = random_seeds['wind_v']

        v_counter = 0

        for key, value in self.wind_v.items():
            # 'a' and 'b' for uniform distribution
            low, high = WIND['v'][key]

            # sample from distribution
            np.random.seed(self.rand_seed_wtc + random_seed_v + v_counter)
            sample = np.random.uniform(low, high, size=self.N)

            self.wind_v[key] = sample

            v_counter += 33

        self.write_to_file(self.wind_v, "V")

    def sample_atc(self):

        random_seed_atc = random_seeds['atc']

        # percentage of flights without ATC instruction
        instr_names = ATC['instr'][0]
        instr_probs = ATC['instr'][1]
        instr_numbs = list(range(0, len(instr_names)))

        # sample each instruction or no instruction according to its probability
        np.random.seed(self.rand_seed_wtc + random_seed_atc)
        instr_distr = np.random.choice(instr_numbs, p=instr_probs, size=self.N, replace=True)

        # create a list of instruction names to be recorded in the output file
        # and randomly choose a value and duration of instruction if applicable
        instr_distr_list = instr_distr.tolist()

        instr_names_list = []
        instr_value_list = []
        instr_dur_list = []

        for val in instr_distr_list:
            rand_seed_counter = 0
            instr_names_list.append(instr_names[val])

            if instr_names[val] in ('direct', '0', 'vs'):
                """ if direct, no instruction or vertical speed is given
                append zeros to instruction value and duration 
                since these types do not use value and duration"""
                instr_value_list.append(0)
                instr_dur_list.append(0)

            elif instr_names[val] == 'heading':
                """ if it is a heading instruction, 
                sample value (left or right turn) and duration """
                # np.random.seed(rand_seed + rand_seed_counter)
                hdg_val = np.random.choice([-20, -15, -10, 10, 15, 20], size=1)

                # np.random.seed(rand_seed + rand_seed_counter)
                hdg_dur = np.random.choice(range(30, 60*4+1), size=1)

                rand_seed_counter += 1
                instr_value_list.append(float(hdg_val))
                instr_dur_list.append(float(hdg_dur))

            elif instr_names[val] == 'speed':
                """ if it is a speed instruction, 
                sample mach value from a distribution and duration """
                # np.random.seed(rand_seed + rand_seed_counter)
                spd_val = np.random.choice(ATC['mach'][0], p=ATC['mach'][1], size=1)

                # np.random.seed(rand_seed + rand_seed_counter)
                spd_dur = np.random.choice(range(30, 60*3+1), size=1)

                rand_seed_counter += 1
                instr_value_list.append(float(spd_val))
                instr_dur_list.append(float(spd_dur))

        self.atc['ATC_instr'] = instr_names_list
        self.atc['ATC_value'] = instr_value_list
        self.atc['ATC_dur'] = instr_dur_list

        # write values into a file
        self.write_to_file(self.atc, '')

    def check_sample(self, sample, distr, params, param_type):

        # determine minimum and maximum parameter values
        # based on a parameter type
        if param_type == 'bank':
            min_val = BANK['min']
            max_val = BANK['max']

        elif param_type == "leveloff dur":
            min_val = TEMP_LEVELOFF['min dur']
            max_val = TEMP_LEVELOFF['max dur']

        elif param_type == "leveloff fl":
            min_val = TEMP_LEVELOFF['min fl']
            max_val = TEMP_LEVELOFF['max fl']

        elif param_type == 'temp':
            min_val = AIR_TEMP['min temp']
            max_val = AIR_TEMP['max temp']

        elif param_type == "lapse":
            min_val = AIR_TEMP['min lapse']
            max_val = AIR_TEMP['max lapse']

        elif param_type == "cas":
            min_val = SPEED_CAS['min'][self.wtc]
            max_val = SPEED_CAS['max'][self.wtc]

        elif param_type == "mach":
            min_val = SPEED_MACH['min'][self.wtc]
            max_val = SPEED_MACH['max'][self.wtc]

        elif param_type == "roc":
            min_val = ROCD['min roc']
            max_val = ROCD['max roc']

        elif param_type == "rod":
            min_val = ROCD['min rod']
            max_val = ROCD['max rod']

        # check for values that are smaller than minimum
        # and larger than maximum and replace them by sampling
        # from a distribution until value is within the limits
        small_values = np.where(sample < min_val)
        large_values = np.where(sample > max_val)

        wrong_values = np.concatenate((small_values[0], large_values[0]))

        if param_type in ['bank', 'leveloff dur', 'temp', 'lapse', 'cas', 'mach', 'roc', 'rod']:
            """ these are continuous variables"""
            # this counter for seed (increase seed) is used to ensure
            # that all values that need to be replaced are different
            # otherwise, if the global seed is used,
            # all replaced values are the same
            # thus creating frequency of values where it should not be
            rand_seed_counter = 0
            for val_index in wrong_values:
                # initialize with a large (dummy) negative number
                # to ensure that this algorithm works
                # for all types of parameters
                new_val = -1000000

                while new_val < min_val or new_val > max_val:
                    new_val = distr.rvs(*params, size=1,
                                        random_state=self.rand_seed_wtc + random_seeds_for_resample[param_type] + rand_seed_counter)
                    rand_seed_counter += 1
                sample[val_index] = new_val

            small_values = np.where(sample < min_val)
            large_values = np.where(sample > max_val)

            return sample.tolist()

        elif param_type in ['leveloff fl']:
            """ these are discrete variables"""
            rand_seed_counter = 0
            for val_index in wrong_values:
                # initialize with a large (dummy) negative number
                # to ensure that this algorithm works
                # for all types of parameters
                new_val = -1000000
                while new_val < min_val or new_val > max_val:
                    np.random.seed(self.rand_seed_wtc + random_seeds_for_resample[param_type] + rand_seed_counter)
                    new_val = np.random.choice(distr, p=params)
                    rand_seed_counter += 1
                sample[val_index] = new_val

            small_values = np.where(sample < min_val)
            large_values = np.where(sample > max_val)

            return sample.tolist()

    def sample_fpl(self):
        np.random.seed(self.rand_seed_wtc + random_seeds_for_resample['fpl'])
        fpl_ids = np.random.choice(range(0, 10), p=[0.1]*10, size=self.N, replace=True)
        self.fpl["FPL"] = fpl_ids.tolist()
        self.write_to_file(self.fpl, "N")

    def write_to_file(self, sample, param_type):

        for key, values in sample.items():

            self.sheet_counter = 0
            self.sheet = self.workbook.get_sheet("MCsim_%s" % self.sheet_counter)

            for r, val in enumerate(values):

                if self.rows_count == 0:
                    self.row = self.sheet.row(self.rows_count)
                    if param_type:
                        self.row.write(self.columns_count, "%s_%s" % (param_type, key))
                    else:
                        self.row.write(self.columns_count, "%s" % key)

                self.rows_count += 1
                self.row = self.sheet.row(self.rows_count)

                if isinstance(val, int) or isinstance(val, float):
                    if param_type == "M":  # mach number
                        self.row.write(self.columns_count, float('%.3f' % float(val)))

                    elif param_type == "ROC" or param_type == "ROD" or param_type == "LVLOFF_FL" or param_type == 'N':
                        self.row.write(self.columns_count, int(val))

                    elif param_type == 'U' or param_type == 'V':
                        self.row.write(self.columns_count, float('%.6f' % float(val)))

                    else:
                        self.row.write(self.columns_count, float('%.1f' % float(val)))

                elif isinstance(val, str):
                    self.row.write(self.columns_count, val)

                if self.rows_count == int(self.N / 10):
                    self.sheet_counter += 1
                    try:
                        self.sheet = self.workbook.add_sheet("MCsim_%s" % self.sheet_counter)
                    except:
                        self.sheet = self.workbook.get_sheet("MCsim_%s" % self.sheet_counter)
                    self.rows_count = 0
                    # self.columns_count = 0

            self.columns_count += 1


if __name__ == "__main__":
    path = (os.path.join(os.path.dirname(__file__), "inputs/"))
    i = InputPreparation(10000, 'L', path)
