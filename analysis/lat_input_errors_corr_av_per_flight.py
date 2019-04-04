""" INPUT-ERRORS CORRELATION

This module allows to:
- visualize the input-error correlation (scatter or 2-d histogram)
- output correlation in a tabular format as Pearson or Spearman coefficient

This is the MAIN VERSION of code used to obtain the correlation between
the inputs and the prediction errors"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import xlwt
# import statsmodels.formula.api as smf
import seaborn as sns
import matplotlib as mpl

label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
title_size = 16
legend_size = 14

# set correlation graphical output
PLOT = True  # plot input-error correlation as scatter or 2-d histogram
SCATTER = False  # scatter or 2-d histogram

# set correlation tabular output
PROB = False  # save p-values of coefficients in a file, otherwise save correlation coefficients
SPEARMAN = False  # Spearman or Pearson coefficient


lats = [5, 10, 15, 20]
# run each WTC and simulation setup one-by-one
wtc = 'L'
sim = '4'

path_to_input_data_file = os.path.dirname(os.path.realpath(__file__)) + '/MCsim_%s.xls' % wtc
path_to_errors = '/media/julia/UserData/MC/lat_errors_%s/sim_res_%s/' % (wtc, sim)
path_to_save_fig = '/media/julia/UserData/MCsim_results/corr_factors_2d_hist/%s_scatter/' % wtc

workbook = xlwt.Workbook()

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


for lat in lats:

    print('\n\nlook-ahead time ', lat)

    for phase_name in ['CLIMB', 'CRUISE', 'DESCENT']:

        print(phase_name)

        errors = {'ATE': [],
                  'CTE': [],
                  'AE': []}

        for error_name, errors_list in errors.items():

            print(error_name, 'error')

            sheet = workbook.add_sheet('%s_%s_%s' % (lat, phase_name, error_name))

            # errors are computed for each route separately
            routes_errors = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: []}

            for route_id in range(0, 10):
                # to collect data per route, initialize the lists here
                input_temp, input_lapse, input_ATC, input_bank, input_vs = [], [], [], [], []
                input_cas, input_mach, input_lvl_fl, input_lvl_dur, input_u, input_v = [], [], [], [], [], []
                errors_array = []

                print('\nroute id', route_id)

                for flight in os.listdir(path_to_errors + '%s/%s/'%(route_id, lat)):

                    try:
                        flight_data = pd.read_csv(path_to_errors + '%s/%s/'%(route_id, lat) + flight,
                                                  header=None, usecols=[0, 1, 2, 3, 4, 5, 6])
                        errors_values = None
                        if error_name == 'ATE':
                            errors_values = flight_data[0].values
                        elif error_name == 'CTE':
                            errors_values = flight_data[1].values
                        elif error_name == 'AE':
                            errors_values = flight_data[2].values
                            errors_values = errors_values / 0.3048 * 3.28084  # correct altitude

                        u = flight_data[3].values
                        v = flight_data[4].values

                        t_phase = flight_data[5].values
                        p_phase = flight_data[6].values

                        phase_ids = np.where(t_phase == phase_name)

                        # !!! compute average error over all values !!!
                        # average of all values (positive and negative)
                        av_errors = np.average(errors_values[phase_ids])

                        # if average is computed over absolute values only,
                        # this results in very weak correlation
                        # average of absolute values
                        # av_errors = np.average(np.abs(errors_values[phase_ids]))

                        nan_values = np.isnan(av_errors)

                        if not np.any(nan_values):

                            errors_array.append(av_errors)

                            u_av = np.average(u[phase_ids])
                            v_av = np.average(v[phase_ids])

                            designator = flight[:-4].split('_')
                            sheet_num = int(designator[1])
                            flight_id = int(designator[3])

                            # extract inputs for this flight
                            row_inputs = sheets[sheet_num].loc[flight_id, :]

                            input_temp.append(row_inputs['TEMP_OFFSET'])
                            input_lapse.append(row_inputs['LAPSE_RATE'])

                            ATC = row_inputs['ATC_instr']
                            if ATC == '0':
                                input_ATC.append(0)
                            elif ATC == 'direct':
                                input_ATC.append(1)
                            elif ATC == 'heading':
                                input_ATC.append(2)
                            elif ATC == 'speed':
                                input_ATC.append(3)
                            elif ATC == 'vs':
                                input_ATC.append(4)
                            input_u.append(u_av)
                            input_v.append(v_av)

                            if phase_name == 'CLIMB':
                                input_bank.append(row_inputs['BANK_Climb'])
                                input_lvl_fl.append(row_inputs['LVLOFF_FL_Climb'])
                                input_lvl_dur.append(row_inputs['LVLOFF_DUR_Climb'])
                                input_cas.append(row_inputs['CAS_Climb'])
                                input_mach.append(row_inputs['M_Climb'])
                                input_vs.append(row_inputs['ROC'])

                            elif phase_name == 'CRUISE':

                                input_bank.append(row_inputs['BANK_Cruise'])
                                input_mach.append(row_inputs['M_Cruise'])

                            elif phase_name == 'DESCENT':
                                input_bank.append(row_inputs['BANK_Descent'])
                                input_lvl_fl.append(row_inputs['LVLOFF_FL_Descent'])
                                input_lvl_dur.append(row_inputs['LVLOFF_DUR_Descent'])
                                input_cas.append(row_inputs['CAS_Descent'])
                                input_mach.append(row_inputs['M_Descent'])
                                input_vs.append(row_inputs['ROD'])

                    except pd.errors.EmptyDataError:
                        pass

                if PLOT:
                    print('plotting')

                    fig, axes = plt.subplots(3, 4, figsize=(12, 9))
                    fig.canvas.set_window_title('%s %s %s %s min LAT route %s' % (wtc, phase_name, error_name, lat, route_id))
                    # plt.suptitle('%s %s %s min LAT route %s' % (wtc, phase_name, lat, route_id))

                    if error_name == 'AE':
                        if np.max(errors_array) > 500:
                            bins_error = np.arange(np.min(errors_array), np.max(errors_array) + 100, 100)
                        else:
                            bins_error = np.arange(np.min(errors_array), np.max(errors_array) + 10, 10)

                    else:
                        if np.max(errors_array) > 50:
                            bins_error = np.arange(np.min(errors_array), np.max(errors_array) + 5, 5)
                        elif np.max(errors_array) < 6:
                            bins_error = np.arange(np.min(errors_array), np.max(errors_array) + 0.1, 0.1)
                        else:
                            bins_error = np.arange(np.min(errors_array), np.max(errors_array) + 1, 1)

                    if SCATTER:  # scatter plots
                        axes[0, 0].scatter(input_temp, errors_array)
                        axes[0, 1].scatter(input_lapse, errors_array)
                        axes[0, 2].scatter(input_u, errors_array)
                        axes[0, 3].scatter(input_v, errors_array)

                        axes[1, 0].scatter(input_bank, errors_array)
                        axes[1, 1].scatter(input_ATC, errors_array)
                        axes[1, 2].scatter(input_mach, errors_array)

                        if phase_name != 'CRUISE':
                            axes[2, 0].scatter(input_cas, errors_array)
                            axes[2, 1].scatter(input_vs, errors_array)
                            axes[2, 2].scatter(input_lvl_fl, errors_array)
                            axes[2, 3].scatter(input_lvl_dur, errors_array)

                    else:  # 2-d histogram

                        axes[0, 0].hist2d(input_temp, errors_array,
                                          bins=(np.arange(np.min(input_temp), np.max(input_temp)+1, 1),
                                                bins_error), normed=True, cmap=plt.cm.Greys)
                        axes[0, 1].hist2d(input_lapse, errors_array,
                                          bins=(np.arange(np.min(input_lapse), np.max(input_lapse)+0.2, 0.2),
                                                bins_error), normed=True, cmap=plt.cm.Greys)
                        axes[0, 2].hist2d(input_u, errors_array,
                                          bins=(np.arange(np.min(input_u), np.max(input_u) + 5, 5),
                                                bins_error), normed=True, cmap=plt.cm.Greys)
                        axes[0, 3].hist2d(input_v, errors_array,
                                          bins=(np.arange(np.min(input_v), np.max(input_v) + 5, 5),
                                                bins_error), normed=True, cmap=plt.cm.Greys)


                        axes[1, 0].hist2d(input_bank, errors_array,
                                          bins=(np.arange(5, 35 + 2.5, 2.5),
                                                bins_error), normed=True, cmap=plt.cm.Greys)
                        axes[1, 1].hist2d(input_ATC, errors_array,
                                          bins=(np.arange(0, 4 + 2, 1),
                                                bins_error), normed=True, cmap=plt.cm.Greys)
                        axes[1, 2].hist2d(input_mach, errors_array,
                                          bins=(np.arange(np.min(input_mach), np.max(input_mach) + 0.005, 0.005),
                                                bins_error), normed=True, cmap=plt.cm.Greys)

                        if phase_name != 'CRUISE':
                            axes[2, 0].hist2d(input_cas, errors_array,
                                              bins=(np.arange(np.min(input_cas), np.max(input_cas) + 5, 5),
                                                    bins_error), normed=True, cmap=plt.cm.Greys)
                            axes[2, 1].hist2d(input_vs, errors_array,
                                              bins=(np.arange(np.min(input_vs), np.max(input_vs) + 100, 100),
                                                    bins_error), normed=True, cmap=plt.cm.Greys)
                            axes[2, 2].hist2d(input_lvl_fl, errors_array,
                                              bins=(np.arange(0, np.max(input_lvl_fl) + 20, 20),
                                                    bins_error), normed=True, cmap=plt.cm.Greys)
                            axes[2, 3].hist2d(input_lvl_dur, errors_array,
                                              bins=(np.arange(0, np.max(input_lvl_dur) + 30, 30),
                                                    bins_error), normed=True, cmap=plt.cm.Greys)

                    # perform fitting to a line
                    fit_temp = np.polyfit(input_temp, errors_array, 1)
                    fit_fn_temp = np.poly1d(fit_temp)

                    fit_lapse = np.polyfit(input_lapse, errors_array, 1)
                    fit_fn_lapse = np.poly1d(fit_lapse)

                    fit_u = np.polyfit(input_u, errors_array, 1)
                    fit_fn_u = np.poly1d(fit_u)

                    fit_v = np.polyfit(input_v, errors_array, 1)
                    fit_fn_v = np.poly1d(fit_v)

                    fit_bank = np.polyfit(input_bank, errors_array, 1)
                    fit_fn_bank = np.poly1d(fit_bank)

                    fit_atc = np.polyfit(input_ATC, errors_array, 1)
                    fit_fn_atc = np.poly1d(fit_atc)

                    fit_mach = np.polyfit(input_mach, errors_array, 1)
                    fit_fn_mach = np.poly1d(fit_mach)

                    if phase_name != 'CRUISE':
                        fit_cas = np.polyfit(input_cas, errors_array, 1)
                        fit_fn_cas = np.poly1d(fit_cas)

                        fit_vs = np.polyfit(input_vs, errors_array, 1)
                        fit_fn_vs = np.poly1d(fit_vs)

                        fit_lvl_fl = np.polyfit(input_lvl_fl, errors_array, 1)
                        fit_fn_lvl_fl = np.poly1d(fit_lvl_fl)

                        fit_lvl_dur = np.polyfit(input_lvl_dur, errors_array, 1)
                        fit_fn_lvl_dur = np.poly1d(fit_lvl_dur)

                    # plot linear regression line
                    axes[0, 0].plot(input_temp, fit_fn_temp(input_temp), '*r', markersize=4, alpha=0.1)
                    axes[0, 1].plot(input_lapse, fit_fn_lapse(input_lapse), '*r', markersize=4, alpha=0.1)
                    axes[0, 2].plot(input_u, fit_fn_u(input_u), '*r', markersize=4, alpha=0.1)
                    axes[0, 3].plot(input_v, fit_fn_v(input_v), '*r', markersize=4, alpha=0.1)

                    axes[1, 0].plot(input_bank, fit_fn_bank(input_bank), '*r', markersize=4, alpha=0.1)
                    axes[1, 1].plot(input_ATC, fit_fn_atc(input_ATC), '*r', markersize=4, alpha=0.1)
                    axes[1, 2].plot(input_mach, fit_fn_mach(input_mach), '*r', markersize=4, alpha=0.1)

                    if phase_name != 'CRUISE':
                        axes[2, 0].plot(input_cas, fit_fn_cas(input_cas), '*r', markersize=4, alpha=0.1)
                        axes[2, 1].plot(input_vs, fit_fn_vs(input_vs), '*r', markersize=4, alpha=0.1)
                        axes[2, 2].plot(input_lvl_fl, fit_fn_lvl_fl(input_lvl_fl), '*r', markersize=4, alpha=0.1)
                        axes[2, 3].plot(input_lvl_dur, fit_fn_lvl_dur(input_lvl_dur), '*r', markersize=4, alpha=0.1)

                    # compute correlation coefficient
                    corr_coeff_temp, p_val_temp = st.pearsonr(input_temp, errors_array)
                    corr_coeff_lapse, p_val_lapse = st.pearsonr(input_lapse, errors_array)
                    corr_coeff_u, p_val_u = st.pearsonr(input_u, errors_array)
                    corr_coeff_v, p_val_v = st.pearsonr(input_v, errors_array)

                    corr_coeff_bank, p_val_bank = st.pearsonr(input_bank, errors_array)
                    corr_coeff_atc, p_val_atc = st.pearsonr(input_ATC, errors_array)
                    corr_coeff_mach, p_val_mach = st.pearsonr(input_mach, errors_array)

                    try:
                        corr_coeff_cas, p_val_cas = st.pearsonr(input_cas, errors_array)
                        corr_coeff_vs, p_val_vs = st.pearsonr(input_vs, errors_array)
                        corr_coeff_lvl_fl, p_val_lvl_fl = st.pearsonr(input_lvl_fl, errors_array)
                        corr_coeff_lvl_dur, p_val_lvl_dur = st.pearsonr(input_lvl_dur, errors_array)
                    except:
                        pass

                    # plot correlation coefficients as titles of each subplot
                    axes[0, 0].set_title('corr %.2f' % corr_coeff_temp, fontsize=title_size)
                    axes[0, 1].set_title('corr %.2f' % corr_coeff_lapse, fontsize=title_size)
                    axes[0, 2].set_title('corr %.2f' % corr_coeff_u, fontsize=title_size)
                    axes[0, 3].set_title('corr %.2f' % corr_coeff_v, fontsize=title_size)

                    axes[1, 0].set_title('corr %.2f' % corr_coeff_bank, fontsize=title_size)
                    axes[1, 1].set_title('corr %.2f' % corr_coeff_atc, fontsize=title_size)
                    axes[1, 2].set_title('corr %.2f' % corr_coeff_mach, fontsize=title_size)

                    try:
                        axes[2, 0].set_title('corr %.2f' % corr_coeff_cas, fontsize=title_size)
                        axes[2, 1].set_title('corr %.2f' % corr_coeff_vs, fontsize=title_size)
                        axes[2, 2].set_title('corr %.2f' % corr_coeff_lvl_fl, fontsize=title_size)
                        axes[2, 3].set_title('corr %.2f' % corr_coeff_lvl_dur, fontsize=title_size)
                    except:
                        pass

                    axes[0, 0].set_xlabel('Temperature offset, K', fontsize=label_size)
                    axes[0, 1].set_xlabel('Lapse rate, K/km', fontsize=label_size)
                    axes[0, 2].set_xlabel('U wind, kts', fontsize=label_size)
                    axes[0, 3].set_xlabel('V wind, kts', fontsize=label_size)

                    axes[1, 0].set_xlabel('Bank angle, deg', fontsize=label_size)
                    axes[1, 1].set_xlabel('ATC, -', fontsize=label_size)
                    axes[1, 1].set_xticks([0, 1, 2, 3, 4])
                    axes[1, 1].set_xticklabels(['0', 'D', 'H', 'S', 'VS'])

                    axes[1, 2].set_xlabel('Mach, -', fontsize=label_size)
                    # axes[1, 2].set_xticklabels(['%.2f'%float(val) for val in axes[1, 2].get_xticklabels()])

                    try:
                        axes[2, 0].set_xlabel('CAS, kts', fontsize=label_size)
                        axes[2, 1].set_xlabel('VS, fpm', fontsize=label_size)
                        axes[2, 2].set_xlabel('Temp level-off FL, -', fontsize=label_size)
                        axes[2, 3].set_xlabel('Temp level-off duration, sec', fontsize=label_size)
                    except:
                        pass

                    for row in range(3):
                        axes[row, 0].set_ylabel('%s, ft' % error_name if error_name == 'AE' else '%s, nmi' % error_name, fontsize=label_size)
                    plt.tight_layout()
                    fig.savefig(path_to_save_fig + '%s %s %s %s min LAT route %s.png' % (wtc, phase_name, error_name, lat, route_id))

                else:  # no plotting, save correlation in a tabular format

                    # compute correlation coefficient
                    if SPEARMAN:
                        corr_coeff_temp, p_val_temp = st.spearmanr(input_temp, errors_array)
                        corr_coeff_lapse, p_val_lapse = st.spearmanr(input_lapse, errors_array)
                        corr_coeff_u, p_val_u = st.spearmanr(input_u, errors_array)
                        corr_coeff_v, p_val_v = st.spearmanr(input_v, errors_array)

                        corr_coeff_bank, p_val_bank = st.spearmanr(input_bank, errors_array)
                        corr_coeff_atc, p_val_atc = st.spearmanr(input_ATC, errors_array)
                        corr_coeff_mach, p_val_mach = st.spearmanr(input_mach, errors_array)

                        if phase_name != 'CRUISE':
                            corr_coeff_cas, p_val_cas = st.spearmanr(input_cas, errors_array)
                            corr_coeff_vs, p_val_vs = st.spearmanr(input_vs, errors_array)
                            corr_coeff_lvl_fl, p_val_lvl_fl = st.spearmanr(input_lvl_fl, errors_array)
                            corr_coeff_lvl_dur, p_val_lvl_dur = st.spearmanr(input_lvl_dur, errors_array)

                    else:  # Pearson
                        corr_coeff_temp, p_val_temp = st.pearsonr(input_temp, errors_array)
                        corr_coeff_lapse, p_val_lapse = st.pearsonr(input_lapse, errors_array)
                        corr_coeff_u, p_val_u = st.pearsonr(input_u, errors_array)
                        corr_coeff_v, p_val_v = st.pearsonr(input_v, errors_array)

                        corr_coeff_bank, p_val_bank = st.pearsonr(input_bank, errors_array)
                        corr_coeff_atc, p_val_atc = st.pearsonr(input_ATC, errors_array)
                        corr_coeff_mach, p_val_mach = st.pearsonr(input_mach, errors_array)

                        if phase_name != 'CRUISE':
                            corr_coeff_cas, p_val_cas = st.pearsonr(input_cas, errors_array)
                            corr_coeff_vs, p_val_vs = st.pearsonr(input_vs, errors_array)
                            corr_coeff_lvl_fl, p_val_lvl_fl = st.pearsonr(input_lvl_fl, errors_array)
                            corr_coeff_lvl_dur, p_val_lvl_dur = st.pearsonr(input_lvl_dur, errors_array)

                    if PROB:  # save p-value
                        routes_errors[route_id].append(p_val_temp)
                        routes_errors[route_id].append(p_val_lapse)
                        routes_errors[route_id].append(p_val_u)
                        routes_errors[route_id].append(p_val_v)
                        routes_errors[route_id].append(p_val_bank)
                        routes_errors[route_id].append(p_val_atc)
                        routes_errors[route_id].append(p_val_mach)

                        if phase_name != 'CRUISE':
                            routes_errors[route_id].append(p_val_cas)
                            routes_errors[route_id].append(p_val_vs)
                            routes_errors[route_id].append(p_val_lvl_fl)
                            routes_errors[route_id].append(p_val_lvl_dur)

                    else:  # save correlation coefficient
                        routes_errors[route_id].append(corr_coeff_temp)
                        routes_errors[route_id].append(corr_coeff_lapse)
                        routes_errors[route_id].append(corr_coeff_u)
                        routes_errors[route_id].append(corr_coeff_v)
                        routes_errors[route_id].append(corr_coeff_bank)
                        routes_errors[route_id].append(corr_coeff_atc)
                        routes_errors[route_id].append(corr_coeff_mach)

                        if phase_name != 'CRUISE':
                            routes_errors[route_id].append(corr_coeff_cas)
                            routes_errors[route_id].append(corr_coeff_vs)
                            routes_errors[route_id].append(corr_coeff_lvl_fl)
                            routes_errors[route_id].append(corr_coeff_lvl_dur)

                    # print('Correlation coefficients')
                    # print('Temperature %.2f \t'%corr_coeff_temp)
                    # print('Lapse rate %.2f \t'%corr_coeff_lapse)
                    # print('Wind U %.2f \t'%corr_coeff_u)
                    # print('Wind V %.2f \t'%corr_coeff_v)
                    # print('Bank %.2f \t'%corr_coeff_bank)
                    # print('ATC %.2f \t'%corr_coeff_atc)
                    # print('Mach %.2f \t'%corr_coeff_mach)
                    #
                    # if phase_name != 'CRUISE':
                    #     print('CAS %.2f \t'%corr_coeff_cas)
                    #     print('VS %.2f \t'%corr_coeff_vs)
                    #     print('LVL FL %.2f \t'%corr_coeff_lvl_fl)
                    #     print('LVL DUR %.2f \t'%corr_coeff_lvl_dur)

            columns_count = 0
            rows_count = 0
            for row_name in [error_name, 'Temperature', 'Lapse rate', 'U wind', 'V wind', 'Bank', 'ATC', 'Mach',
                             'CAS', 'VS', 'LVL FL', 'LVL DUR']:
                row = sheet.row(rows_count)
                row.write(columns_count, row_name)
                rows_count += 1

            # compute average of absolute values for each input
            # and record as the last column -> column 10
            # (i.e. over all routes)

            for j in range(len(routes_errors[0])):
                sum_abs_values = 0
                for i in range(10):
                    sum_abs_values += abs(routes_errors[i][j])

                av_abs_error = sum_abs_values / 10.
                routes_errors[10].append(av_abs_error)

            columns_count += 1

            for key, coeff_values in routes_errors.items():
                rows_count = 0

                row = sheet.row(rows_count)

                row.write(columns_count, int(key))

                for val in coeff_values:
                    rows_count += 1
                    row = sheet.row(rows_count)
                    row.write(columns_count, float('%.2f' % val))

                columns_count += 1

if SPEARMAN:
    xls_file = os.path.dirname(os.path.realpath(__file__)) + '/%s_spearman_corr_coeff_sim_%s.xls' % (wtc, sim)
else:
    xls_file = os.path.dirname(os.path.realpath(__file__)) + '/%s_pearson_corr_coeff_sim_%s.xls' % (wtc, sim)
workbook.save(xls_file)

if PLOT:
    plt.show()
