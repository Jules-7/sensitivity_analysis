""" 1 SIM SETTING: ERRORS DISTRIBUTION WITH LOOK-AHEAD TIME PER FLIGHT PHASE

This module reads the errors computed for 5, 10, 15 and 20 minutes
look-ahead time for 1 simulation setting
and plots histogram per phase of a real flight

This module is the MAIN VERSION OF BOX PLOTS FOR LOOK-AHEAD TIMES"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
title_size = 18
label_size = 16
legend_size = 14

wtc = 'M'

path_to_errors = '/media/julia/UserData/MC/lat_errors_%s/' % wtc


lats = [5, 10, 15, 20]
max_cte = 20
max_ate = 20

for phase_name in ['CLIMB', 'CRUISE', 'DESCENT']:

    print('\n', phase_name)

    ate_mean, cte_mean, alte_mean = [], [], []
    ate_median, cte_median, alte_median = [], [], []

    for lat in lats:

        ATE, CTE, ALTE = [], [], []

        print('%s lat' % lat)

        counter = 0
        ate_list, cte_list, alte_list = [], [], []

        for route_id in range(0, 10):

            try:
                for flight in os.listdir(path_to_errors + '%s/%s/' % (route_id, lat)):

                    try:
                        flight_data = pd.read_csv(path_to_errors + '%s/%s/' % (route_id, lat) + flight,
                                                  header=None, usecols=[0, 1, 2, 3, 4, 5, 6])

                        ate = flight_data[0].values
                        cte = flight_data[1].values
                        alte = flight_data[2].values
                        u = flight_data[3].values
                        v = flight_data[4].values
                        t_phase = flight_data[5].values
                        p_phase = flight_data[6].values

                        phase_ids = np.where(t_phase == phase_name)
                        # cruise_phase = np.where(t_phase == 'CRUISE')
                        # descent_phase = np.where(t_phase == 'DESCENT')

                        ate_phase_errors = ate[phase_ids]
                        cte_phase_errors = cte[phase_ids]
                        alte_phase_errors = alte[phase_ids]

                        # correct altitude error
                        # it was wrongly computed initially
                        alte_phase_errors = alte_phase_errors / 0.3048 * 3.28084

                        # if phase_name == 'CLIMB' and np.any(np.where(np.abs(ate_phase_errors)>max_ate)):
                        #     counter += 1
                        #     print(flight)

                        ATE += ate_phase_errors.tolist()
                        CTE += cte_phase_errors.tolist()
                        ALTE += alte_phase_errors.tolist()

                    except pd.errors.EmptyDataError:
                        pass
            except FileNotFoundError:
                pass

            # print("\n%s values are larger than %s nm\n\n"%(counter, max_cte))

        print('plotting')

        fig, axes = plt.subplots(3, 1, figsize=(6, 4))

        axes[0].hist(ATE)
        axes[1].hist(CTE)
        axes[2].hist(ALTE)

        axes[0].set_xlabel('ATE, \nnmi', fontsize=label_size)
        axes[1].set_xlabel('CTE, \nnmi', fontsize=label_size)
        axes[2].set_xlabel('AE, ft', fontsize=label_size)

        # axes[0].set_xlabel('Look-ahead time, minutes')
        # axes[1].set_xlabel('Look-ahead time, minutes')
        # axes[2].set_xlabel('Look-ahead time, minutes', fontsize=label_size)

        # axes[0].set_xlim(0, lats[-1]+5)
        # axes[1].set_xlim(0, lats[-1]+5)
        # axes[2].set_xlim(0, lats[-1]+5)

        # axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45)

        # axes[0].grid(True)
        # axes[1].grid(True)
        # axes[2].grid(True)

        plt.suptitle('%s %s %s lat' % (phase_name, wtc, lat))

plt.show()
