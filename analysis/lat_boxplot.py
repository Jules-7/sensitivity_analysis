""" 4 SIM SETTINGS: ERRORS DISTRIBUTION WITH LOOK-AHEAD TIME PER FLIGHT PHASE

This module reads the errors computed for 5, 10, 15 and 20 minutes
look-ahead time for 1 to 4 simulation settings
and plots boxplots per phase of a real flight

This module is the MAIN VERSION OF BOX PLOTS FOR LOOK-AHEAD TIMES"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
title_size = 20
label_size = 18
legend_size = 14

wtc = 'L'
OUTLIERS = False  # plot outliers
PLOT = True  # save plots

path_to_errors_sim = '/media/julia/UserData/MC/lat_errors_%s/' % wtc

lats = [5, 10, 15, 20]
# lats = [20]

# simulation settings
sims = [1, 2, 3, 4]
# sims = [2]
# sims = [1, 2]
# max_cte = 20
# max_ate = 20

for phase_name in ['CLIMB', 'CRUISE', 'DESCENT']:

    fig, axes = plt.subplots(3, 1, figsize=(6, 4), sharex=True)

    print('\n', phase_name)

    for lat in lats:

        ATE_sim1, CTE_sim1, ALTE_sim1 = [], [], []

        print('%s lat'%lat)

        for sim in sims:
            # for sim in [1, 2]:

            print('sim', sim)

            counter = 0
            ate_list, cte_list, alte_list = [], [], []

            for route_id in range(0, 10):

                try:
                    for flight in os.listdir(path_to_errors_sim + 'sim_res_%s/%s/%s/' % (sim, route_id, lat)):

                        try:
                            flight_data = pd.read_csv(path_to_errors_sim + 'sim_res_%s/%s/%s/' % (sim, route_id, lat) + flight,
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
                            # it was computed initially wrong
                            alte_phase_errors = alte_phase_errors / 0.3048 * 3.28084

                            # if phase_name == 'CLIMB' and np.any(np.where(np.abs(ate_phase_errors)>max_ate)):
                            #     counter += 1
                            #     print(flight)

                            ate_list += ate_phase_errors.tolist()
                            cte_list += cte_phase_errors.tolist()
                            alte_list += alte_phase_errors.tolist()

                        except pd.errors.EmptyDataError:
                            pass
                except FileNotFoundError:
                    pass

            # print("\n%s values are larger than %s nm\n\n"%(counter, max_cte))

            ATE_sim1.append(ate_list)
            CTE_sim1.append(cte_list)
            ALTE_sim1.append(alte_list)

            # ate_mean.append(np.mean(np.array(ate_list)))
            # cte_mean.append(np.mean(np.array(cte_list)))
            # alte_mean.append(np.mean(np.array(alte_list)))

            # ate_median.append(np.median(np.array(ate_list)))
            # cte_median.append(np.median(np.array(cte_list)))
            # alte_median.append(np.median(np.array(alte_list)))

        # print('\nate mean\n', ate_mean)
        # print('\nate median\n', ate_median)

        # print('\ncte mean\n', cte_mean)
        # print('\ncte median\n', cte_median)

        # print('\nalte mean\n', alte_mean)
        # print('\nalte median\n', alte_median)

        # print("plotting")

        pos = [lat - 1.5, lat - 0.5, lat + 0.5, lat + 1.5][:len(sims)]
        width = [0.8, 0.8, 0.8, 0.8][:len(sims)]  # tuple([4] * len(lats))

        if not OUTLIERS:
            # no outliers
            bp1 = axes[0].boxplot(ATE_sim1, positions=pos, widths=width, sym='')
            bp2 = axes[1].boxplot(CTE_sim1, positions=pos, widths=width, sym='')
            bp3 = axes[2].boxplot(ALTE_sim1, positions=pos, widths=width, sym='')
        else:
            # with outliers
            bp1 = axes[0].boxplot(ATE_sim1, positions=pos, widths=width)
            bp2 = axes[1].boxplot(CTE_sim1, positions=pos, widths=width)
            bp3 = axes[2].boxplot(ALTE_sim1, positions=pos, widths=width)

        if lat == 20:
            # print the values of the whiskers
            print('%s %s 20 lat' % (phase_name, wtc))
            print('ATE whiskers values', [item.get_ydata()[1] for item in bp1['whiskers']])
            print('CTE whiskers values', [item.get_ydata()[1] for item in bp2['whiskers']])
            print('AE whiskers values', [item.get_ydata()[1] for item in bp3['whiskers']])

        # change color and linewidth of the medians
        for median in bp1['medians']:
            median.set(color='r', linewidth=1)
        for median in bp2['medians']:
            median.set(color='r', linewidth=1)
        for median in bp3['medians']:
            median.set(color='r', linewidth=1)

        # ## change the style of fliers and their fill
        # for flier in bp1['fliers']:
        #     flier.set(marker='+', color='b', alpha=0.1)
        # for flier in bp2['fliers']:
        #     flier.set(marker='+', color='b', alpha=0.1)
        # for flier in bp3['fliers']:
        #     flier.set(marker='+', color='b', alpha=0.1)

        axes[0].set_ylabel('ATE, \nnmi', fontsize=label_size)
        axes[1].set_ylabel('CTE, \nnmi', fontsize=label_size)
        axes[2].set_ylabel('AE, ft', fontsize=label_size)

        # axes[0].set_xlabel('Look-ahead time, minutes')
        # axes[1].set_xlabel('Look-ahead time, minutes')
        axes[2].set_xlabel('Look-ahead time, minutes', fontsize=label_size)

        axes[0].set_xlim(2, lats[-1]+3)
        axes[1].set_xlim(2, lats[-1]+3)
        axes[2].set_xlim(2, lats[-1]+3)

        axes[2].set_xticks([5, 10, 15, 20])
        axes[2].set_xticklabels(['5', '10', '15', '20'])

        # axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45)

        # axes[0].grid(True)
        # axes[1].grid(True)
        # axes[2].grid(True)

        axes[0].yaxis.grid(True)
        axes[1].yaxis.grid(True)
        axes[2].yaxis.grid(True)

    # plt.suptitle("%s %s"% (phase_name, wtc))
    if PLOT:
        plt.tight_layout()
        if OUTLIERS:
            plt.savefig('%s_%s_outliers.png' % (phase_name, wtc))
        else:
            plt.savefig('%s_%s.png'%(phase_name, wtc))
    # plt.show()

