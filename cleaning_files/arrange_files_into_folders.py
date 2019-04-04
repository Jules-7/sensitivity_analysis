""" Move flights per WTC into folders according to the route number """

import os

wtc = 'L'

flights_dir = '/media/julia/UserData/MC/MC_%s/flight' % wtc
target_dir = '/media/julia/UserData/MC/true_flights/%s' % wtc

for filename in os.listdir(flights_dir):
    route_num = filename[-5]
    os.rename('%s/%s' % (flights_dir, filename), '%s/%s/%s' % (target_dir, route_num, filename))
