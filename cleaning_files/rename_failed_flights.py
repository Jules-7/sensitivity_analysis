""" Rename the failed flights listed in the log file per WTC """

import os

wtc = 'L'

flights_dir = '/media/julia/UserData/MC/%s' % wtc
log_file_path = os.path.dirname(os.path.realpath(__file__)) + '/res_log_%s.txt' % wtc
counter = 1

failed_flights = open(log_file_path).readlines()
failed_flights = [value.split()[0] for value in failed_flights]

for dirname in os.listdir(flights_dir):

    for filename in os.listdir(flights_dir + '/' + dirname):
        check_name = filename[:-9]

        if check_name in failed_flights:
            new_filename = '__' + filename
            os.rename('%s/%s' % (flights_dir + '/' + dirname, filename), '%s/%s' % (flights_dir+"/"+dirname, new_filename))
            counter += 1

print(counter, 'failed flights')
