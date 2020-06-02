import numpy as np
import argparse

parser = argparse.ArgumentParser()
# parser.add_argument('-zone',help='Time zone difference from UT',type=int,default=0)
parser.add_argument('-table',help='Supply table for calculations for multiple nights')
args = parser.parse_args()

def convert_to_datetime(date_string,time_string):

    date_time = date_string[-4:]+'-'+date_string[:2]+'-'+date_string[3:5]+'T'+time_string

    return np.datetime64(date_time)

print("\n")

if args.table is not None:
    evening_twilight_date,evening_twilight_time,morning_twilight_date,morning_twilight_time,zone = np.loadtxt(args.table,dtype=str,unpack=True)

    ntransits = len(evening_twilight_date)

    for i in range(ntransits):

        q1_start = convert_to_datetime(evening_twilight_date[i],evening_twilight_time[i]) - np.timedelta64(zone[i],'h')
        q4_end = convert_to_datetime(morning_twilight_date[i],morning_twilight_time[i]) - np.timedelta64(zone[i],'h')

        quarter_length = (q4_end-q1_start)/4

        q2_start = q1_start+quarter_length
        q3_start = q2_start+quarter_length
        q4_start = q3_start+quarter_length

        print("Q1 start = %s ; Q2 start = %s ; Q3 start = %s ; Q4 start = %s ; Q4 end = %s \n"%(q1_start,q2_start,q3_start,q4_start,q4_end))
