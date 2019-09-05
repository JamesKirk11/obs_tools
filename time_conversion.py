import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-dur',help='Transit duration in hours',type=float)
parser.add_argument('-mid',help='Transit midpoint in UT',type=str)
parser.add_argument('-zone',help='Time zone difference from UT',type=int)
parser.add_argument('-date',help='Date as "yyyy-mm-dd format"',type=str)
parser.add_argument('-os',help='Start of target observability in UT',type=str)
parser.add_argument('-oe',help='End of target observability in UT',type=str)
args = parser.parse_args()


obs_length = args.dur*2+0.5

print("\nFull (desired) observation duration = %d minutes (%.2f hours) "%(obs_length*60,obs_length))

mid_point_ut = np.datetime64("T".join([args.date,args.mid]))
start_time_ut = mid_point_ut - np.timedelta64(int(np.round((args.dur+0.5)*60)),'m')

# print(start_time_ut)

end_time_ut = mid_point_ut + np.timedelta64(int(np.round((args.dur)*60)),'m')

# print(end_time_ut)

try:
	observability_start = np.datetime64("T".join([args.date,args.os]))
except:
	observability_start = np.datetime64(args.os)

try:
	observability_end = np.datetime64("T".join([args.date,args.oe]))
except:
	observability_end = np.datetime64(args.oe)


print("Full observation duration allowed = %s"%(observability_end-observability_start+np.timedelta64(30,'m')))


if start_time_ut + np.timedelta64(30,'m') < observability_start:
	start_time_ut = observability_start - np.timedelta64(30,'m')
	end_time_ut = start_time_ut + np.timedelta64(int(np.round(obs_length*60)),'m')
	
if end_time_ut > observability_end:
	end_time_ut = observability_end

if start_time_ut + np.timedelta64(int(np.round(obs_length*60)),'m') < end_time_ut:
	end_time_ut = start_time_ut + np.timedelta64(int(np.round(obs_length*60)),'m')
	
if end_time_ut - np.timedelta64(int(np.round(obs_length*60)),'m') < start_time_ut and end_time_ut - np.timedelta64(int(np.round(obs_length*60)),'m') > observability_start:
	start_time_ut = end_time_ut - np.timedelta64(int(np.round(obs_length*60)),'m')
	

obs_length = end_time_ut - start_time_ut
 
mid_point_local = mid_point_ut - np.timedelta64(args.zone,'h')
start_time_local = start_time_ut - np.timedelta64(args.zone,'h')
end_time_local = end_time_ut - np.timedelta64(args.zone,'h')

 
print("\nObs start (UT) = %s ; Obs mid (UT) = %s ; Obs end (UT) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_ut,mid_point_ut,end_time_ut,obs_length,obs_length.astype(float)/60))

print("Obs start (local) = %s ; Obs mid (local) = %s ; Obs end (local) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_local,mid_point_local,end_time_local,obs_length,obs_length.astype(float)/60))
