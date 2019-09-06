import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-dur',help='Transit duration in hours',type=float)
parser.add_argument('-mid',help='Transit midpoint in UT',type=str)
parser.add_argument('-zone',help='Time zone difference from UT',type=int)
parser.add_argument('-date',help='Date as "yyyy-mm-dd format"',type=str)
parser.add_argument('-os',help='Start of target observability in UT',type=str)
parser.add_argument('-oe',help='End of target observability in UT',type=str)
parser.add_argument('-table',help='Supply table for calculations if doing for multiple nights')
parser.add_argument('-save',help='save output to table?',action="store_true")
args = parser.parse_args()


def convert_to_datetime(date_string,time_string):

    date_time = date_string[-4:]+'-'+date_string[:2]+'-'+date_string[3:5]+'T'+time_string
    
    return np.datetime64(date_time)


if args.table is not None:
	planet,_,transit_duration, mid_point_date, mid_point_time, observability_start_date, observability_start_time, observability_end_date, observability_end_time,zone = np.loadtxt(args.table,dtype=str,unpack=True)
	
	ntransits = len(np.atleast_1d(transit_duration))
	
	if args.save:
		new_tab = open(args.table[:-4]+"_transit_times_and_obs_durations.txt",'w')
	
	for i in range(ntransits):
		
		obs_length_hours = float(transit_duration[i])*2+0.5
		obs_length_mins = int(np.round(obs_length_hours*60))
		transit_duration_mins = int(float(transit_duration[i])*60)
		
		print('\n##############\n')
		print(planet[i],mid_point_date[i])

		
		print("\nFull (desired) observation duration = %d minutes (%.2f hours) "%(obs_length_mins,obs_length_hours))
		
		mid_point_ut = convert_to_datetime(mid_point_date[i],mid_point_time[i])
		
		start_time_ut = mid_point_ut - np.timedelta64(transit_duration_mins+30,'m')
		
		end_time_ut = mid_point_ut + np.timedelta64(transit_duration_mins,'m')
		
		observability_start = convert_to_datetime(observability_start_date[i], observability_start_time[i])

		observability_end = convert_to_datetime(observability_end_date[i], observability_end_time[i])
		
		
		print("Full observation duration allowed = %s minutes (%.2f hours)"%(observability_end-observability_start+np.timedelta64(30,'m'),(observability_end-observability_start+np.timedelta64(30,'m')).astype(float)/60))
		
		
		if start_time_ut + np.timedelta64(30,'m') < observability_start:
			# print('start time constrained by observability start')
			start_time_ut = observability_start - np.timedelta64(30,'m')
			end_time_ut = start_time_ut + np.timedelta64(obs_length_mins,'m')
			
		if end_time_ut > observability_end:
			# print('end time constrained by observability end')
			end_time_ut = observability_end
		
		if start_time_ut + np.timedelta64(obs_length_mins,'m') < end_time_ut:
			# print('updating end time')
			end_time_ut = start_time_ut + np.timedelta64(obs_length_mins,'m')
			
		if end_time_ut - np.timedelta64(obs_length_mins,'m') < start_time_ut and end_time_ut - np.timedelta64(obs_length_mins,'m') > observability_start:
			start_time_ut = end_time_ut - np.timedelta64(obs_length_mins,'m')
			
		obs_length = end_time_ut - start_time_ut
		 
		mid_point_local = mid_point_ut - np.timedelta64(zone[i],'h')
		start_time_local = start_time_ut - np.timedelta64(zone[i],'h')
		end_time_local = end_time_ut - np.timedelta64(zone[i],'h')
				
		print("\nObs start (UT) = %s ; Obs mid (UT) = %s ; Obs end (UT) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_ut,mid_point_ut,end_time_ut,obs_length,obs_length.astype(float)/60))
		
		print("Obs start (local) = %s ; Obs mid (local) = %s ; Obs end (local) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_local,mid_point_local,end_time_local,obs_length,obs_length.astype(float)/60))
		
		if args.save:
			new_tab.write('\n##############\n')
			new_tab.write('%s ; %s'%(planet[i],mid_point_date[i]))
			new_tab.write("\nFull (desired) observation duration = %d minutes (%.2f hours) "%(obs_length_mins,obs_length_hours))
			new_tab.write("\nFull observation duration allowed = %s minutes (%.2f hours)"%(observability_end-observability_start+np.timedelta64(30,'m'),(observability_end-observability_start+np.timedelta64(30,'m')).astype(float)/60))
			new_tab.write("\nObs start (UT) = %s ; Obs mid (UT) = %s ; Obs end (UT) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_ut,mid_point_ut,end_time_ut,obs_length,obs_length.astype(float)/60))
			new_tab.write("Obs start (local) = %s ; Obs mid (local) = %s ; Obs end (local) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_local,mid_point_local,end_time_local,obs_length,obs_length.astype(float)/60))
		
	if args.save:
		new_tab.close()
	
