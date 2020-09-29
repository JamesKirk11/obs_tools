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
parser.add_argument('-use_quarters',help='split the observations neatly into quarters for Keck obs?',action="store_true")
parser.add_argument('-tel',help='Which telescope are we observing with? If Keck (K), baseline is transit +2 (baseline) +0.5 (extra duration) +0.5 (setup) = +3 hours. If Magellan (M), baseline is transit +2 (baseline) +0.25 (acquisition) +0.5 (standard) = +2.75 hours')
parser.add_argument('-hard',help="Use a hard limit on the out of transit baseline? i.e. don't maximise out of transit baseline by having imbalanced pre- and post-transit. This is most time economical",action="store_true")
args = parser.parse_args()


def convert_to_datetime(date_string,time_string):

    date_time = date_string[-4:]+'-'+date_string[:2]+'-'+date_string[3:5]+'T'+time_string

    return np.datetime64(date_time)


if args.table is not None:
    if args.use_quarters:
        planet,_,transit_duration, mid_point_date, mid_point_time, observability_start_date, observability_start_time, observability_end_date, observability_end_time,zone,evening_twilight_date,evening_twilight_time,morning_twilight_date,morning_twilight_time = np.loadtxt(args.table,dtype=str,unpack=True)

    else:
        planet,_,transit_duration, mid_point_date, mid_point_time, observability_start_date, observability_start_time, observability_end_date, observability_end_time,zone = np.loadtxt(args.table,dtype=str,unpack=True)

    ntransits = len(np.atleast_1d(transit_duration))

    if args.save:
        new_tab = open(args.table[:-4]+"_transit_times_and_obs_durations.txt",'w')

    for i in range(ntransits):

        if args.tel.lower() == 'k':
            obs_length_hours = float(transit_duration[i])+3#2+0.5

        elif args.tel.lower() == 'm':
            obs_length_hours = float(transit_duration[i])+2.75

        else:
            raise NameError("have to define -tel as 'm' (Magellan) or 'k' (Keck)")

        obs_length_mins = int(np.round(obs_length_hours*60))
        transit_duration_mins = int(np.round(float(transit_duration[i])*60))

        if args.use_quarters:
            q1_start = convert_to_datetime(evening_twilight_date[i],evening_twilight_time[i]) - np.timedelta64(zone[i],'h')
            q4_end = convert_to_datetime(morning_twilight_date[i],morning_twilight_time[i]) - np.timedelta64(zone[i],'h')

            quarter_length = (q4_end-q1_start)/4

            q2_start = q1_start+quarter_length
            q3_start = q2_start+quarter_length
            q4_start = q3_start+quarter_length




        print('\n##############\n')
        print(planet[i],mid_point_date[i])


        print("\nFull (desired) observation duration = %d minutes (%.2f hours) "%(obs_length_mins,obs_length_hours))

        mid_point_ut = convert_to_datetime(mid_point_date[i],mid_point_time[i])
        transit_start_ut = mid_point_ut - np.timedelta64(int(np.round(transit_duration_mins/2)),'m')
        transit_end_ut = mid_point_ut + np.timedelta64(int(np.round(transit_duration_mins/2)),'m')

        if args.tel.lower() == 'k':
            # start_time_ut = mid_point_ut - np.timedelta64(int(np.round(transit_duration_mins/2))+90,'m')
            # end_time_ut = mid_point_ut + np.timedelta64(int(np.round(transit_duration_mins/2+60)),'m')
            start_time_ut = mid_point_ut - np.timedelta64(int(np.round(transit_duration_mins/2)+15+60+30),'m')
            end_time_ut = mid_point_ut + np.timedelta64(int(np.round(transit_duration_mins/2)+15+60),'m')

        if args.tel.lower() == 'm':
            # start_time_ut = mid_point_ut - np.timedelta64(transit_duration_mins+30,'m')
            # end_time_ut = mid_point_ut + np.timedelta64(transit_duration_mins,'m')
            start_time_ut = mid_point_ut - np.timedelta64(int(np.round(transit_duration_mins/2)+15+60+30),'m')
            end_time_ut = mid_point_ut + np.timedelta64(int(np.round(transit_duration_mins/2)+60),'m')


        observability_start = convert_to_datetime(observability_start_date[i], observability_start_time[i])
        observability_end = convert_to_datetime(observability_end_date[i], observability_end_time[i])

        max_obs_allowed = observability_end-observability_start+np.timedelta64(30,'m')
        max_obs_allowed_hours = max_obs_allowed.astype(float)/60

        print("Full observation duration allowed = %s (%.2f hours)"%(max_obs_allowed,max_obs_allowed_hours))

        if (max_obs_allowed_hours - obs_length_hours) < -1:
            print("WARNING!! Max allowed observation time 1 hour less than desired observation time!")

        if not args.hard:
            if start_time_ut + np.timedelta64(30,'m') < observability_start:
                start_time_ut = observability_start - np.timedelta64(30,'m')
                end_time_ut = start_time_ut + np.timedelta64(obs_length_mins,'m')

            if end_time_ut > observability_end:
                end_time_ut = observability_end

            if start_time_ut + np.timedelta64(obs_length_mins,'m') < end_time_ut:
                end_time_ut = start_time_ut + np.timedelta64(obs_length_mins,'m')

            if end_time_ut - np.timedelta64(obs_length_mins,'m') < start_time_ut: #and
                if end_time_ut - np.timedelta64(obs_length_mins,'m') > observability_start:
                    start_time_ut = end_time_ut - np.timedelta64(obs_length_mins,'m')
                if end_time_ut - np.timedelta64(obs_length_mins,'m') < observability_start - np.timedelta64(30,'m'):
                    start_time_ut = observability_start - np.timedelta64(30,'m')

        if args.hard:
            if start_time_ut + np.timedelta64(30,'m') < observability_start:
                start_time_ut = observability_start - np.timedelta64(30,'m')
            if end_time_ut > observability_end:
                end_time_ut = observability_end

        obs_length = end_time_ut - start_time_ut

        mid_point_local = mid_point_ut - np.timedelta64(zone[i],'h')
        start_time_local = start_time_ut - np.timedelta64(zone[i],'h')
        end_time_local = end_time_ut - np.timedelta64(zone[i],'h')

        pre_baseline = (transit_start_ut-(start_time_ut + np.timedelta64(30,'m')))/np.timedelta64(1, 'm')
        print(pre_baseline)
        if pre_baseline < 0:
            pre_baseline = 0

        post_baseline = (end_time_ut-transit_end_ut)/np.timedelta64(1, 'm')
        if post_baseline < 0:
            post_baseline = 0

        print("\nObs start (UT) = %s ; Obs mid (UT) = %s ; Obs end (UT) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_ut,mid_point_ut,end_time_ut,obs_length,obs_length.astype(float)/60))

        print("Obs start (local) = %s ; Obs mid (local) = %s ; Obs end (local) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_local,mid_point_local,end_time_local,obs_length,obs_length.astype(float)/60))

        print("Pre transit baseline = %s mins ; Post transit baseline = %s mins \n"%(pre_baseline,post_baseline))

        if args.use_quarters:
            print("\n## Quarters (local time) \n")
            quarters = np.array([q1_start,q2_start,q3_start,q4_start,q4_end])

            transit_start_quarter_idx = np.where((transit_start_ut - np.timedelta64(zone[i],'h')) > quarters)[0].max()
            transit_start_quarter = quarters[transit_start_quarter_idx]

            transit_end_quarter_idx = np.where((transit_end_ut - np.timedelta64(zone[i],'h')) > quarters)[0].max()
            transit_end_quarter = quarters[transit_end_quarter_idx]
            transit_end_quarter_end = quarters[transit_end_quarter_idx+1] # end of the quarter that the transit ends in

            print("Q1 start = %s ; Q2 start = %s ; Q3 start = %s ; Q4 start = %s ; Q4 end = %s \n"%(q1_start,q2_start,q3_start,q4_start,q4_end))
            print("Transit starts at %s in Q%d ; transit ends at %s in Q%d \n"%(transit_start_ut - np.timedelta64(zone[i],'h'),transit_start_quarter_idx+1,transit_end_ut - np.timedelta64(zone[i],'h'),transit_end_quarter_idx+1))
            print("%d mins before transit in Q%d ; %d mins after transit in Q%d \n"%((transit_start_ut- np.timedelta64(zone[i],'h')-transit_start_quarter)/np.timedelta64(1, 'm'),transit_start_quarter_idx+1,\
            (transit_end_quarter_end-(transit_end_ut- np.timedelta64(zone[i],'h')))/np.timedelta64(1, 'm'),transit_end_quarter_idx+1))

            # calculate amount of time in quarter before and after the target is observable
            quarter_dead_time_pre = (observability_start- np.timedelta64(zone[i],'h')-transit_start_quarter)/np.timedelta64(1, 'm')
            quarter_dead_time_post = (transit_end_quarter_end - (observability_end- np.timedelta64(zone[i],'h')))/np.timedelta64(1, 'm')

            if quarter_dead_time_pre < 0:
                quarter_dead_time_pre = 0

            if quarter_dead_time_post < 0:
                quarter_dead_time_post = 0

            print("%d mins deadtime pre-transit in Q%d ; %d mins deadtime post-transit in Q%d \n"%(quarter_dead_time_pre,transit_start_quarter_idx+1,\
            quarter_dead_time_post,transit_end_quarter_idx+1))
            # raise SystemExit
            # for i in q
            # transit_start_quarter = transit_start_ut >
            # print("Transit starts in %s"%(transit_start_ut>))

        if args.save:
            new_tab.write('\n##############\n')
            new_tab.write('%s ; %s'%(planet[i],mid_point_date[i]))
            new_tab.write("\nFull (desired) observation duration = %d minutes (%.2f hours) "%(obs_length_mins,obs_length_hours))
            new_tab.write("\nFull observation duration allowed = %s (%.2f hours)"%(observability_end-observability_start+np.timedelta64(30,'m'),(observability_end-observability_start+np.timedelta64(30,'m')).astype(float)/60))
            new_tab.write("\nObs start (UT) = %s ; Obs mid (UT) = %s ; Obs end (UT) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_ut,mid_point_ut,end_time_ut,obs_length,obs_length.astype(float)/60))
            new_tab.write("Obs start (local) = %s ; Obs mid (local) = %s ; Obs end (local) = %s ; Obs dur = %s (%.2f hours) \n"%(start_time_local,mid_point_local,end_time_local,obs_length,obs_length.astype(float)/60))
            new_tab.write("Pre transit baseline = %s mins ; Post transit baseline = %s mins \n"%(pre_baseline,post_baseline))

            if args.use_quarters:
                new_tab.write("\n## Quarters (local time) \n")
                new_tab.write("Q1 start = %s ; Q2 start = %s ; Q3 start = %s ; Q4 start = %s ; Q4 end = %s \n"%(q1_start,q2_start,q3_start,q4_start,q4_end))
                new_tab.write("Transit starts at %s in Q%d ; transit ends at %s in Q%d \n"%(transit_start_ut- np.timedelta64(zone[i],'h'),transit_start_quarter_idx+1,transit_end_ut- np.timedelta64(zone[i],'h'),transit_end_quarter_idx+1))
                new_tab.write("%d mins before transit in Q%d ; %d mins after transit in Q%d \n"%((transit_start_ut- np.timedelta64(zone[i],'h')-transit_start_quarter)/np.timedelta64(1, 'm'),transit_start_quarter_idx+1,\
                (transit_end_quarter_end-(transit_end_ut- np.timedelta64(zone[i],'h')))/np.timedelta64(1, 'm'),transit_end_quarter_idx+1))
                new_tab.write("%d mins deadtime pre-transit in Q%d ; %d mins deadtime post-transit in Q%d \n"%(quarter_dead_time_pre,transit_start_quarter_idx+1,\
                quarter_dead_time_post,transit_end_quarter_idx+1))


            if (max_obs_allowed_hours - obs_length.astype(float)/60) < -1:
                new_tab.write("WARNING!! Max allowed observation time 1 hour less than desired observation time!\n")

    if args.save:
        new_tab.close()
