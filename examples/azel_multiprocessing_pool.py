from __future__ import print_function    # (at top of module)
# A script to calaculate Rotation Measures at transit for elevations
#   in steps of 5 degrees in elevation
# Parallel processing is used

from multiprocessing import Pool
from RMextract.getRM import getRM
import RMextract.PosTools as PosTools
import numpy 
import time
import math
import os
import future
#from string import split, strip

###########################################################################

# worker function for multi processing
def iono_worker(*args, **kwargs):
    result = getRM(*args, **kwargs) 
    return result

###########################################################################

# generate output report in ASCII format
def output_report(log, stat_pos, results):
 log.write ('RMextract report \n')
 log.write ('columns explanation (zero relative):\n')
 log.write ('0 - sequence number (zero relative)\n')
 log.write ('1 - separator colon\n')
 log.write('2 - value of 0 - valid data/calculation\n')
 log.write('    value of 1 - invalid data/calculation (usually, but not always, because elevation < 0 )\n')
 log.write ('3 - relative time in seconds relative to reference time\n')
 log.write ('4 - time step between sequence calculations (default is 300 sec)\n')
 log.write ('5 - elevation in degrees\n')
 log.write ('6 - azimuth in degrees\n')
 log.write ('7 - TEC (in tec units) in the current azimuth/elevation direction\n')
 log.write ('8 - Rotation Measure (radians/m^2) in the current azimuth/elevation direction\n')
 log.write ('9 - correction factor to convert STEC to value at the zenith\n')
 log.write ('     \n')

 for i in range(len(results)):
  result = results[str(i)]  
  stat_pos = result['stat_pos']
  stat_name = result['station_names']
  shape  = result['times'].shape
  timerange=[result['times'][0],result['times'][shape[0]-1]]
  reference_time=result['reference_time'] 
  timegrid=result['times']
  timerange=[result['times'][0],result['times'][shape[0]-1]]
  reference_time=result['reference_time']
  str_start_time=PosTools.obtain_observation_year_month_day_hms(reference_time)
  log.write ('start and end times %s %s \n' % (timerange[0], timerange[1]))
  log.write ('reference time for rel_time=0: year month day hr min sec %s %s %s %s %s %s \n' % str_start_time)

  k = 0
  for key in result['station_names']:
    seq_no = 0
    log.write ('az, el data for station  at position %s \n' % (stat_pos))
    log.write ('\nseq  rel_time time_width El  Az    STEC           RM (rad/m2)     VTEC factor  \n')

    for i in range (timegrid.shape[0]):
      el = math.degrees(result['elev'][key][i])
      # handle round-off error near zero
      if el < 1.0e-13:
        el = 0.0
      az = math.degrees(result['azimuth'][key][i])
      stec =result['STEC'][key][i]
      rm = result['RM'][key][i]
      vtec_factor = 1.0 / result['AirMass'][key][i]
      rel_time = timegrid[i] - reference_time
      if i  == 0:
        time_width = reference_time - timegrid[i] 
      else:
        time_width = timegrid[i] - timegrid[i-1]
      if el < 0 :
        ok = 1
        stec = 0.0
        rm = 0.0
        vtec_factor = 1.0
      else:
        ok = 0
      log.write("%s : %s %s %s %s %s %s %s %s\n" % (seq_no, ok, rel_time, time_width, el, az, stec, rm, vtec_factor))
      seq_no = seq_no + 1
    log.write (' \n')

#########################################################################
#########################################################################

# Start of actual processing

startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
process_start = time.time()
print ("Start at %s" % startime)

# find number of cores available for multiprocessing
num_processors = 1
if num_processors <= 2:
    try:
      import multiprocessing
      processors =  multiprocessing.cpu_count()
      if processors > num_processors:
        num_processors = processors
        print ('*** setting number of processors to',num_processors)
    except:
      pass

OBJECT="DRAO_transit_az_el"
START_TIME="2017/11/09 00:00:00"
END_TIME="2017/11/09 23:59:59"
TIME_STEP=300.0
TIME_OFFSET=120.0
use_azel = True
out_file='RMextract_report_' + OBJECT
multi = False
multi = True

DRAO_position = [[-2059666.914945, -3621225.951266, 4814123.619458]]
stat_pos = DRAO_position
az = math.radians(180.0)
max_el = math.radians(90.0)
delta_el = math.radians(5.0)


# required parameters for RMextract
args = (None,'ftp://ftp.aiub.unibe.ch/CODE/','CODG','IONEXdata/',0,0,use_azel,-1000,None)


results_dict = {}

if multi:
  el = 0.0
  pointing=numpy.array([az, el])
# optional keywords for RMextract
  kw = {}       
  kw['radec']=pointing
  kw['object']=OBJECT
  kw['start_time']=START_TIME
  kw['end_time']=END_TIME 
  kw['timestep']=TIME_STEP
  kw['useEMM']=True
  kw['TIME_OFFSET']=TIME_OFFSET
  kw['stat_positions']= [stat_pos[0]] 
  kw['stat_names']=  ['station-number_' + str(0)]
# we need to make a single non-parallel call to getRM in order to get IONEX
# data if it has not been previously collected
  res = getRM(*args, **kw) 

  station = res['station_names'][0]
  index = station.split('_')
  key = index[1]
  results_dict[key] = res

  results = []
  P = Pool(processes=num_processors)
  start_iter = 1
else:
  start_iter = 0
  el = math.radians(-5.0)

for i in range(start_iter,50):
  el = el + delta_el
  if el > max_el:
    delta_el = -1.0* delta_el
    el = max_el + delta_el
    az = 0.0
  if el < 0.0:
      break

  pointing=numpy.array([az, el])
  print ('pointing', pointing)

  # for some wierd reason I have to instantiate this dict`
  # every time I go around this loop or else one task
  # gets dropped somewhere in Pool-land
  kw = {}       
  kw['object']=OBJECT
  kw['start_time']=START_TIME
  kw['end_time']=END_TIME 
  kw['timestep']=TIME_STEP
  kw['useEMM']=True
  kw['TIME_OFFSET']=TIME_OFFSET
  kw['radec']=pointing
  kw['stat_positions']= [stat_pos[0]] 
  kw['stat_names']=  ['station-number_' + str(i)]
  if multi:
    results.append(P.apply_async(iono_worker, args, kw))
  else:
    res = getRM(*args, **kw) 
    station = res['station_names'][0]
    index = station.split('_')
    key = index[1]
    results_dict[key] = res

if multi:
  print ('number of pool  tasks ', len(results))
  P.close()
  P.join()
  for s in results:
    res = s.get(timeout = 1 )
    station = res['station_names'][0]
    index = station.split('_')
    key = index[1]
    results_dict[key] = res

if os.path.exists(out_file):
  os.remove(out_file)
log = open(out_file, 'a')
output_report(log, stat_pos, results_dict)
log.close()

endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print ("End at %s" % endtime)
process_end = time.time()
duration = (process_end - process_start)/60.0
print ("Total run time: %7.2f minutes" % duration)
