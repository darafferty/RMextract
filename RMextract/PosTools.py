#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 11:48:10 2018

@author: mevius
"""
import os
import numpy as np
import logging


HAS_PYRAP = True
try:
  from pyrap import tables as tab
  import pyrap.quanta as qu
  from pyrap.measures import measures
  logging.info ('pyrap will be used to compute positions')
except:
  logging.info ('We will need PyEphem to perform calculations!')
  logging.info ('the accuracy of results might decease a bit')
  HAS_PYRAP = False

HAS_EPHEM = True
try:
  import ephem
  if not HAS_PYRAP:
    logging.info ('PyEphem will be used to perform calculations!')
except:
  HAS_EPHEM = False


# height of thin layer ionosphere
ION_HEIGHT=450.e3

# some constants
R_earth=6364.62e3
earth_ellipsoid_a = 6378137.0
earth_ellipsoid_a2 = earth_ellipsoid_a*earth_ellipsoid_a
earth_ellipsoid_b = 6356752.3142
earth_ellipsoid_b2 = earth_ellipsoid_b*earth_ellipsoid_b
earth_ellipsoid_e2 = (earth_ellipsoid_a2 - earth_ellipsoid_b2) / earth_ellipsoid_a2
# default station position, center of LOFAR
posCS002=[3826577.1095, 461022.900196, 5064892.758]
########################################
sm_a = 6378137.0
invf = 298.257223563
f = 1.0 / invf

# The following two functions were provided by Emil Lenc
# Convert latitude (radians), longitude (radians) and elevation (metres) 
# to ITRF XYZ and vice versa
def WGS84ToITRF(lats, lons, h): 
    '''WGS-84 to ITRF (input in radians).
    Convert latitude (radians), longitude (radians) and elevation (metres) to ITRF XYZ
     Args:
        lats (np.array) : angles in radians
        lons (np.array) : angles in radians
        h(float) : elevation (meters)
    Returns:
        Tuple[np.array, np.array, np.array] : arrays 
        of ITRF x, y, z positions    
    '''
    SINK = np.sin(lats)
    COSK = np.cos(lats)
    e2 = 2.0 * f - f * f
    v = sm_a / np.sqrt(1.0 - e2 * SINK * SINK)
    x = (v + h) * COSK * np.cos(lons)
    y = (v + h) * COSK * np.sin(lons)
    z = ((1 - e2) * v + h) * SINK
    return x, y, z

def ITRFToWGS84(x, y, z):
    '''ITRF to WGS-84 (input in meters).
    Convert ITRF XYZ to latitude (radians), longitude (radians) and elevation (metres) 
     Args:
        x (np.array) : x positions 
        y (np.array) : y positions
        z (np.array) : z positions
    Returns:
        Tuple[np.array, np.array, np.array] : arrays of 
          latitude (radians), longitude (radians)
          and elevation (metres)  
    '''
    e2 = 2.0 * f - f * f
    E = e2 / (1.0 - e2)
    b = sm_a * (1.0 - f)
    p = np.sqrt(x * x + y * y)
    q = np.atan2(z * sm_a, (p * b))
    lat = np.atan2((z + E * b * np.sin(q) * np.sin(q) \
                      * np.sin(q)), (p - e2 * sm_a * np.cos(q)\
                              *np.cos(q) * np.cos(q)))
    v = sm_a / np.sqrt(1.0 - e2 * np.sin(lat) * np.sin(lat))
    lon = np.atan2(y, x)
    h = (p / np.cos(lat)) - v
    #lat = np.degrees(lat)
    #lon = np.degrees(lon)
    return lat, lon, h

def GeodeticToGeocentricLat(geodetic_lat, height):
    '''convert geodetic latitude to geocentric latitude.
    Args:
        geodetic_lat (np.array) : geodetic latitude in radians
        height (float) : height
    Returns:
        np.array : geocentric latitudes in radians
        '''
    l_sin = np.sin(geodetic_lat)
    e2 = 2.0 * f - f * f
    div = np.sqrt(1 - e2 * l_sin**2)
    rn =  sm_a / div 
    rn_div = rn + height
    ratio = 1 - e2 * rn / rn_div
    tan_geocentric_lat = ratio * np.tan(geodetic_lat) 
    geocentric_lat = np.atan(tan_geocentric_lat)
    return geocentric_lat

########################################
# see http://stackoverflow.com/questions/15954978/ecef-from-azimuth-elevation-range-and-observer-lat-lon-alt

# input expected in degrees here
def aer2ecef(azimuthDeg, elevationDeg, slantRange, obs_lat, obs_lon, obs_alt):
    '''Convert az,el, range to ECEF positions.
    
    see http://stackoverflow.com/questions/15954978/ecef-from-azimuth-elevation-range-and-observer-lat-lon-alt.
    
    Args:
       azimuthDeg (np.array) : azimuths (degrees) 
       elevationDeg (np.array) : elevations (degrees)
       slantRange (float) : range
       obs_lat (float) : latitude of observer (degrees)
       obs_lon (float) : longitude of observer (degrees)
       obs_alt(float) : altitude of observer (metres)
    Returns : 
        Tuple[np.array,np.array,np.array] : latitude, longitude, height
    '''
    obs_lat_r = np.radians(obs_lat)
    obs_lon_r = np.radians(obs_lon)
    sitex, sitey, sitez = WGS84ToITRF(obs_lat_r,obs_lon_r,obs_alt)

    #some needed calculations
    slat = np.sin(obs_lat_r)
    slon = np.sin(obs_lon_r)
    clat = np.cos(obs_lat_r)
    clon = np.cos(obs_lon_r)

    azRad = np.radians(azimuthDeg)
    elRad = np.radians(elevationDeg)

    # az,el,range to sez convertion
    south  = -slantRange * np.cos(elRad) * np.cos(azRad)
    east   =  slantRange * np.cos(elRad) * np.sin(azRad)
    zenith =  slantRange * np.sin(elRad)


    x = ( slat * clon * south) + (-slon * east) + (clat * clon * zenith) + sitex
    y = ( slat * slon * south) + ( clon * east) + (clat * slon * zenith) + sitey
    z = (-clat *        south) + ( slat * zenith) + sitez
    lat,lon,h = ITRFToWGS84(x, y, z)
    return lat, lon, h  # data are in units of degrees

################ from JMA_tools in ALBUS ########################
def get_day_of_year(year, month, day):
    """Get the day-of-year integer from the year/month/day
    
    Args:
        year (int) : YYYY
        month (int) : MM starting from January == 1
        day (int) :  DD starting from 1 == 1
    """
    day_of_year_list = [0,31,59,90,120,151,181,212,243,273,304,334]
    doy = day_of_year_list[month-1] + day
    if(month>2):
        if((year&0x3)==0):
            if((year % 100 != 0) or (year % 400 == 0)):
                doy = doy+1
    return doy
################ from JMA_tools in ALBUS ########################
def get_ymdf_from_JD(JD):
    """get the year, month, day, day_fraction from an MJD 
    
    Taken from _Astronomical Algorithms_, Meeus, 1991
    Args:
        JD (float) : MJD
    Returns:
        Tuple[int,int,int,float] : year, month, day, day_fraction
    """
    JD2 = JD + 0.5
    Z = int(JD2)
    F = JD2 - Z
    A = Z
    if(Z>= 2299161):
        alpha = int((Z-1867216.25)/36524.25)
        A = Z + 1 + alpha - int(alpha/4)
    B = A + 1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)
    day_total = B - D - int(30.6001*E) + F
    month = E - 1
    if(E >= 14): month = E - 13
    year = C - 4716
    if(month <= 2): year += 1
    day = int(day_total)
    day_fraction = day_total - day
    return year, month, day, day_fraction

################################################################################
def get_hms_from_frac(day_fraction):
    """get hours, minues, seconds from a fractional day.

Does not worry about leap seconds.
"""
    h = day_fraction * 24.0
    hour = int(h+2E-13)
    m = (h - hour) * 60.0
    minute = int(m+1E-11)
    second = (m - minute) * 60.0
    return hour, minute, second


################################################################################
def get_ymdh_from_JD(JD):
    """get hours, minues, seconds from a fractional day.
"""
    year, month, day, day_fraction = get_ymdf_from_JD(JD)
    hour, minute, second = get_hms_from_frac(day_fraction)
    return year, month, day, hour, minute, second

################################################################################
def obtain_observation_year_month_day_fraction(start_time):
    julian_day = (start_time / 86400.0) + 2400000.5
    result = get_ymdf_from_JD(julian_day)
    return result

################################################################################
# Get the day of year from the Year, month, day for the start of observations
def obtain_observation_year_month_day_hms(start_time):
    if HAS_PYRAP:
      #logging.info ("getting time", str(start_time)+'s')
      date_list = qu.quantity(str(start_time)+'s').formatted("YMD").split("/")
      year = int(date_list[0])
      month = int(date_list[1])
      day = int(date_list[2])
      time_list=date_list[3].split(":")
      return (year, month, day,int(time_list[0]),int(time_list[1]),float(time_list[2]))
    else:
      julian_day = (start_time / 86400.0) + 2400000.5
      year, month, day, hour, minute, second  = get_ymdh_from_JD(julian_day)
      return (year, month, day,hour, minute, second)





def getMSinfo(MS):
    '''Get metadata fromMeasurement set
    
    Args: 
        MS (string) : path to the MS
    Returns: 
        Tuple[list,float,np.array,np.array,np.array] : 
            timerange,timestep,pointing (ra,dec in radians), stations,
            station positions (x,y, in meters)''' 
    if not HAS_PYRAP:
        logging.info ("Install pyrap to be able to extract info from MS")
        return

        
    if os.path.isdir(MS):
        myMS=tab.table(MS)
    else:
        logging.error("Do not understand the format of MS  %s "%MS)
        return
    timerange=[np.amin(myMS.getcol('TIME_CENTROID')),
               np.amax(myMS.getcol('TIME_CENTROID'))]
    timestep=myMS.getcell('INTERVAL',0)
    
    pointing= tab.table(myMS.getkeyword('FIELD')).getcell('PHASE_DIR',0)    
    stations = tab.table(myMS.getkeyword('ANTENNA')).getcol('NAME')
    station_pos = tab.table(myMS.getkeyword('ANTENNA')).getcol('POSITION')

    return timerange,timestep,pointing.flatten(),stations,station_pos


def getPPsimple(height=[ION_HEIGHT,], mPosition=[0.,0.,0.],
                direction=[0.,0.,0.]):
    '''get ITRF xyz position of piercepoints 
    
    for antenna position mPosition in m, direction ITRF in m on unit sphere 
    and for array of heights, assuming a spherical Earth
    
    Args:
        height (list) : array of altitudes in meters
        mPosition (list) : x,y,z position of stations (meters)
        directions (list) : x,y,z of unit vector of ITRF direction
    Returns:
        Tuple[np.array,np.array] : array of piercepoints (nr heights x xyz)
        array of airmass factors
    '''
    height=np.array(height)
    stX=mPosition[0]
    stY=mPosition[1]
    stZ=mPosition[2]
    x=np.divide(stX,(R_earth+height))
    y=np.divide(stY,(R_earth+height))
    z=np.divide(stZ,(R_earth+height))
    
    c = x*x + y*y + z*z - 1.0
   
    dx=np.divide(direction[0],(R_earth+height))
    dy=np.divide(direction[1],(R_earth+height))
    dz=np.divide(direction[2],(R_earth+height))

    a = dx*dx + dy*dy + dz*dz
    b = x*dx + y*dy  + z*dz

    alpha = (-b + np.sqrt(b*b - a*c))/a


    pp=np.zeros(height.shape+(3,))
    pp[:,0]=stX+alpha*direction[0]
    pp[:,1]=stY+alpha*direction[1]
    pp[:,2]=stZ+alpha*direction[2]

    am=np.divide(1.,pp[:,0]*dx+pp[:,1]*dy+pp[:,2]*dz)
    return pp,am

def getPPsimpleAngle(height=[ION_HEIGHT,], mPosition=[0.,0.,0.],
                     direction=[0.,0.,0.]):
    '''get (lon,lat,h values) of piercepoints 
    
    for antenna position mPosition in m, direction ITRF in m on unit sphere
    and for array of heights, assuming a spherical Earth
    Args:
        height (list) : array of altitudes in meters
        mPosition (list) : x,y,z position of stations (meters)
        directions (list) : x,y,z of unit vector of ITRF direction
    Returns:
        Tuple[np.array,np.array] : array of piercepoints (nr heights x 
        lon,lat,h) array of airmass factors
    
    '''
    pp,am = getPPsimple(height, mPosition, direction)
    ppl=np.zeros(height.shape+(3,))
    ppl[:,0]=np.arctan2(pp[:,1],pp[:,0])
    ppl[:,1]=np.arctan2(pp[:,2],np.sqrt(pp[:,0]*pp[:,0]+pp[:,1]*pp[:,1]))
    ppl[:,2]=height

    return ppl,am
    

def getPP(h=ION_HEIGHT,mPosition=[0.,0.,0.],direction=[0.,0.,0.]):
    '''get ITRF xyz position of piercepoint
    
    for antenna position mPosition in m, direction ITRF in m on unit sphere 
    and for single height, taking into account elipsoidal shape of Earth    
    Args:
        height (float) : altitude in meters
        mPosition (list) : x,y,z position of stations (meters)
        directions (list) : x,y,z of unit vector of ITRF direction
    Returns:
        Tuple[np.array,float] : array of piercepoint  xyz, airmass factor
    '''
   
    stationX = mPosition[0]
    stationY = mPosition[1]
    stationZ = mPosition[2]

    ion_ellipsoid_a = earth_ellipsoid_a + h
    ion_ellipsoid_a2_inv = 1.0 / (ion_ellipsoid_a * ion_ellipsoid_a)
    ion_ellipsoid_b = earth_ellipsoid_b + h
    ion_ellipsoid_b2_inv = 1.0 / (ion_ellipsoid_b * ion_ellipsoid_b)
    
    x = stationX/ion_ellipsoid_a
    y = stationY/ion_ellipsoid_a
    z = stationZ/ion_ellipsoid_b
    c = x*x + y*y + z*z - 1.0

    dx = direction [0]/ ion_ellipsoid_a
    dy = direction [1] / ion_ellipsoid_a
    dz = direction [2] / ion_ellipsoid_b

    a = dx*dx + dy*dy + dz*dz
    b = x*dx + y*dy  + z*dz
    alpha = (-b + np.sqrt(b*b - a*c))/a
    pp_x = stationX + alpha*direction[0]
    pp_y = stationY + alpha*direction[1]
    pp_z = stationZ + alpha*direction[2]

    normal_x = pp_x * ion_ellipsoid_a2_inv
    normal_y = pp_y * ion_ellipsoid_a2_inv
    normal_z = pp_z * ion_ellipsoid_b2_inv
    norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z
    norm_normal = np.sqrt(norm_normal2)
    sin_lat2 = normal_z*normal_z / norm_normal2

 
    g = 1.0 - earth_ellipsoid_e2*sin_lat2
    sqrt_g = np.sqrt(g)

    M = earth_ellipsoid_b2 / ( earth_ellipsoid_a * g * sqrt_g )
    N = earth_ellipsoid_a / sqrt_g

    local_ion_ellipsoid_e2 = (M-N) / ((M+h)*sin_lat2 - N - h)
    local_ion_ellipsoid_a = (N+h) * np.sqrt(
            1.0 - local_ion_ellipsoid_e2*sin_lat2)
    local_ion_ellipsoid_b = local_ion_ellipsoid_a*np.sqrt(
            1.0 - local_ion_ellipsoid_e2)

    z_offset = ((1.0-earth_ellipsoid_e2)*N + h \
                - (1.0-local_ion_ellipsoid_e2)*(N+h)) * np.sqrt(sin_lat2)

    x1 = stationX/local_ion_ellipsoid_a
    y1 = stationY/local_ion_ellipsoid_a
    z1 = (stationZ-z_offset)/local_ion_ellipsoid_b
    c1 = x1*x1 + y1*y1 + z1*z1 - 1.0

    dx = direction[0] / local_ion_ellipsoid_a
    dy = direction[1] / local_ion_ellipsoid_a
    dz = direction[2] / local_ion_ellipsoid_b
    a = dx*dx + dy*dy + dz*dz
    b = x1*dx + y1*dy  + z1*dz
    alpha = (-b + np.sqrt(b*b - a*c1))/a

    pp_x = stationX + alpha*direction[0]
    pp_y = stationY + alpha*direction[1]
    pp_z = stationZ + alpha*direction[2]

    normal_x = pp_x / (local_ion_ellipsoid_a * local_ion_ellipsoid_a)
    normal_y = pp_y / (local_ion_ellipsoid_a * local_ion_ellipsoid_a)
    normal_z = (pp_z-z_offset) / (local_ion_ellipsoid_b * local_ion_ellipsoid_b)

    norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z
    norm_normal = np.sqrt(norm_normal2)
    
    pp_airmass = norm_normal / (direction[0]*normal_x 
                                + direction[1]*normal_y 
                                + direction[2]*normal_z)

    return (pp_x,pp_y,pp_z,pp_airmass)


def getLonLat(pos):
    '''converts ITRF pos in xyz to lon lat
    '''
    me=measures()
    a=me.measure(me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m'),"ITRF")
    return (a['m0']['value'],a['m1']['value'])



def getLonLatStation(az=0,el=0,pos=posCS002):
    '''converts local station direction to ITRF lon/lat'''
    if not isinstance(az,str):
        az=str(az)+'rad'
    if not isinstance(el,str):
        el=str(el)+'rad'
    me=measures()
    me.do_frame(me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m'))
    direction=me.direction("AZELGEO",az,el)
    return me.measure(direction,"ITRF")


def radec2azel(ra,dec,time, pos):
    '''convert ra,dec toazimuth elevation'''
    me=measures()
    if type(ra)!=str:
        ra=str(ra)+'rad'
    if type(dec)!=str:
        dec=str(dec)+'rad'
    phasedir=me.direction('J2000',ra,dec)
    t=me.epoch("UTC",qu.quantity(time))
    me.do_frame(t)

    p = me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
    me.do_frame(p)
    azel = me.measure(phasedir,'azel')
    return azel

def azel2radec(az,el,time, pos):
    '''convert azimuth elevation to ra dec'''
    me=measures()
    if type(az)!=str:
        az=str(az)+'rad'
    if type(el)!=str:
        el=str(el)+'rad'
    phasedir=me.direction('AZELGEO',az,el)
    t=me.epoch("UTC",qu.quantity(time))
    me.do_frame(t)

    p = me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
    me.do_frame(p)

    radec = me.measure(phasedir,'RADEC')
    return radec


def getuvw(ra,dec,time, pos1,pos2):
    me=measures()
    if type(ra)!=str:
        ra=str(ra)+'rad'
    if type(dec)!=str:
        dec=str(dec)+'rad'
    phasedir=me.direction('J2000',ra,dec)
    me.do_frame(phasedir)
    t=me.epoch("UTC",qu.quantity(time))
    me.do_frame(t)
    
    p = me.position('ITRF',str(pos1[0])+'m',str(pos1[1])+'m',str(pos1[2])+'m')
    bl = me.baseline('ITRF',str(pos2[0]-pos1[0])+'m',str(pos2[1]-pos1[1])+'m',str(pos2[2]-pos1[2])+'m')
    logging.info (bl)
    me.do_frame(p)
    logging.info (me.to_uvw(bl)['xyz'])
    #return uvw
     
def getIONEXtimerange(timerange,timestep):
    '''Get list with array of time per day
    IONEX files go per day, check if more than one file is  needed.
        
    Args:
        timerange (list) : start time, end  time in MJD seconds
        timestep (float) : time stein secionds
    Returns:
        List[np.array] : list of np.array of times in seconds per day
    '''
    times=[]
    oldtimerange=-100
    while timerange[0]< timerange[1] and timerange[0]>oldtimerange:
      oldtimerange=timerange[0]
      #logging.info (timerange)
      result =  obtain_observation_year_month_day_fraction(timerange[0])
      part_of_day = result[3]
      result2 =  obtain_observation_year_month_day_fraction(timerange[1])
      if result2[2]==result[2]:  #sameday
          #make sure to include the last timestep        
          times.append(np.arange(timerange[0],timerange[1]+timestep,timestep)) 
      else:
        nr_remaining_seconds=(1.-part_of_day) * 24.0 * 60. * 60.
        #logging.info ("new day",nr_remaining_seconds)
        times.append(np.arange(timerange[0],timerange[0]+nr_remaining_seconds,timestep))
      if len(times[-1]):
        timerange[0]=times[-1][-1]+timestep
    return times

def getlonlatheight(az,el,position,h=ION_HEIGHT):
    '''get piercepoint coordinates
    
    Args:
        az(float) : azimuth (radians)
        el(float) : elevation(radians)
        position (np.array) : xy ITRF station position in metres
        h (float) : ionospheric layer height in metres
    Returns:
        Tuple[float,float,float,float,float,float] : lat,lon,h piercepoint,
                                                     lon,lat station,
                                                     airmass factor
        
        '''
    if HAS_PYRAP:
      lonlat=getLonLatStation(az,el,pos=position)

      lon=lonlat['m0']['value']
      lat=lonlat['m1']['value']
      # convert to itrf coordinates on sphere with radius 1
      diritrf=[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]
      # calculate piercepoint in xyz(code from Bas vd Tol)

      (ppx1,ppy1,ppz1,am1) = getPP(h=h, mPosition=position,
                                   direction=diritrf)

      #get pp in lon,lat h
      me=measures()
      pp1position=me.position("ITRF", str(ppx1)+'m', str(ppy1)+'m', 
                              str(ppz1)+'m')

      lonpp = np.degrees(pp1position['m0']['value'])
      latpp = np.degrees(pp1position['m1']['value'])
      height = pp1position['m2']['value'] / 1.0e3
    elif HAS_EPHEM:
      slantRange = 3.0e20    # seem to need a large value to
                               # get near the WNB measures equivalent
                               #             lat,lon,ht = aer2ecef(np.degrees(body.az), np.degrees(body.alt), slantRange, np.degrees(geocentric_latitude), location_lon, ION_HEIGHT)
      location_lat, location_lon, location_height = ITRFToWGS84(position[0], position[1], position[2])
      lat,lon,ht = aer2ecef(np.degrees(az), np.degrees(el), slantRange, location_lat, location_lon, ION_HEIGHT)
      # convert to itrf coordinates on sphere with radius 1
      lat = np.radians(lat)
      lon = np.radians(lon)
      diritrf=[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]
      # calculate piercepoint in xyz(code from Bas vd Tol)

      (ppx1,ppy1,ppz1,am1) = getPP(h=ION_HEIGHT, mPosition=position,
                                   direction=diritrf)
      #get pp in lon,lat h
      latpp, lonpp, height =  ITRFToWGS84(ppx1,ppy1,ppz1)
    else:
      logging.info ('unable to compute position parameters - exiting!')
      return -1,-1,-1
    return latpp, lonpp, height,lon,lat,am1




def getAzEl(pointing,time,position,ha_limit=-1000):
    '''Get azimuth Elevation
    
    Args:
        pointing (list) : ra,dec  in radians
        time (float) : time in MJD seconds
        position (list) : ITRF xyz station coordinates
        ha_limit (float) : optional limit on hour angle
    Returns:
        Tuple[float,float] : azimuth elevation (radians)
        
    '''

    if HAS_PYRAP:
        if ha_limit==-1000:
            azel=radec2azel(pointing[0],pointing[1],time=str(time)+'s',pos=position)
            az=azel['m0']['value']
            el=azel['m1']['value']
        else:
            me=measures()
            p=me.position("ITRF",str(position[0])+'m',str(position[1])+'m',str(position[2])+'m')
            t=me.epoch("UTC",qu.quantity(str(time)+'s'))
            phasedir=me.direction('J2000',str(pointing[0])+'rad',str(pointing[1])+'rad')
            me.doframe(p)
            me.doframe(t)
            hadec=me.measure(phasedir,"HADEC")
            if abs(hadec['m0']['value'])>ha_limit:
                logging.info ("below horizon %s %f %f"%(
                        tab.taql('calc ctod($time s)')[0], 
                        np.degrees(hadec['m0']['value']), 
                        np.degrees(hadec['m1']['value'])))
                return 999,999
            else:
                azel=me.measure(phasedir,"AZELGEO")
  
                az=azel['m0']['value']
                el=azel['m1']['value']
    elif HAS_EPHEM:
        if ha_limit!=-1000:
            logging.info ("limiting on HA/DEC not implemented for PyEphem yet, ignoring")
        location_lat, location_lon, location_height = ITRFToWGS84(position[0], position[1], position[2])
        location = ephem.Observer()
        # convert geodetic latitude to geocentric
        # flattening, f, defined above for WGS84 stuff
        geodet_lat = np.radians(location_lat)
        #tan_geocentric_latitude =  np.tan(geodet_lat) * (1 - f) **2
        geocentric_latitude = GeodeticToGeocentricLat(geodet_lat,
                                                      location_height)
        location.lat = geocentric_latitude
        location.lon = np.radians(location_lon)
        location.elevation = location_height
        location.pressure = 0.0
        #  convert to Dublin Julian Date for PyEphem
        location.date =  time/86400.0 - 15019.5
        #lst = location.sidereal_time()
        equatorial = ephem.Equatorial(str(12.0 * np.degrees(pointing[0])/
                                          180),str(np.degrees(pointing[1])))
        body = ephem.FixedBody()
        body._ra = equatorial.ra
        body._dec = equatorial.dec
        body._epoch = equatorial.epoch
        body.compute(location)
        az = np.degrees(body.az)   * np.pi / 180.0
        el = np.degrees(body.alt) * np.pi / 180.0
    else: 
      logging.info ('failure to get azimuth and elevation! Exiting!')
      return -1,-1
    return az,el


def get_time_range(start_time,end_time,timestep,time_in_sec,TIME_OFFSET=0):
    '''get  time range in MJD seconds for measures like time specification'''
    if HAS_PYRAP:
      try:
        start_time = qu.quantity(start_time).get_value()
        end_time = qu.quantity(end_time).get_value()
        logging.info ('**** specified start and end time ', start_time, end_time)
        reference_time = start_time * 86400.0 - TIME_OFFSET
        st = reference_time - timestep 
        et = end_time * 86400.0 +  timestep + TIME_OFFSET
        timerange= [st, et]
      except:
        logging.info ('no time range given')
        logging.info ('exiting')
        return -1,-1,-1
    elif HAS_EPHEM:
      if time_in_sec:
        dublin_start = start_time / 86400.0 -15019.5
        dublin_end = end_time / 86400.0 -15019.5
        start_time = ephem.julian_date(ephem.Date(dublin_start)) - 2400000.5
        end_time = ephem.julian_date(ephem.Date(dublin_end)) - 2400000.5
      else:
        start_time = ephem.julian_date(ephem.Date(start_time)) - 2400000.5
        end_time = ephem.julian_date(ephem.Date(end_time)) - 2400000.5
      logging.info ('ephem start and end time ', start_time, end_time)
      reference_time = start_time * 86400.0 - TIME_OFFSET
      st = reference_time - timestep 
      et = end_time * 86400.0 + timestep + TIME_OFFSET
      timerange = [st, et]
    else:
      logging.info ('unable to get time range so exiting!')
      return -1,-1,-1
    return timerange,reference_time
