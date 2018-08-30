import RMextract.getRM as gt
import pyrap.tables as tab
from pylab import *

rc("lines",lw=3)
rc('font', family='serif', size=16)


#a=tab.taql('calc MJD("2013/03/02/17:02:54")')[0]*3600*24
#b=tab.taql('calc MJD("2013/03/03/01:02:50")')[0]*3600*24
a=tab.taql('calc MJD("2014/10/16/17:01:41")')[0]*3600*24
b=tab.taql('calc MJD("2014/10/17/05:35:33")')[0]*3600*24

statpos=gt.PosTools.posCS002
pointing=array([ 0.,  0.5*np.pi]) #NCP
pointing=array([ .15374123,  .8415521 ]) #3C196
#tec3=tec2=tec=0
#pointing=array([3.7146860578925645, 0.9111636804140731  ]) #3C295
#tec=gt.getRM(ionexPath='/home/mevius/IONEXdata/',earth_rot=0,ha_limit=1*np.pi,radec=pointing,timestep=300, timerange = [a, b],stat_positions=[statpos,],server="ftp://gnss.oma.be/gnss/products/IONEX/",prefix="ROBR")
#tec2=gt.getRM(ionexPath='/home/mevius/IONEXdata/',earth_rot=0,ha_limit=1*np.pi,radec=pointing,timestep=300, timerange = [a, b],stat_positions=[statpos,],prefix='iltg')
tec3=gt.getRM(ionexPath='/home/mevius/IONEXdata/',server="ftp://igs-final.man.olsztyn.pl/pub/gps_data/GPS_IONO/cmpcmb/",prefix="igsg",earth_rot=1,radec=pointing,timestep=150, timerange = [a, b+150],stat_positions=[statpos,])
#tec2=gt.getRM(ionexPath='/home/mevius/IONEXdata/',server='ftp://213.184.6.172/home/gnss/products/ionex/',earth_rot=0.5,ha_limit=1*np.pi,radec=pointing,timestep=150, timerange = [a, b+150],stat_positions=[statpos,],prefix='ILTF')
#tec=gt.getRM(ionexPath='/home/mevius/IONEXdata/',server='ftp://213.184.6.172/home/gnss/products/ionex/',earth_rot=0.5,ha_limit=1*np.pi,radec=pointing,timestep=150, timerange = [a, b+150],stat_positions=[statpos,],prefix='ILTQ')
#tec3=gt.getRM(ionexPath='/home/mevius/IONEXdata/',earth_rot=1,ha_limit=1*np.pi,radec=pointing,timestep=100, timerange = [a, b],stat_positions=[statpos,])
times=tec3['times']
flags=tec3['flags']['st1']
timess=[tm/(3600*24.) for tm in times]
dates=tab.taql('calc ctod($timess)')



if 'array' in dates.keys():
    dates=dates['array']
else:
    dates=dates[dates.keys()[0]] #backward compatibility with older casacore vresions

format="%Y/%m/%d/%H:%M:%S.%f"
mydatetimes=[datetime.datetime.strptime(mydate,format) for mydate in dates]
#maskeddata=np.ma.array(tec2['RM']['st1'],mask=np.logical_not(flags))
#plot_date(mydatetimes,maskeddata,'-')
#maskeddata=np.ma.array(tec['RM']['st1'],mask=np.logical_not(flags))
#plot_date(mydatetimes,maskeddata,'-')
maskeddata=np.ma.array(tec3['RM']['st1'],mask=np.logical_not(flags))
plot_date(mydatetimes,maskeddata,'-')
plt.gcf().autofmt_xdate()
ylabel("RM (rad/m^2)")
#title("RM variation direction 3C196 2015/07/22")
title("RM variation direction 3C196 2013/03/02")
show()
