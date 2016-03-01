import scipy.io.netcdf as NC
import numpy as np
import os
from commons.time_interval import TimeInterval
from basins.region import Region, Rectangle

from instruments import bio_float
from instruments.var_conversions import FLOATVARS

TI = TimeInterval('20150520','20150830','%Y%m%d')
reg = Rectangle(-6,36,30,46)

modelvarname='P_i'

Profilelist=bio_float.FloatSelector(FLOATVARS[modelvarname],TI,reg)


#T2 = TimeInterval('20150501','20150601','%Y%m%d')
T2 = TimeInterval('20150501','20160215','%Y%m%d')
reg2=Rectangle(10,20,30,40)
Profilelist_2 = bio_float.FloatSelector('CHLA',T2,reg2)
nP = len(Profilelist_2)
Lon = np.zeros((nP,), np.float64)
Lat = np.zeros((nP,), np.float64)

for ip, p in enumerate(Profilelist_2):
    Lon[ip] = p.lon
    Lat[ip] = p.lat
    F=p._my_float
    print F.filename
    Pres,Prof,Profile_adj,Qc = F.read_very_raw('CHLA')


# unique delle stringhe wmo = p.name()
S=set()
for p in Profilelist_2: S.add( p.name() )
#SET TO LIST
Slist=list(S)


# Find the maximum and its depth [DCM]:
MAX1 = max(Prof)
DCM1 = Pres[np.argmax(Prof)]

# Plot on a map:
import pylab as pl
coastline=np.load('Coastline.npy')
c_lon=coastline['Lon']
c_lat=coastline['Lat']
pl.plot(c_lon,c_lat,'k')
pl.plot(Lon,Lat,'ro')
pl.show(block=False)

############TEST
lonrad = np.zeros((len(Profilelist_2),), np.float64)
latrad = np.zeros((len(Profilelist_2),), np.float64)
dist = np.zeros((len(Profilelist_2)-1,), np.float64)

#S2=set()
Cyclelist=[]
Timelist=[]
LonlistR=[]
LatlistR=[]
LonlistD=[]
LatlistD=[]

for ind in range(0, len(Profilelist_2)):
    wmolist =  Profilelist_2[ind]._my_float.wmo
#    if wmolist == Slist[2] :
    if wmolist == '6901768' :
#       for p in Profilelist_2: S2.add( Profilelist_2[ind]._my_float.cycle() );
#       S2.add( Profilelist_2[ind]._my_float.cycle() )
       Cyclelist.append(Profilelist_2[ind]._my_float.cycle) 
       print Profilelist_2[ind]._my_float.wmo
       print Profilelist_2[ind]._my_float.cycle
       print Profilelist_2[ind]._my_float.time
       print Profilelist_2[ind]._my_float.lon, Profilelist_2[ind]._my_float.lat
       lonrad[ind]=np.pi/180.*Profilelist_2[ind]._my_float.lon
       latrad[ind]=np.pi/180.*Profilelist_2[ind]._my_float.lat
#       dist[ind]=np.arccos(np.sin
       print lonrad[ind], latrad[ind]
       Timelist.append(Profilelist_2[ind]._my_float.time)
       LonlistR.append(lonrad[ind])
       LatlistR.append(latrad[ind])
       LonlistD.append(Profilelist_2[ind]._my_float.lon)
       LatlistD.append(Profilelist_2[ind]._my_float.lat)

np.asarray(Timelist)
np.asarray(LonlistR)
np.asarray(LatlistR)
np.asarray(LonlistD)
np.asarray(LatlistD)

# CALCULATE DISTANCES
npoints=len(LonlistR)
ndist=npoints-1

Lon1= LonlistR[0]
Lat1=LatlistR[0]
Lon2=LonlistR[1]
Lat2=LatlistR[1]
#for k in np.linspace(0,ndist,npoint):
dist=np.arccos(np.sin(Lat1)*np.sin(Lat2) + np.cos(Lat1)*np.cos(Lat2)* np.cos(np.abs(Lon2-Lon1)))*6371
#Out[61]: 5.0835963434982814
Dist= np.arccos(np.sin(Lat1)*np.sin(Lat2) + np.cos(Lat1)*np.cos(Lat2)* np.cos(np.abs(Lon2-Lon1)))*6371

# Radius of Earth: mean = 6371 Km; otherwise:
RE=6378.1370 # EQUATORIAL RADIUS
RP=6356.7523 # POLAR RADIUS

RR=np.sqrt(((RE*RE*np.cos(Lat1))**2 + (RP*RP*np.sin(Lat1))**2)/((RE*np.cos(Lat1))**2 + (RP*np.sin(Lat1))**2))
M=((RP*RE)**2 / ((RE*np.cos(Lat1))**2 + (RP*np.sin(Lat1))**2)**(3/2) )

Llat = LatlistR
Llon = LonlistR
#Llat = np.zeros(npoints, np.float64)
#Llon = np.zeros(npoints, np.float64)
Ddist = np.zeros(ndist, np.float64)
#Dtime = np.zeros(ndist, np.integer)
Dt = np.zeros(ndist, np.float64)

for k in range(0,ndist): 
#    Ddist[k] = np.arccos(np.sin(Llat[k])*np.sin(Llat[k+1]) + np.cos(Llat[k])*np.cos(Llat[k+1])* np.cos(np.abs(Llon[k+1]-Llon[k])))*6371
    Ddist[k] = np.arccos(np.sin(Llat[k])*np.sin(Llat[k+1]) + np.cos(Llat[k])*np.cos(Llat[k+1])* np.cos(np.abs(Llon[k+1]-Llon[k])))*RR

#    Dtime[k] = Timelist[k+1].day-Timelist[k].day
    Dtime = Timelist[k+1]-Timelist[k]
#    Dt[k]=(Dtime.days*86400 + Dtime.seconds)/86400.
    Dt[k]=(Dtime.days*86400 + Dtime.seconds)/3600.
out_file = open("test.txt","w")
out_file.write("distance in km: , Ddist[k], and hours between the two records: ,Dt[k] , Lon=,LonlistD[k], to ,LonlistD[k+1],; Lat=,LatlistD[k], to ,LatlistD[k+1], \n")
out_file.close()
#    print "distance in km: ", Ddist[k], "and hours between the two records: ",Dt[k] , "Lon=",LonlistD[k]," to ",LonlistD[k+1],"; Lat=",LatlistD[k]," to ",LatlistD[k+1]

print "distance in km: ", Ddist, "and days between the two records: ",Dt


for pp, p in enumerate(Profilelist_2): print pp, p._my_float.wmo , p._my_float.cycle

for pp, p in enumerate(Profilelist_2): print pp, p._my_float.wmo , p._my_float.cycle ,p._my_float.lon , p._my_float.lat , p._my_float.time
# 102 6901768 62 18.6155 37.6098 2016-01-14 10:34:00
F1=Profilelist[102]._my_float
Pres,Prof,Prof_adj,Qc=F1.read_very_raw('CHLA')
F1.plot(Pres,Prof)
F1.plot(Pres,Prof_adj)


