import os,sys
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Region, Rectangle
from layer_integral import coastline
import instruments
from instruments import bio_float
from instruments.var_conversions import FLOATVARS
import scipy.io.netcdf as NC
import numpy as np
from commons.utils import addsep
import pylab as pl
from basins import OGS
from validation.online.profileplotter import figure_generator, ncwriter, add_metadata
import matplotlib
from mhelpers.pgmean import PLGaussianMean

TI1 = TimeInterval('20130101','20160301','%Y%m%d')
reg1 = [OGS.med]
reg_sn = ['med']

varname = ['CHLA','DOXY','NITRATE','TEMP','PSAL']
plotvarname = [r'Chl $[mg/m^3]$',r'Oxy $[mmol/m^3]$',r'Nitr $[mmol/m^3]$',r'Temp $[^\circ C]$','Sal']
read_adjusted = [True,False,False,False,False]
mapgraph = [3,4,5,1,2]

meanObj11 = PLGaussianMean(11,1.0)
Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)

# SELECT SINGLE FLOATS:
# IDENTIFICA IL FLOAT ID PER STAGIONE E SALVALO NELLA LISTA S1:
S1=set()
for k in Profilelist_1: S1.add( k.name() )
S1list=list(S1)
nS1=len(S1list)
snS1=str(nS1)
S1l=np.asarray(S1list)
nP = len(Profilelist_1)

# PICK ALL PROFILES FOR ONE SINGLE FLOAT
for ind, fl in enumerate(S1list) :
	wmolist = S1l[ind]
	print ind, " - WMOLIST: ", wmolist
#	fig = pl.figure(ind)

# CREATE THE LIST OF PROFILES FOR A SINGLE FLOAT
# SALVA NELLA LISTA "Float" i profili di un determinato FLOAT ID:
        Float = []
        for ip, pp in enumerate(Profilelist_1):
	   if wmolist == pp._my_float.wmo :	   
	   	print pp._my_float.wmo , pp._my_float.cycle
	   	Float.append(pp._my_float)
	print "FLOAT ",wmolist,": number of profiles: ",len(Float)
	All_Profiles = np.zeros((len(Float),201,len(varname)), np.float64)
###########
# FIGURE SETUP:

	fig, axs = pl.subplots(2,3, facecolor='w', edgecolor='k')
	fig_name = ''.join(['FLOAT_',wmolist,'.png'])
    	hsize=10
    	vsize=12
    	fig.set_size_inches(hsize,vsize)
    	fig.subplots_adjust(hspace = 0.15, wspace=0.3)
    	axs = axs.ravel()



###########
# SELECT THE VARIABLE:
	for ipp, ppp in enumerate(Float) : 

#########
# TRAJECTORY FIGURE:
    	    ax = axs[0]
    	    c_lon, c_lat=coastline.get()
    	    ax.plot(c_lon,c_lat, color='#000000',linewidth=0.5)
    	    ax.plot(ppp.lon,ppp.lat,'ro')
    	    ax.set_xticks(np.arange(-6,36,2))
    	    ax.set_yticks(np.arange(0,100,2))
     	    ax.set_xlabel("lon")
    	    ax.set_ylabel("lat")

	    extent=10 #degrees
	    ax.set_xlim([ppp.lon -extent/2, ppp.lon+extent/2])
	    ax.set_ylim([ppp.lat -extent/2, ppp.lat+extent/2])
	    bbox=ax.get_position()

	    deltax, _ =bbox.size
	    new_deltay = deltax* hsize/vsize
	    bottom = bbox.ymax - new_deltay
	    ax.set_position([bbox.xmin, bottom, deltax, new_deltay])
#	    floatlabel = 'Float '+ str(ppp.wmo) + ' \n' #+ season_text2
	    floatlabel = 'Float '+ str(ppp.wmo) + ' \n' + 'n.profiles: ' + str(len(Float))
	    ax.set_title(floatlabel)

#########
#########
	    for ax in axs[1:]:
              	 ax.set_ylim(0,200)
              	 ax.locator_params(axis='x',nbins=4)
              	 ax.yaxis.grid()

            for ax in [axs[2], axs[4], axs[5]]:
                 ax.set_yticklabels([])

#########################

  	    MEAN1 = np.zeros((201,len(varname)), np.float64)
            STD1 = np.zeros((201,len(varname)), np.float64)
	    DCM1 = np.zeros((201,len(varname)), np.float64)
            MAX1 = np.zeros((201,len(varname)), np.float64)
	    for i_var, var in enumerate(varname) :

		NewPres_1m=np.linspace(0,200,201)
		if var in ppp.available_params.split(" ") :		
			Pres,Prof,Qc = ppp.read(var,read_adjusted=read_adjusted[i_var])

# SELECT JUST DEPTHS BETWEEN 0-200M
# SELEZIONA SOLO LE PROFONDITA' FINO A 200M:
            		ii200 = Pres<200 ;
            		if len(Prof[ii200])>0 :
# INTERPOLATE THE CHLA DATA TO 1m VERTICAL RESOLUTION:
              			NewProf_1m = np.interp(NewPres_1m,Pres[ii200],Prof[ii200])
				Prof_smooth = meanObj11.compute(NewProf_1m,NewPres_1m)
				All_Profiles[ipp,:,i_var] = Prof_smooth

	for i_var in range(0,len(varname)):
		MEAN1[:,i_var]= np.mean(All_Profiles[:,:,i_var],0)
		STD1[:,i_var] = np.std(All_Profiles[:,:,i_var],0)

		ax=axs[mapgraph[i_var]] 
		ax.plot(MEAN1[:,i_var],NewPres_1m,'b')
		ax.plot(MEAN1[:,i_var]+STD1[:,i_var],NewPres_1m,'b:')
		ax.plot(MEAN1[:,i_var]-STD1[:,i_var],NewPres_1m,'b:')
		ax.invert_yaxis()
		title_str = ''.join([plotvarname[i_var]])
		ax.set_title(title_str)

	pl.show(block=False)
	pl.savefig(fig_name)
	pl.close(fig)

sys.exit()
#--------------------------------------

