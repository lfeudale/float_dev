import scipy.io.netcdf as NC
import numpy as np
import os
from commons.time_interval import TimeInterval
from commons.layer import Layer
import pylab as pl

from profiler import *
import basins.OGS as OGS
from instruments import bio_float
from instruments.var_conversions import FLOATVARS
from commons.layer import Layer

# QUI SETTI IL RANGE IN CUI VUOI CREARE I DATI DI MISFIT: AL MOMENTO SONO SOLO PER L'ANNO 2015
# DA ESTENDERE DAL 20130101 NONAPPENA SONO DISPONIBILI I DATI DI ANALISI DEL MODELLO OPERATIVO 
TI1 = TimeInterval('20150101','20160101','%Y%m%d')

M = Matchup_Manager(TI1,INPUTDIR,BASEDIR)

maskfile    = os.getenv("MASKFILE"); 
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

layer=Layer(0,200)
modelvarname='P_i'

# GENERA UN VETTORE PER LA PROFONDITA' DA 0 A 200M CON 1M DI RISOLUZIONE, SU CUI
# INTERPOLARE I DATI:
NewPres=np.linspace(0,200,201)

# ESTRAI I PROFILI DI CLOROFILLA CAMPIONATI NELL'INTERVALLO DI TEMPO TI1
# PER L'INTERA AREA MEDITERRANEA
Profilelist_1=bio_float.FloatSelector(FLOATVARS[modelvarname],TI1,OGS.med)

# INIZIALIZZA LE QUATTRO STAGIONI CREANDO DELLE LISTE VUOTE DOVE VERRANNO 
# ARCHIVIATI I PROFILI CHE CADONO IN OGNUNA DI QUESTE STAGIONI
Winter_list = []
Spring_list = []
Summer_list = []
Autumn_list = []

for w in Profilelist_1 :
	if (w.time.month >=1) & (w.time.month <=3) :
		Winter_list.append(w)
	if (w.time.month >=4) & (w.time.month <=6) :
                Spring_list.append(w)
        if (w.time.month >=7) & (w.time.month <=9) :
                Summer_list.append(w)
        if (w.time.month >=10) & (w.time.month <=12) :
                Autumn_list.append(w)

# INIZIA IL LOOP SULLE STAGIONI
for LOOP in range(1,5) :
    if LOOP == 1:
       season_list = Winter_list
       print 'WINTER'
       season_text = "wn"
       season_text2 = "[JAN-MAR]"
    if LOOP == 2:
       season_list = Spring_list
       print 'SPRING' 
       season_text = "sp"
       season_text2 = "[APR-JUN]"
    if LOOP == 3:
       season_list = Summer_list
       season_text = "sm"
       print 'SUMMER'
       season_text2 = "[JUL-SEP]"
    if LOOP == 4:
       season_list = Autumn_list
       season_text = "at"
       print 'AUTUMN'
       season_text2 = "[OCT-DEC]"

# IDENTIFICA IL FLOAT ID PER STAGIONE E SALVALO NELLA LISTA S1:
    S1=set()
    for k in season_list: S1.add( k.name() )
    S1list=list(S1)
    nS1=len(S1list)
    snS1=str(nS1)
    S1l=np.asarray(S1list)
    nP = len(season_list)

# PER SINGOLO FLOAT E PER STAGIONE, FARE L'ANALISI SUI PROFILI
    for ind, fl in enumerate(S1list) :
#      if ind == 0 :
	wmolist = S1l[ind]
	print ind, " - WMOLIST: ", wmolist
	fig = pl.figure(ind)

# SALVA NELLA LISTA "Float" i profili di un determinato FLOAT ID 
# PER UNA DELLE 4 STAGIONI: 
        Float = []
        for i, p in enumerate(season_list):
	  if wmolist == p._my_float.wmo :
	   
	   	print p._my_float.wmo , p._my_float.cycle
	   	Float.append(p)
	
	AllProfiles = np.zeros((len(Float),201),np.float64)
# ESTRAI I PROFILI DI CLOROFILLA DEL BIOFLOAT E QUELLI ANALOGHI DEL MODELLO:
 	for ip, pp in enumerate(Float) 
	    singlefloatmatchup = M.getMatchups([pp], nav_lev, modelvarname)
	    s200 = singlefloatmatchup.subset(layer)
	    print s200.bias()
	    if np.invert(np.isnan(s200.bias())) :
	    #print 'DIFF:'
	    #print s200.diff()
# EFFETTUA LI'INTERPOLAZIONE A 1M DEI VALORI DI DIFFERENZIALE DI CLOROFILLA
	    	s200int=np.interp(NewPres,s200.Depth,s200.diff())

# SALVA TUTTI I MATCHUP DEI PROFILI (PER FLOAT ID E STAGIONE) IN UN UNICA MATRICE TEMPORANEA
	    	AllProfiles[ip,:] = s200int

# EFFETTUA CALCOLI DI VALORE MEDIO E STD DEL MODELLO-OBS
	pmean=(np.mean(AllProfiles,0))
	pstd=(np.std(AllProfiles,0))
	pl.plot(pmean,NewPres,'b')
	pl.plot(pmean+pstd,NewPres,'b:')
	pl.plot(pmean-pstd,NewPres,'b:')
	pl.gca().invert_yaxis()
	floatlabel = 'Float '+ str(wmolist) + ' ' + season_text2 + ' \n' + 'n.profiles: ' + str(len(Float))
	fig_name = ''.join(['DIFF_Obs_Mod_',wmolist,'_',season_text,'.png'])
	fig.suptitle(floatlabel)
	pl.xlabel(r'Chl $[mg/m^3]$')
	pl.ylabel(r'Depth $[m]$')
	pl.savefig(fig_name)
	pl.show(block=False)
	pl.close(fig)

