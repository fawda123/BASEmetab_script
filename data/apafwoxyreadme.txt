README 

2022-10-17

r000.csv

Synthetic dataset to be used in EBASE development 
Created by Maria Herrmann
~/magic_nsf/forward_oxy_model/foxym_v3/src/main_foxym using
inputs from Apalachicola CatPoint 2021 from Jill Arriola
---------------------------------------------------------------------------

Varables:

datet	date and time	
oxy	    dissolved oxygen, mmol/m3
troc	time-rate-of-change, mmol/m2/d
gasex	gase exchange (positive is up), mmol/m2/d
gpp	    Gross Primary Production, mmol/m2/d
er	    Ecosystem Respiration (positive is loss), mmol/m2/d
ht	    Height of the water column	m
aparam  light use efficiency parameter for gpp, (mmolO2/m2/d)/(W/m2)
oxysu	oxygen supersaturation (O2sat - O2obs), mmmol/m3
wspd2	Wind speed squared, m2/s2
sc	    Schmidt number for oxygen, dimensionless
kw	    Gas transfer velocity (Wanninkhof 2014), m/s
par	    PAR,W/m2
temp	temperature,degC
salt	salinity, ppt

Mass balance: troc = gpp - er - gasex