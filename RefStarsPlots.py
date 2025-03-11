import numpy as np
import matplotlib.pyplot as plt
import astropy
import pandas
from astroquery.simbad import Simbad
import time
import astropy.units as u
from astropy.coordinates import SkyCoord
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch,Rectangle
from astropy.coordinates import HeliocentricTrueEcliptic
from astropy.coordinates import HeliocentricMeanEcliptic
from astropy.coordinates import Galactic
from astropy.coordinates import ICRS
from astropy.coordinates import GeocentricTrueEcliptic
from astropy.time import Time as astroTime
from astropy.table import vstack
simbad = Simbad()
simbad.add_votable_fields("pmra",
							"pmdec",
							"sp_type",
							"U",
							"B",
							"V",
							"G",
							"R",
							"I",
							"J",
							"H",
							"K",
							"plx_value",
							"rvz_radvel")

def queryTap_byname(name):
	#Return a TAP query of a single object.
	#Inputs:
	#	Name - Valid Simbad ID of an object
	#Returns:
	#	The astropy table of the query result
	example_base = """SELECT TOP 1 oid, main_id, ra, dec, coo_err_maj, coo_err_min, coo_err_angle, pmra, pmdec, plx_value, rvz_radvel, otypes, V, G
						FROM basic JOIN alltypes ON basic.oid = alltypes.oidref JOIN allfluxes ON basic.oid = allfluxes.oidref JOIN ident ON basic.oid = ident.oidref
						WHERE id = '{insertname}'
  						AND ra IS NOT NULL
  						AND dec IS NOT NULL;"""
	example = example_base.format(insertname=name)

	queryresult = simbad.query_tap(example)

	return queryresult

def queryTap_byVmagAsym(VmagLow,VmagHigh):
	#Return a TAP query that satisfies the criteria of being within a certain range of V magnitudes.
	#Inputs:
	#	VmagLow - Lower bound of possible target V mag
	#	VmagHigh - Upper bound of possible target V mag
	#Returns:
	#	The astropy table of the query result
	example_base = """SELECT TOP 50 oid, main_id, ra, dec, coo_err_maj, coo_err_min, coo_err_angle, pmra, pmdec, plx_value, rvz_radvel, otype, otypes, sp_type, U, B, V, G, R, I, J, H, K
						FROM basic JOIN allfluxes ON basic.oid = allfluxes.oidref JOIN alltypes ON basic.oid = alltypes.oidref
						WHERE (otype = '*') AND (otype != 'V*..') AND (otype != '**..') AND (V <= {VmagUp}) AND (V >= {VmagDown})
						AND (otypes) NOT LIKE ('%V*%')
						AND (otypes) NOT LIKE ('%**%')
  						AND ra IS NOT NULL
  						AND dec IS NOT NULL
  						AND plx_value IS NOT NULL
						AND pmra IS NOT NULL
						AND pmdec IS NOT NULL
						AND rvz_radvel IS NOT NULL
  						;"""
	example = example_base.format(VmagUp=str(VmagHigh),VmagDown=str(VmagLow))

	queryresult = simbad.query_tap(example)

	return queryresult

def queryTap_byVmagUp(VmagHigh):
	#Return a TAP query that satisfies the criteria of being within a certain range of V magnitudes.
	#Inputs:
	#	VmagLow - Lower bound of possible target V mag
	#	VmagHigh - Upper bound of possible target V mag
	#Returns:
	#	The astropy table of the query result
	#example_base = """SELECT TOP 50 oid, main_id, ra, dec, coo_err_maj, coo_err_min, coo_err_angle, pmra, pmdec, plx_value, rvz_radvel, otype, otypes, sp_type, U, B, V, G, R, I, J, H, K
	#					FROM basic JOIN allfluxes ON basic.oid = allfluxes.oidref JOIN alltypes ON basic.oid = alltypes.oidref
	#					WHERE (otype = '*') AND (V <= {VmagUp})
  	#					AND ra IS NOT NULL
  	#					AND dec IS NOT NULL
  	#					AND plx_value IS NOT NULL
	#					AND pmra IS NOT NULL
	#					AND pmdec IS NOT NULL
	#					AND rvz_radvel IS NOT NULL
  	#					;"""
	example_base = """SELECT TOP 200 oid, main_id, ra, dec, coo_err_maj, coo_err_min, coo_err_angle, pmra, pmdec, plx_value, rvz_radvel, otype, otypes, sp_type, V
						FROM basic JOIN allfluxes ON basic.oid = allfluxes.oidref JOIN alltypes ON basic.oid = alltypes.oidref
						WHERE (otypes) LIKE ('*|%') AND (V <= {VmagUp})
  						AND ra IS NOT NULL
  						AND dec IS NOT NULL
  						;"""
	example = example_base.format(VmagUp=str(VmagHigh))

	queryresult = simbad.query_tap(example)

	print(queryresult)

	return queryresult

def simplePlot1():
	mag855 = queryTap_byVmagUp(3.0)

	mag855Dist = 1000./mag855["plx_value"].value

	#print(mag855)

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	ics855 = SkyCoord(ra=mag855['ra'],dec=mag855['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag855Dist*u.pc,pm_ra_cosdec=mag855["pmra"],pm_dec=mag855["pmdec"],
		radial_velocity=mag855["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	newICS855 = ics855.apply_space_motion(new_obstime=epochtime)

	ra_rad855 = newICS855.ra.wrap_at(180 * u.deg).radian
	dec_rad855 = newICS855.dec.radian

	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'},dpi=1000)

	ax1.grid(True)

	ax1.plot(ra_rad855,dec_rad855,marker='o',linestyle='None',color='navy',markersize=10,label='V<3')

	#ax1.set_title('Core Throughput Standards, Heliocentric True Ecliptic, J2027, eq=2027',fontsize=18)

	ax1.set_xlabel('ICRS RA [deg]',fontsize=18)

	ax1.set_ylabel('ICRS Dec [deg]',fontsize=18)

	#ax1.legend(bbox_to_anchor=(1.2,0),loc='lower right',fontsize=16)

	fig1.tight_layout()

	fig1.savefig('ReferenceStarsV_v2.jpg')

	plt.show()

def simplePlot2():

	referencelist = pandas.read_csv('ReferenceList.csv')

	PSFstatus = referencelist["PSFstan?"].values

	names = referencelist["Name"].values

	polstatus = referencelist["Pol Measurement"].values

	fixnames = []

	for i in np.arange(len(names)):
		fixnames.append('* '+names[i])

	print(fixnames)

	#allstars = simbad.query_objects(fixnames[:42])
	allstars = simbad.query_objects(fixnames)

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	#print(allstars["rvz_radvel"])

	ics = SkyCoord(ra=allstars['ra'],dec=allstars['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=1000./allstars["plx_value"]*u.pc,pm_ra_cosdec=allstars["pmra"],pm_dec=allstars["pmdec"],
		radial_velocity=allstars["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	#eclT = ics.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))
	eclT = ics.transform_to(GeocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))
	#print(eclT)
	#print(ics.transform_to(ICRS(obstime=epochtime)))

	newICS = ics.apply_space_motion(new_obstime=epochtime)

	#neweclT = eclT.apply_space_motion(new_obstime=epochtime)

	#print(ics)

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	#newECL = newICS.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))

	newECL = newICS.transform_to(GeocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad = newECL.lon.wrap_at(180 * u.deg).radian
	dec_rad = newECL.lat.radian

	#ra_rad = newICS.ra.wrap_at(180 * u.deg).radian
	#dec_rad = newICS.dec.radian
	print(len(ra_rad))

	#legend_elements1 = [Patch(color='lime', label='Good'),
    #					Patch(color='gold', label='Good?'),
    #					Patch(color='indianred', label='Bad?')]#,
                   #Line2D([0], [0], marker='^', color='black', label='Pol Measured',
                   #       markerfacecolor='black', markersize=10),
                   #Line2D([0], [0], marker='o', color='black', label='Not Measured',
                   #       markerfacecolor='black', markersize=10)]
	legend_elements1 = [Line2D([0], [0], marker='s', color='w', label='Grade A',
                          markerfacecolor='lime', markersize=10,linestyle='None',markeredgecolor='black'),
    					Line2D([0], [0], marker='o', color='w', label='Grade B',
                          markerfacecolor='gold', markersize=10,linestyle='None',markeredgecolor='black'),
    					Line2D([0], [0], marker='^', color='w', label='Grade C',
                          markerfacecolor='indianred', markersize=10,linestyle='None',markeredgecolor='black')]

	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'})#,dpi=1000)

	ax1.grid(True)

	#for i in np.arange(len(fixnames[:42])):
	for i in np.arange(len(fixnames[:41])):
		if PSFstatus[i] == 'good':
			PSFcol = 'lime'
			PSFmark = 's'
			#PSFcol = 'navy'
		elif PSFstatus[i] == 'good?':
			PSFcol = 'gold'
			PSFmark = 'o'
			#PSFcol = 'navy'
		elif PSFstatus[i] == 'bad?':
			PSFcol = 'indianred'
			PSFmark = '^'
			#PSFcol = 'navy'


		#PSFmark = 'o'
		#PSFcol = 'navy'

		if PSFstatus[i] != 'bad':
			ax1.plot(ra_rad[i],dec_rad[i],marker=PSFmark,linestyle='None',color=PSFcol,markersize=10,markeredgecolor="black")
			if fixnames[i] == '* eta Cen' or fixnames[i] == '* alf Lep' or fixnames[i] == '* bet Eri' or fixnames[i] == '* del Leo' or fixnames[i] == '* eps UMa' or fixnames[i] == '* eps CMa':
				ax1.text(ra_rad[i]-np.deg2rad(3),dec_rad[i],fixnames[i],fontweight='extra bold',horizontalalignment='right',verticalalignment='top')
			else:
				ax1.text(ra_rad[i]+np.deg2rad(3),dec_rad[i],fixnames[i],fontweight='extra bold',horizontalalignment='left',verticalalignment='top')

	ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan')

	ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')

	#ax1.legend(handles=legend_elements1,loc=4,fontsize=16)#,bbox_to_anchor=(1.0, 0.0, 1., 0.1))
	ax1.legend(handles=legend_elements1,fontsize=16,bbox_to_anchor=(1.0, 1.0),bbox_transform=fig1.transFigure)

	#ax1.set_title('Current Reference Stars',fontsize=18)

	ax1.set_xlabel('l [deg]',fontsize=18)
	ax1.set_ylabel('b [deg]',fontsize=18)

	#ax1.set_xlabel('ICRS RA [deg]',fontsize=18)

	#ax1.set_ylabel('ICRS Dec [deg]',fontsize=18)

	fig1.tight_layout()

	#fig1.savefig('ReferenceStarsVD_v2.jpg')

	fig1.savefig('RefStar_Prop_Named.jpg')

	plt.show()

def GapID(Ecl_RA):
	diff_arr = np.array([]) #array of gap edges
	for i in np.arange(len(Ecl_RA)):
#		if i == 0:
#			bottom_end = Ecl_RA[i]-5
		if i>0:
			difference = Ecl_RA[i]-Ecl_RA[i-1]
			if difference > 10*u.deg:
				centralDIFFedgeR = Ecl_RA[i]-(5*u.deg)
				centralDIFFedgeL = Ecl_RA[i-1]+(5*u.deg)
				diff_arr = np.append(diff_arr,centralDIFFedgeL.value)
				diff_arr = np.append(diff_arr,centralDIFFedgeR.value)
		if Ecl_RA[i] == Ecl_RA[-1]:
			difference = (Ecl_RA[0]+(180*u.deg))-Ecl_RA[i]
			if difference > 10*u.deg:
				centralDIFFedgeR = Ecl_RA[0]+(180-5)*u.deg
				centralDIFFedgeL = Ecl_RA[i]+(5*u.deg)
				diff_arr = np.append(diff_arr,centralDIFFedgeL.value)
				diff_arr = np.append(diff_arr,centralDIFFedgeR.value)

	return diff_arr

def skyCoverageplots():
	referencelist = pandas.read_csv('ReferenceList.csv')

	PSFstatus = referencelist["PSFstan?"].values

	names = referencelist["Name"].values

	polstatus = referencelist["Pol Measurement"].values

	fixnames = []

	for i in np.arange(len(names)):
		fixnames.append('* '+names[i])

	GoodsInds = []
	for i in np.arange(len(fixnames[:42])):
		if PSFstatus[i] == 'good' or PSFstatus[i] == 'good?':
			GoodsInds.append(i)

	#print(GoodsInds)

	GoodInds = []
	for i in np.arange(len(fixnames[:42])):
		if PSFstatus[i] == 'good':
			GoodInds.append(i)

	allstars = simbad.query_objects(fixnames[:42])
	#allstars = simbad.query_objects(fixnames)

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	#print(allstars["rvz_radvel"])

	ics = SkyCoord(ra=allstars['ra'],dec=allstars['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=1000./allstars["plx_value"]*u.pc,pm_ra_cosdec=allstars["pmra"],pm_dec=allstars["pmdec"],
		radial_velocity=allstars["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	eclT = ics.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))
	#print(eclT)
	#print(ics.transform_to(ICRS(obstime=epochtime)))

	newICS = ics.apply_space_motion(new_obstime=epochtime)

	#neweclT = eclT.apply_space_motion(new_obstime=epochtime)

	#print(ics)

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	newECL = newICS.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	lon_rad = newECL.lon.wrap_at(180 * u.deg).radian
	lat_rad = newECL.lat.radian

	lon_all = newECL.lon.wrap_at(180 * u.deg)

	lon_goods = newECL.lon.wrap_at(180*u.deg)[GoodsInds]

	lon_good = newECL.lon.wrap_at(180*u.deg)[GoodInds]

	Ecl_lon_ALLFOLD = np.copy(newECL.lon.degree)*u.deg
	Ecl_lon_GYFOLD = np.copy(newECL.lon.degree[GoodsInds])*u.deg
	Ecl_lon_GFOLD = np.copy(newECL.lon.degree[GoodInds])*u.deg
	print(Ecl_lon_ALLFOLD)

	for i in np.arange(len(Ecl_lon_ALLFOLD)):
		if Ecl_lon_ALLFOLD[i]>180*u.deg:
			Ecl_lon_ALLFOLD[i] = Ecl_lon_ALLFOLD[i]-180*u.deg
	for i in np.arange(len(Ecl_lon_GYFOLD)):
		if Ecl_lon_GYFOLD[i]>180*u.deg:
			Ecl_lon_GYFOLD[i] = Ecl_lon_GYFOLD[i]-180*u.deg
	for i in np.arange(len(Ecl_lon_GFOLD)):
		if Ecl_lon_GFOLD[i]>180*u.deg:
			Ecl_lon_GFOLD[i] = Ecl_lon_GFOLD[i]-180*u.deg

	Ecl_ALL = np.sort(Ecl_lon_ALLFOLD)
	Ecl_GY = np.sort(Ecl_lon_GYFOLD)
	Ecl_G = np.sort(Ecl_lon_GFOLD)

	#print('All ',Ecl_ALL)
	#print('GY ',Ecl_GY)
	#print('G ',Ecl_G)

	#print(lon_goods)

	#print('Longitudes',lon_all)

	#AllGaps = GapID(np.sort(lon_all))

	#GoodsGaps = GapID(np.sort(lon_goods))

	#GoodGaps = GapID(np.sort(lon_good))

	AllGaps = GapID(Ecl_ALL)

	GoodsGaps = GapID(Ecl_GY)

	GoodGaps = GapID(Ecl_G)

	print('Gaps',AllGaps)

	legend_elements1 = [Patch(color='lime', label='Good'),
    					Patch(color='gold', label='Good?'),
    					Patch(color='indianred', label='Bad?')]#,
                   #Line2D([0], [0], marker='^', color='black', label='Pol Measured',
                   #       markerfacecolor='black', markersize=10),
                   #Line2D([0], [0], marker='o', color='black', label='Not Measured',
                   #       markerfacecolor='black', markersize=10)]

	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'})

	ax1.grid(True)

	for i in np.arange(len(fixnames[:42])):
		if PSFstatus[i] == 'good':
			PSFcol = 'lime'
		elif PSFstatus[i] == 'good?':
			PSFcol = 'gold'
		elif PSFstatus[i] == 'bad?':
			PSFcol = 'indianred'

		PSFmark = 'o'

		if PSFstatus[i] != 'bad':
			ax1.plot(lon_rad[i],lat_rad[i],marker=PSFmark,linestyle='None',color=PSFcol,markersize=10)

	#ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan')

	#ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')

	ax1.legend(handles=legend_elements1,loc=4,fontsize=16)

	#ax1.set_title('Current Reference Stars',fontsize=18)

	ax1.set_xlabel('l [deg]',fontsize=18)
	ax1.set_ylabel('b [deg]',fontsize=18)

	height = 180

	



	for i in np.arange(0,len(AllGaps),2):
		leftboundary = AllGaps[i]
		rightboundary = AllGaps[i+1]
		width = rightboundary-leftboundary

		ax1.axvspan(np.deg2rad(leftboundary),np.deg2rad(rightboundary),alpha=0.5,color='indianred')
		ax1.axvspan(np.deg2rad(leftboundary-180), np.deg2rad(rightboundary-180), alpha=0.5, color='indianred')

		if leftboundary>180 and rightboundary>180:
			ax1.axvspan(np.deg2rad(leftboundary-360), np.deg2rad(rightboundary-360), alpha=0.5, color='indianred')
		elif leftboundary<180 and rightboundary>180:
			ax1.axvspan(np.deg2rad(leftboundary-360), np.deg2rad(rightboundary-360), alpha=0.5, color='indianred')
		
		#region = Rectangle((leftboundary,-90),width,height,color='indianred',alpha=0.5)
		#ax1.add_patch(region)
		#if GoodGaps[i]>180:

		#	leftboundary = GoodGaps[i]-180

		#	rightboundary = GoodGaps[i+1]-180
		#	width = rightboundary-leftboundary

		#	ax1.axvspan(np.deg2rad(rightboundary),np.deg2rad(leftboundary),alpha=0.5,color='indianred')

			#region = Rectangle((leftboundary,-90),width,height,color='indianred',alpha=0.5)
			#ax1.add_patch(region)


	fig1.tight_layout()

	fig1.savefig('AllStarsGood.jpg')

	plt.show()

#Make the proposed multi-instrument detection plot

	#ra_rad = newICS.ra.wrap_at(180 * u.deg).radian
	#dec_rad = newICS.dec.radian





