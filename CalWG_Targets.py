import numpy as np
import matplotlib.pyplot as plt
import astropy
import pandas
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
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
from astropy.time import Time as astroTime
from astropy.table import vstack
from astroquery.xmatch import XMatch
from astropy.table import Table, QTable, Column, MaskedColumn
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

def absFluxCoverage(catsave=False,plotsave=False):
	#This function queries Simbad for the properties of faint and bright absolute flux standards. It creates a sky coverage map after applying sky motion and catalogs the properties of the objects in a csv format.
	#Inputs:
	#	catsave: Save the queried target properties to a csv file, Default: False
	#	plotsave: Save the created sky coverage map as a jpg, Default: False
	#
	#
	#Returns:
	#	None
	#Absolute flux standards are pre-defined choices.
	dimStandards = ['2MASS J18083474+6927286',
					'2MASS J18120957+6329423',
					'2MASS J17571324+6703409',
					'2MASS J18052927+6427520',
					'2MASS J14515797+7143173',
					'2MASS J17551622+6610116',
					'2MASS J18022716+6043356',
					'2MASS J17325264+7104431',
					'2MASS J16313382+3008465',
					'2MASS J17403468+6527148']

	truncDim = ['J1808+6927',
					'J1812+6329',
					'J1757+6703',
					'J1805+6427',
					'J1451+7143',
					'J1755+6610',
					'J1802+6043',
					'J1732+7104',
					'J1631+3008',
					'J1740+6527']

	brightStandards = ['* 109 Vir',
						'* alf Lyr',
						'* eta UMa',
						'* ksi02 Cet']

	truncBright = ['109 Vir',
						'Vega',
						'eta UMa',
						'ksi02 Cet']

	#Use astroquery's query_objects to query Simbad for dim standards

	dimstars = simbad.query_objects(dimStandards)

	#Add a column to the astropy data table to specify these are dim standards

	dimstars['CalType'] = ['DimAbsFlux' for i in np.arange(len(dimStandards))]

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	#Establish SkyCoord objects to apply space motion

	dimics = SkyCoord(ra=dimstars['ra'],dec=dimstars['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=1000./dimstars["plx_value"]*u.pc,pm_ra_cosdec=dimstars["pmra"],pm_dec=dimstars["pmdec"],
		radial_velocity=dimstars["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))


	dimeclT = dimics.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))

	newdimICS = dimics.apply_space_motion(new_obstime=epochtime)

	#Transform to Heliocentric Truc Ecliptic

	newdimECL = newdimICS.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	dimra_rad = newdimECL.lon.wrap_at(180 * u.deg).radian
	dimdec_rad = newdimECL.lat.radian

	#Add in a wait so that Simbad doesn't blacklist your IP for 30 minutes

	time.sleep(1.0)

	#Query Simbad for Bright standards

	brightstars = simbad.query_objects(brightStandards)

	brightstars['CalType'] = ['BriAbsFlux' for i in np.arange(len(brightStandards))]

	#Combine the dim and bright standards list into one table

	combinedStandards = vstack([dimstars,brightstars])

	#Save the catalog to a csv if desired

	if catsave == True:
		#Remove extraneous column
		combinedStandards.remove_column('object_number_id')
		df = combinedStandards.to_pandas()
		df.to_csv('AbsoluteFluxStandards.csv')

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	brightics = SkyCoord(ra=brightstars['ra'],dec=brightstars['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=1000./brightstars["plx_value"]*u.pc,pm_ra_cosdec=brightstars["pmra"],pm_dec=brightstars["pmdec"],
		radial_velocity=brightstars["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	brighteclT = brightics.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))

	newbrightICS = brightics.apply_space_motion(new_obstime=epochtime)

	newbrightECL = newbrightICS.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	brightra_rad = newbrightECL.lon.wrap_at(180 * u.deg).radian
	brightdec_rad = newbrightECL.lat.radian

	#Make sky coverage plot

	fig1,ax1 = plt.subplots(1,1,figsize=(10,6),subplot_kw={'projection': 'mollweide'})

	ax1.grid(True)

	#Plot CVZs


	ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan')

	ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')

	ax1.plot(dimra_rad,dimdec_rad,marker='o',linestyle='None',color='indianred',markersize=10,label='Faint Standard')

	ax1.plot(brightra_rad,brightdec_rad,marker='*',linestyle='None',color='goldenrod',markersize=10,label='Bright Standard')

	#Label each standard's name

	for i in np.arange(len(truncDim)):
		ax1.text(dimra_rad[i],dimdec_rad[i],truncDim[i],fontweight='extra bold',horizontalalignment='center',verticalalignment='top')
	for i in np.arange(len(truncBright)):
		ax1.text(brightra_rad[i],brightdec_rad[i],truncBright[i],fontweight='extra bold',horizontalalignment='center',verticalalignment='top')

	ax1.legend(loc='best')

	ax1.set_xlabel('l [deg]',fontsize=18)
	ax1.set_ylabel('b [deg]',fontsize=18)

	ax1.set_title('Absolute Flux Standards Epoch J2027, Equinox 2027',fontsize=18)

	fig1.tight_layout()

	if plotsave == True:

		fig1.savefig('AbsFluxStandards.jpg')

	plt.show()

def queryTap_byVmag(Vmag,cone=0.1):
	#Return a TAP query that satisfies the criteria of being within a certain range of a V magnitude.
	#Inputs:
	#	Vmag - Target V magnitude you want to query
	#	cone - Range around Vmag that is acceptable for the query. Default = 0.1
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
	example = example_base.format(VmagUp=str(Vmag+cone),VmagDown=str(Vmag-cone))

	queryresult = simbad.query_tap(example)

	return queryresult

def queryTap_byname_old(name):
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

def queryTap_byname(name):
	example_base = """SELECT TOP 1 oid, main_id, ids, ra, dec, coo_err_maj, coo_err_min, coo_err_angle, coo_qual, coo_bibcode,
						pmra, pmdec, pm_err_maj, pm_err_min, pm_err_angle, pm_qual, pm_bibcode, plx_value, plx_err, plx_qual, plx_bibcode,
						rvz_radvel, rvz_err, rvz_bibcode,
						otype, otypes, sp_type, sp_bibcode, 
						U, B, V, G, R, I, J, H, K
						FROM basic JOIN alltypes ON basic.oid = alltypes.oidref JOIN allfluxes ON basic.oid = allfluxes.oidref JOIN flux ON basic.oid = flux.oidref
						JOIN ident ON basic.oid = ident.oidref JOIN ids ON basic.oid = ids.oidref
						WHERE main_id = '{insertname}'
  						AND ra IS NOT NULL
  						AND dec IS NOT NULL;"""
	#example_base = """SELECT TOP 1 basic.oid, basic.main_id, ids.ids, basic.ra, basic.dec, basic.coo_err_maj, basic.coo_err_min, basic.coo_err_angle, basic.coo_qual, basic.coo_bibcode,
	#					basic.pmra, basic.pmdec, basic.pm_err_maj, basic.pm_err_min, basic.pm_err_angle, basic.pm_qual, basic.pm_bibcode, basic.plx_value, basic.plx_err, basic.plx_qual, basic.plx_bibcode,
	#					basic.rvz_radvel, basic.rvz_err, basic.rvz_bibcode,
	#					basic.otype, alltypes.otypes, basic.sp_type, basic.sp_bibcode, 
	#					allfluxes.U, allfluxes.B, allfluxes.V, allfluxes.G, allfluxes.R, allfluxes.I, allfluxes.J, allfluxes.H, allfluxes.K,
	#					mesVar.vartyp,mesVar.lowVmax,mesVar.vmax,mesVar.r_vmax,mesVar.magtyp,mesVar.uppVmin,mesVar.upperiod,mesVar.period,mesVar.r_period,mesVar.epoch,mesVar.r_epoch
	#					FROM basic JOIN alltypes ON basic.oid = alltypes.oidref JOIN allfluxes ON basic.oid = allfluxes.oidref JOIN flux ON basic.oid = flux.oidref
	#					JOIN ident ON basic.oid = ident.oidref JOIN ids ON basic.oid = ids.oidref JOIN mesVar ON basic.oid = mesVar.oidref
	#					WHERE basic.main_id = '{insertname}'
  	#					AND ra IS NOT NULL
  	#					AND dec IS NOT NULL;"""
	example = example_base.format(insertname=name)

	queryresult = simbad.query_tap(example)

	time.sleep(1)

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

def coreThruput(catsave=False,plotsave=False):
	#This function queries Simbad for the properties of core throughput calibration targets. It creates a sky coverage map after applying sky motion and catalogs the properties of the objects in a csv format.
	#Inputs:
	#	catsave: Save the queried target properties to a csv file, Default: False
	#	plotsave: Save the created sky coverage map as a jpg, Default: False
	#
	#
	#Returns:
	#	None

	#Roman-Coronagraph core throughput standards cannot use the ND filter. Therefore, they must be fainter than V=10.9 or else they will saturate. An upper bound of 12 is set to restrict the range of possible targets.
	mag855 = queryTap_byVmagAsym(10.7,12.0)

	#Fill in missing values necessary for Sky motion. The assumption here is that missing values are essentially minimal.

	#mag855["plx_value"] = mag855["plx_value"].filled(0.00001)
	#mag855["pmra"] = mag855["pmra"].filled(0.0)
	#mag855["pmdec"] = mag855["pmdec"].filled(0.0)
	#mag855["rvz_radvel"] = mag855["rvz_radvel"].filled(0.0)

	mag855Dist = 1000./mag855["plx_value"].value

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	mag855['CalType'] = ['CoreThruput' for i in np.arange(50)]

	#Save to csv if desired

	if catsave == True:
		df = mag855.to_pandas()
		df.to_csv('CoreThruputTargs.csv')


	#Establish SkyCoord objects

	ics855 = SkyCoord(ra=mag855['ra'],dec=mag855['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag855Dist*u.pc,pm_ra_cosdec=mag855["pmra"],pm_dec=mag855["pmdec"],
		radial_velocity=mag855["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	eclT855 = ics855.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))

	newICS855 = ics855.apply_space_motion(new_obstime=epochtime)

	#Transform to Heliocentric True Ecliptic

	newECL855 = newICS855.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad855 = newECL855.lon.wrap_at(180 * u.deg).radian
	dec_rad855 = newECL855.lat.radian

	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'})

	ax1.grid(True)

	#Plot CVZs

	ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan')

	ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')

	ax1.plot(ra_rad855,dec_rad855,marker='o',linestyle='None',color='indianred',markersize=10,label='V>10.7')

	ax1.set_title('Core Throughput Standards, Heliocentric True Ecliptic, J2027, eq=2027',fontsize=18)

	ax1.set_xlabel('l [deg]',fontsize=18)

	ax1.set_ylabel('b [deg]',fontsize=18)

	ax1.legend(bbox_to_anchor=(1.2,0),loc='lower right',fontsize=16)

	fig1.tight_layout()

	if plotsave == True:

		fig1.savefig('CoreThruputStandards.jpg')

	plt.show()

def commissioningFlats(catsave=False,plotsave=False):
	#This function queries Simbad for the properties of commissioning flat calibration targets. It creates a sky coverage map after applying sky motion and catalogs the properties of the objects in a csv format.
	#Inputs:
	#	catsave: Save the queried target properties to a csv file, Default: False
	#	plotsave: Save the created sky coverage map as a jpg, Default: False
	#
	#
	#Returns:
	#	None
	mag855 = queryTap_byVmag(8.26,cone=0.05)

	#mag855["plx_value"] = mag855["plx_value"].filled(0.00001)
	#mag855["pmra"] = mag855["pmra"].filled(0.0)
	#mag855["pmdec"] = mag855["pmdec"].filled(0.0)
	#mag855["rvz_radvel"] = mag855["rvz_radvel"].filled(0.0)

	mag855Dist = 1000./mag855["plx_value"].value

	mag855['CalType'] = ['CommFlat' for i in np.arange(50)]

	#Save to csv if desired

	if catsave == True:
		df = mag855.to_pandas()
		df.to_csv('CommissioningFlats.csv')

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	ics855 = SkyCoord(ra=mag855['ra'],dec=mag855['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag855Dist*u.pc,pm_ra_cosdec=mag855["pmra"],pm_dec=mag855["pmdec"],
		radial_velocity=mag855["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	eclT855 = ics855.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))

	newICS855 = ics855.apply_space_motion(new_obstime=epochtime)

	#Transform to Heliocentric True Ecliptic

	newECL855 = newICS855.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad855 = newECL855.lon.wrap_at(180 * u.deg).radian
	dec_rad855 = newECL855.lat.radian

	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'})

	ax1.grid(True)

	ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan')

	ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')

	ax1.plot(ra_rad855,dec_rad855,marker='o',linestyle='None',color='indianred',markersize=10,label='V~8.26')

	ax1.set_title('Commissioning Flats, Heliocentric True Ecliptic, J2027, eq=2027',fontsize=18)

	ax1.set_xlabel('l [deg]',fontsize=18)

	ax1.set_ylabel('b [deg]',fontsize=18)

	ax1.legend(bbox_to_anchor=(1.2,0),loc='lower right',fontsize=16)

	fig1.tight_layout()

	if plotsave == True:

		fig1.savefig('CommissioningFlatStandards.jpg')

	plt.show()

def phaseRetrieval(catsave=False,plotsave=False):
	#This function queries Simbad for the properties of commissioning flat calibration targets. It creates a sky coverage map after applying sky motion and catalogs the properties of the objects in a csv format.
	#Inputs:
	#	catsave: Save the queried target properties to a csv file, Default: False
	#	plotsave: Save the created sky coverage map as a jpg, Default: False
	#
	#
	#Returns:
	#	None
	mag855 = queryTap_byVmag(5.0,cone=0.1)

	#mag855["plx_value"] = mag855["plx_value"].filled(0.00001)
	#mag855["pmra"] = mag855["pmra"].filled(0.0)
	#mag855["pmdec"] = mag855["pmdec"].filled(0.0)
	#mag855["rvz_radvel"] = mag855["rvz_radvel"].filled(0.0)

	mag855Dist = 1000./mag855["plx_value"].value

	mag855['CalType'] = ['PhaseRet' for i in np.arange(50)]

	#Save to csv if desired

	if catsave == True:
		df = mag855.to_pandas()
		df.to_csv('PhaseRetrieval.csv')

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	ics855 = SkyCoord(ra=mag855['ra'],dec=mag855['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag855Dist*u.pc,pm_ra_cosdec=mag855["pmra"],pm_dec=mag855["pmdec"],
		radial_velocity=mag855["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	eclT855 = ics855.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))

	newICS855 = ics855.apply_space_motion(new_obstime=epochtime)

	#Transform to Heliocentric True Ecliptic

	newECL855 = newICS855.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad855 = newECL855.lon.wrap_at(180 * u.deg).radian
	dec_rad855 = newECL855.lat.radian

	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'})

	ax1.grid(True)

	ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan')

	ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')

	ax1.plot(ra_rad855,dec_rad855,marker='o',linestyle='None',color='indianred',markersize=10,label='V~5')

	ax1.set_title('Phase Retrieval Stars, Heliocentric True Ecliptic, J2027, eq=2027',fontsize=18)

	ax1.set_xlabel('l [deg]',fontsize=18)

	ax1.set_ylabel('b [deg]',fontsize=18)

	ax1.legend(bbox_to_anchor=(1.2,0),loc='lower right',fontsize=16)

	fig1.tight_layout()

	if plotsave == True:

		fig1.savefig('PhaseRetrieval.jpg')

	plt.show()

def imageCorrections(catsave=False,plotsave=False):
	#This function queries Simbad for the properties of image correction calibration targets. It creates a sky coverage map after applying sky motion and catalogs the properties of the objects in a csv format.
	#Inputs:
	#	catsave: Save the queried target properties to a csv file, Default: False
	#	plotsave: Save the created sky coverage map as a jpg, Default: False
	#
	#
	#Returns:
	#	None

	magList = np.array([4.1,7.5,9.7,11.8]) #Array of V magnitudes that we want to query
	markers_arr = ['o','^','*','x'] #Plot markers for the sky map
	markers_col = ['indianred','goldenrod','navy','magenta'] #Marker colors for plotted points
	#Setting up the plot
	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'}) #This creates a Mollweide projection map
	ax1.grid(True) #Grid lines
	ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan') #Plot out CVZs as a shaded region
	ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')
	ax1.set_title('Image Corrections, Heliocentric True Ecliptic, J2027, eq=2027',fontsize=18) #Set a plot title
	ax1.set_xlabel('l [deg]',fontsize=18) #Set axes labels
	ax1.set_ylabel('b [deg]',fontsize=18)
	epochtime = astroTime('J2027',format='jyear_str') #Define epoch and equinox times for the ecliptic transformation
	equitime = astroTime(2027.0,format='decimalyear')
	for i in np.arange(len(magList)):
		if i != 0:
			mag = queryTap_byVmag(magList[i],cone=0.1)
		elif i == 0:
			mag = queryTap_byVmagAsym(magList[i]-0.1,magList[i]+0.5)
		print(mag)
    	#The following lines fill any missing values in your table query with the numbers in parentheses
    	#mag["plx_value"] = mag["plx_value"].filled(0.00001)
    	#mag["pmra"] = mag["pmra"].filled(0.0)
    	#mag["pmdec"] = mag["pmdec"].filled(0.0)
    	#mag["rvz_radvel"] = mag["rvz_radvel"].filled(0.0)
		magDist = 1000./mag["plx_value"].value #Calculate the system distance in pc
		ics = SkyCoord(ra=mag['ra'],dec=mag['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=magDist*u.pc,pm_ra_cosdec=mag["pmra"],pm_dec=mag["pmdec"],
		radial_velocity=mag["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str')) #Set up a SkyCoord object. Requires the values from Simbad Query.

		mag['CalType'] = ['ImagCorr_'+str(magList[i]) for j in np.arange(len(mag))]

		if i == 0:
			fullTable = mag
		else:
			fullTable = vstack([fullTable,mag])

		newICS = ics.apply_space_motion(new_obstime=epochtime) #Apply proper and radial motion to the new observing time of 2027

		newECL = newICS.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str'))) #Transform your proper motion applied ICRS coordinates to ecliptic coordinates

		ra_rad = newECL.lon.wrap_at(180 * u.deg).radian #Express your ecliptic coordinates (longitude and latitude) in radians for plotting, allow it to wrap around the projection to the other side
		dec_rad = newECL.lat.radian

    	#Plot your result
		ax1.plot(ra_rad,dec_rad,marker=markers_arr[i],linestyle='None',color=markers_col[i],markersize=10,label='V~'+str(magList[i]))

		time.sleep(5) #I put this step in because if you make too many simbad queries in a given time interval, your IP gets blocked for a period of time. So every time you complete a query, wait 5 seconds before starting the next one.
		print('Completed Mag '+str(magList[i]))

	#Save catalog if desired
	if catsave == True:
		df = fullTable.to_pandas()
		df.to_csv('ImageCorrections.csv')

	#Create a plot legend
	ax1.legend(bbox_to_anchor=(1.2,0),loc='lower right',fontsize=16)
	fig1.tight_layout() #Clean it up

	#Save plot if desired
	if plotsave == True:
		fig1.savefig('ImageCorrectStandards.jpg')

	plt.show()


def imageCorrections_old():
	#4 sets of targets. V = 5.5, 8.5, 11.5, and 13.5 +/- 0.5 mag. Ideally they are close together.

	magList = np.array([5.5,8.5,11.5,13.5])

	mag55 = queryTap_byVmag(5.5)

	mag55["plx_value"] = mag55["plx_value"].filled(0.00001)
	mag55["pmra"] = mag55["pmra"].filled(0.0)
	mag55["pmdec"] = mag55["pmdec"].filled(0.0)
	mag55["rvz_radvel"] = mag55["rvz_radvel"].filled(0.0)

	mag55Dist = 1000./mag55["plx_value"].value
	print(len(mag55Dist))

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	#print(allstars["rvz_radvel"])

	ics55 = SkyCoord(ra=mag55['ra'],dec=mag55['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag55Dist*u.pc,pm_ra_cosdec=mag55["pmra"],pm_dec=mag55["pmdec"],
		radial_velocity=mag55["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	eclT55 = ics55.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))
	#print(eclT)
	#print(ics.transform_to(ICRS(obstime=epochtime)))

	newICS55 = ics55.apply_space_motion(new_obstime=epochtime)

	#neweclT = eclT.apply_space_motion(new_obstime=epochtime)

	#print(ics)

	#print(ics.transform_to(HeliocentricMeanEcliptic(equinox=equitime,obstime=epochtime)))

	newECL55 = newICS55.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad55 = newECL55.lon.wrap_at(180 * u.deg).radian
	dec_rad55 = newECL55.lat.radian

	time.sleep(5)

	mag85 = queryTap_byVmag(8.5)

	mag85["plx_value"] = mag85["plx_value"].filled(0.00001)
	mag85["pmra"] = mag85["pmra"].filled(0.0)
	mag85["pmdec"] = mag85["pmdec"].filled(0.0)
	mag85["rvz_radvel"] = mag85["rvz_radvel"].filled(0.0)

	mag85Dist = 1000./mag85["plx_value"].value
	print(len(mag85Dist))

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	#print(allstars["rvz_radvel"])

	ics85 = SkyCoord(ra=mag85['ra'],dec=mag85['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag85Dist*u.pc,pm_ra_cosdec=mag85["pmra"],pm_dec=mag85["pmdec"],
		radial_velocity=mag85["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	eclT85 = ics85.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))
	#print(eclT)
	#print(ics.transform_to(ICRS(obstime=epochtime)))

	newICS85 = ics85.apply_space_motion(new_obstime=epochtime)

	#neweclT = eclT.apply_space_motion(new_obstime=epochtime)

	#print(ics)

	#print(ics.transform_to(HeliocentricMeanEcliptic(equinox=equitime,obstime=epochtime)))

	newECL85 = newICS85.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad85 = newECL85.lon.wrap_at(180 * u.deg).radian
	dec_rad85 = newECL85.lat.radian

	time.sleep(5)

	mag115 = queryTap_byVmag(11.5)

	mag115["plx_value"] = mag115["plx_value"].filled(0.00001)
	mag115["pmra"] = mag115["pmra"].filled(0.0)
	mag115["pmdec"] = mag115["pmdec"].filled(0.0)
	mag115["rvz_radvel"] = mag115["rvz_radvel"].filled(0.0)

	mag115Dist = 1000./mag115["plx_value"].value
	print(len(mag115Dist))

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	#print(allstars["rvz_radvel"])

	ics115 = SkyCoord(ra=mag115['ra'],dec=mag115['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag115Dist*u.pc,pm_ra_cosdec=mag115["pmra"],pm_dec=mag115["pmdec"],
		radial_velocity=mag115["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	eclT115 = ics115.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))
	#print(eclT)
	#print(ics.transform_to(ICRS(obstime=epochtime)))

	newICS115 = ics115.apply_space_motion(new_obstime=epochtime)

	#neweclT = eclT.apply_space_motion(new_obstime=epochtime)

	#print(ics)

	#print(ics.transform_to(HeliocentricMeanEcliptic(equinox=equitime,obstime=epochtime)))

	newECL115 = newICS115.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad115 = newECL115.lon.wrap_at(180 * u.deg).radian
	dec_rad115 = newECL115.lat.radian

	time.sleep(5)

	mag135 = queryTap_byVmag(13.5)

	mag135["plx_value"] = mag135["plx_value"].filled(0.00001)
	mag135["pmra"] = mag135["pmra"].filled(0.0)
	mag135["pmdec"] = mag135["pmdec"].filled(0.0)
	mag135["rvz_radvel"] = mag135["rvz_radvel"].filled(0.0)

	mag135Dist = 1000./mag135["plx_value"].value
	print(len(mag135Dist))

	epochtime = astroTime('J2027',format='jyear_str')
	equitime = astroTime(2027.0,format='decimalyear')

	#print(allstars["rvz_radvel"])

	ics135 = SkyCoord(ra=mag135['ra'],dec=mag135['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=mag135Dist*u.pc,pm_ra_cosdec=mag135["pmra"],pm_dec=mag135["pmdec"],
		radial_velocity=mag135["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str'))

	#print(ics.transform_to(HeliocentricTrueEcliptic(equinox=equitime,obstime=epochtime)))

	eclT135 = ics135.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2000.0,format='decimalyear')))
	#print(eclT)
	#print(ics.transform_to(ICRS(obstime=epochtime)))

	newICS135 = ics135.apply_space_motion(new_obstime=epochtime)

	#neweclT = eclT.apply_space_motion(new_obstime=epochtime)

	#print(ics)

	#print(ics.transform_to(HeliocentricMeanEcliptic(equinox=equitime,obstime=epochtime)))

	newECL135 = newICS135.transform_to(HeliocentricTrueEcliptic(equinox=astroTime(2027.0,format='decimalyear'),obstime=astroTime('J2027',format='jyear_str')))


	ra_rad135 = newECL135.lon.wrap_at(180 * u.deg).radian
	dec_rad135 = newECL135.lat.radian

	fig1,ax1 = plt.subplots(1,1,figsize=(12,6),subplot_kw={'projection': 'mollweide'})

	ax1.grid(True)

	ax1.axhspan(np.deg2rad(54), np.deg2rad(90), alpha=0.5, color='cyan')

	ax1.axhspan(np.deg2rad(-90), np.deg2rad(-54), alpha=0.5, color='cyan')

	ax1.plot(ra_rad55,dec_rad55,marker='o',linestyle='None',color='indianred',markersize=10,label='V~5.5')

	ax1.plot(ra_rad85,dec_rad85,marker='^',linestyle='None',color='goldenrod',markersize=10,label='V~8.5')

	ax1.plot(ra_rad115,dec_rad115,marker='*',linestyle='None',color='navy',markersize=10,label='V~11.5')

	ax1.plot(ra_rad135,dec_rad135,marker='x',linestyle='None',color='magenta',markersize=10,label='V~13.5')

	ax1.set_title('Image Corrections, Heliocentric Mean Ecliptic, J2027, eq=2027',fontsize=18)

	ax1.set_xlabel('l [deg]',fontsize=18)

	ax1.set_ylabel('b [deg]',fontsize=18)

	ax1.legend(bbox_to_anchor=(1.2,0),loc='lower right',fontsize=16)

	fig1.tight_layout()

	fig1.savefig('ImageCorrectStandards.jpg')

	plt.show()

def GaiaXMatchFilter():
	#Note: angDist keyword is the separation
	input_targetlist = 'CoreThruputTargs.csv'
	input_table = Table.read(input_targetlist)
	table_xmatch = XMatch.query(cat1=input_table,cat2='vizier:I/355/gaiadr3',max_distance=1*u.arcsec,colRA1='ra',colDec1='dec')
	print(table_xmatch[np.where(table_xmatch['angDist']<=0.5)])

def LoadExistingDatabase(filename='AbsoluteFluxStandards.csv'):
	loadedCSV = pandas.read_csv(filename)

	return loadedCSV

def CatalogBuild(filebase):
	prebuiltTable = LoadExistingDatabase(filename=filebase+'.csv')
	#prebuiltTable = LoadExistingDatabase(filename='UnpolStandards_Verified_IDs.csv')

	print(prebuiltTable)

	names = prebuiltTable["main_id"].values

	Gaia3_Keys = ["ra","ra_error","dec","dec_error","pmra","pmra_error",
					"pmdec","pmdec_error","parallax","parallax_error",
					"radial_velocity","radial_velocity_error","phot_g_mean_mag","ruwe"]
	Gaia2_Keys = ["ra","ra_error","dec","dec_error","pmra","pmra_error",
					"pmdec","pmdec_error","parallax","parallax_error",
					"radial_velocity","radial_velocity_error","phot_g_mean_mag"]
	Simbad_Gaia3Keys = ["ra_2016","ra_2016_err","dec_2016","dec_2016_err","pmra_2016","pmra_2016_err",
						"pmdec_2016","pmdec_2016_err","plx_value_2016","plx_value_err_2016",
						"rvz_radvel_2016","rvz_radvel_err_2016","Gmag_2016","RUWE_2016"]
	Simbad_Gaia2Keys = ["ra_2015.5","ra_2015.5_err","dec_2015.5","dec_2015.5_err","pmra_2015.5","pmra_2015.5_err",
						"pmdec_2015.5","pmdec_2015.5_err","plx_value_2015.5","plx_value_err_2015.5",
						"rvz_radvel_2015.5","rvz_radvel_err_2015.5","Gmag_2015.5"]

	DiameterKeys = ['LDD','e_LDD','LDDCHI','UDDB','UDDV','UDDR','UDDI','UDDJ','UDDH','UDDK']

	JMMCKeys = ['UDdiam','LDdiam','e_UDdiam','e_LDdiam','UDband','Band','UDDmeas_ref','BibCode']

	MasterJMMCKeys = ['UDDmeas','LDDmeas','e_UDDmeas','e_LDDmeas','UDDmeas_Band','LDDmeas_Band','UDDmeas_ref','LDDmeas_ref']

	Kervella22Keys = ['PlxH2','e_PlxH2','PlxG3','e_PlxG3','RUWE','snrPMaHG1','BinHG1','pmRAH2G2','e_pmRAH2G2','pmDEH2G2','e_pmDEH2G2','PMaRAH2G2',
					'e_PMaRAH2G2','pmDEH2G2','e_pmDEH2G2','PMaRAH2G2','e_PMaRAH2G2','PMaDEH2G2','e_PMaDEH2G2','snrPMaH2G2','BinH2G2','pmRAH2EG3a',
					'e_PMaRAH2EG3a','PMaDEH2EG3a','e_PMaDEH2EG3a','snrPMaH2EG3a','BinH2EG3a','pmRAH2EG3b','e_pmRAH2EG3b','pmDEH2EG3b','e_pmDEH2EG3b',
					'PMaRAH2EG3b','e_PMaRAH2EG3b','PMaDEH2EG3b','e_PMaDEH2EG3b','snrPMaH2EG3b','BinH2EG3b','dVt','e_dVt','dVtPA','e_dVtPA','M1','e_M1',
					'r_M1']
	MasterKervella22Keys = ['PlxH2','e_PlxH2','PlxG3','e_PlxG3','RUWE_EDR3','snrPMaHG1','BinHG1','pmRAH2G2','e_pmRAH2G2','pmDEH2G2','e_pmDEH2G2','PMaRAH2G2',
					'e_PMaRAH2G2','pmDEH2G2','e_pmDEH2G2','PMaRAH2G2','e_PMaRAH2G2','PMaDEH2G2','e_PMaDEH2G2','snrPMaH2G2','BinH2G2','pmRAH2EG3a',
					'e_PMaRAH2EG3a','PMaDEH2EG3a','e_PMaDEH2EG3a','snrPMaH2EG3a','BinH2EG3a','pmRAH2EG3b','e_pmRAH2EG3b','pmDEH2EG3b','e_pmDEH2EG3b',
					'PMaRAH2EG3b','e_PMaRAH2EG3b','PMaDEH2EG3b','e_PMaDEH2EG3b','snrPMaH2EG3b','BinH2EG3b','dVt_K22','e_dVt_K22','dVtPA_K22','e_dVtPA_K22','M1','e_M1',
					'r_M1']

	Kervella19Keys = ['RUWE','dVt','e_dVt','dVtPA','e_dVtPA']

	MasterKervella19Keys = ['RUWE_DR2','dVt_K19','e_dVt_K19','dVtPA_K19','e_dVtPA_K19']

	for i in np.arange(len(names)):
		if names[i][:4] != 'Gaia':
			Simbadresult = queryTap_byname(names[i])

			print(Simbadresult)

		#Simbadresult["PSF_Status"] = PSFstatus[i]

			Simbadresult["HD_id"] = HD_id(Simbadresult["ids"][0],startid='HD')

			Simbadresult["HIP_id"] = HD_id(Simbadresult["ids"][0],startid='HIP')

			Simbadresult["TIC_id"] = HD_id(Simbadresult["ids"][0],startid='TIC')

			Simbadresult["HR_id"] = HD_id(Simbadresult["ids"][0],startid='HR')

			Simbadresult["WDS_id"] = HD_id(Simbadresult["ids"][0],startid='WDS')

			Simbadresult["SBC9_id"] = HD_id(Simbadresult["ids"][0],startid='SBC9')

			Simbadresult = Table(Simbadresult,masked=True,copy=False)

			#check for masks

			dr2name = HD_id(Simbadresult["ids"][0],startid='Gaia DR2')

			if dr2name != str(np.ma.masked):

			#Simbadresult["Gaia_DR2"] = HD_id(Simbadresult["ids"][0],startid='Gaia DR2')
			#Simbadresult["Gaia_DR2"] = MaskedColumn([dr2name],mask=[False])
				Simbadresult["Gaia_DR2"] = str(dr2name)
				Simbadresult["Gaia_DR2"].mask = [False]
			#print(Simbadresult)
			else:
				Simbadresult["Gaia_DR2"] = dr2name
				Simbadresult["Gaia_DR2"].mask = [True]

			#print(Simbadresult)

		#Simbadresult["Gaia_EDR3"] = HD_id(Simbadresult["ids"][0],startid='Gaia EDR3')

			dr3name = HD_id(Simbadresult["ids"][0],startid='Gaia DR3')

			if dr3name != str(np.ma.masked):

			#Simbadresult["Gaia_DR2"] = HD_id(Simbadresult["ids"][0],startid='Gaia DR2')
			#Simbadresult["Gaia_DR3"] = MaskedColumn([dr3name],mask=[False])
				Simbadresult["Gaia_DR3"] = dr3name
				Simbadresult["Gaia_DR3"].mask = [False]
			#print(Simbadresult)
			else:
				Simbadresult["Gaia_DR3"] = dr3name
				Simbadresult["Gaia_DR3"].mask = [True]
			#print(Simbadresult)

		#Simbadresult["Gaia_DR3"] = HD_id(Simbadresult["ids"][0],startid='Gaia DR3')

		#Grab information from Gaia

		#if Simbadresult["Gaia_DR3"] != 'None':
			if Simbadresult["Gaia_DR3"].mask == False:
				print('Gaia DR3 '+Simbadresult["Gaia_DR3"][0])
				GaiaResult = queryGaiaDR3_byID('Gaia DR3 '+Simbadresult["Gaia_DR3"][0])

			#print(GaiaResult)

				for key in np.arange(len(Gaia3_Keys)):
					if GaiaResult[Gaia3_Keys[key]].mask == False:
						Simbadresult[Simbad_Gaia3Keys[key]] = GaiaResult[Gaia3_Keys[key]]
		#if Simbadresult["Gaia_DR2"] != 'None':
			if Simbadresult["Gaia_DR2"].mask == False:
				print('Gaia DR2 '+Simbadresult["Gaia_DR2"][0])
				GaiaResult = queryGaiaDR2_byID('Gaia DR2 '+Simbadresult["Gaia_DR2"][0])

			#print(GaiaResult)

				for key in np.arange(len(Gaia2_Keys)):
					if GaiaResult[Gaia2_Keys[key]].mask == False:
						Simbadresult[Simbad_Gaia2Keys[key]] = GaiaResult[Gaia2_Keys[key]]

		#Grab stellar diameters from JSDCv2

		#print(names[i])
			hdexample = 'HD'+Simbadresult["HD_id"][0]
			try:
				jsdc = query_JSDC('* '+names[i],hdexample.replace(" ",""))
				time.sleep(1)
			except:
				time.sleep(1)
			#jsdc = Table(names=('Name','Vmag','LDD','e_LDD','LDDCHI','UDDB','UDDV','UDDR','UDDI','UDDJ','UDDH','UDDK'),
			#			dtype=(str,float,float,float,float,float,float,float,float,float,float,float))
				NameMask = MaskedColumn(['not found'],name='Name',mask=[True])
				jsdc = Table([NameMask])

			if jsdc['Name'].mask == False:
				for diam in np.arange(len(DiameterKeys)):
					if jsdc[DiameterKeys[diam]].mask == False:
						Simbadresult[DiameterKeys[diam]] = jsdc[DiameterKeys[diam]]


		#Grab JMMC stellar diameters

			try:
				hdexample = 'HD'+Simbadresult["HD_id"][0]
			#print(hdexample.replace(" ",""))
				jmmc = query_JMMC(hdexample.replace(" ",""))
				time.sleep(1)
				print('JMMC Query Successful')
			except:
				time.sleep(1)
			#jsdc = Table(names=('Name','Vmag','LDD','e_LDD','LDDCHI','UDDB','UDDV','UDDR','UDDI','UDDJ','UDDH','UDDK'),
			#			dtype=(str,float,float,float,float,float,float,float,float,float,float,float))
				NameMask = MaskedColumn(['not found'],name='ID1',mask=[True])
				jmmc = Table([NameMask])

			if jmmc['ID1'].mask == False:
				print('Adding in JMMC')
				jmmc = Table(jmmc,masked=True,copy=False)
				for meas in np.arange(len(JMMCKeys)):
					if jmmc[JMMCKeys[meas]].mask == False:
						Simbadresult[MasterJMMCKeys[meas]] = jmmc[JMMCKeys[meas]]

		#Grab astrometry table info from Kervella+22

			hipexample = Simbadresult["HIP_id"][0]
			try:
				k22 = query_K22(hipexample)
				print('K22 Query Successful')
			except:
				time.sleep(1)
				NameMaskK22 = MaskedColumn(['not found'],name='HIP',mask=[True])
				k22 = Table([NameMaskK22])

			if k22['HIP'].mask == False:
				print('Adding in K22')
				k22 = Table(k22,masked=True,copy=False)
				for K22key in np.arange(len(Kervella22Keys)):
					if k22[Kervella22Keys[K22key]].mask == False:
						Simbadresult[MasterKervella22Keys[K22key]] = k22[Kervella22Keys[K22key]]

			try:
				k19 = query_K19(hipexample)
				print('K19 Query Successful')
			except:
				time.sleep(1)
				NameMaskK19 = MaskedColumn(['not found'],name='HIP',mask=[True])
				k19 = Table([NameMaskK19])

			if k19['HIP'].mask == False:
				print('Adding in K19')
				k19 = Table(k19,masked=True,copy=False)
				for K19key in np.arange(len(Kervella19Keys)):
					if k19[Kervella19Keys[K19key]].mask == False:
						Simbadresult[MasterKervella19Keys[K19key]] = k19[Kervella19Keys[K19key]]

		#print(Simbadresult["UDDmeas"])
		#print(Simbadresult["LDDmeas"])
		elif names[i][:4] == 'Gaia':
			#If you only have a Gaia DR3 ID query Gaia DR3 only, no stellar diameters or Kervella catalogs
			GaiaResult = queryGaiaDR3_byID(names[i])
			Simbadresult = GaiaResult
			#Simbadresult["main_id"] = names[i]

			for key in np.arange(len(Gaia3_Keys)):
					if GaiaResult[Gaia3_Keys[key]].mask == False:
						Simbadresult[Simbad_Gaia3Keys[key]] = GaiaResult[Gaia3_Keys[key]]





		#Add in the polarization columns to "Simbadresult"
		#Pol standards columns

		#PreBuiltColumnNames = ["P% (U)","P% error (U)","PA (U)","PA error (U)",
		#						"P% (B)","P% error (B)","PA (B)","PA error (B)",
		#						"P% (V)","P% error (V)","PA (V)","PA error (V)",
		#						"P% (R)","P% error (R)","PA (R)","PA error (R)",
		#						"P% (I)","P% error (I)","PA (I)","PA error (I)",
		#						"P% (J)","P% error (J)","PA (J)","PA error (J)",
		#						"P% (H)","P% error (H)","PA (H)","PA error (H)",
		#						"P% (K)","P% error (K)","PA (K)","PA error (K)",
		#						"P% (Unknown)","P% error (Unknown)","PA (Unknown)","PA error (Unknown)",
		#						"Catalog Name","Catalog Reference","Heiles IDCAT",
		#						"Heiles RA","Heiles Vmag","Heiles SpTy",
		#						"Blinov RA","Blinov Dec","Blinov G Mag","Low SNR Flag"
		#						]

		#Unpol standards columns
		#PreBuiltColumnNames = ["P% (R)","P% error (R)","PA (R)","PA error (R)",
		#						"P% (Unknown)","P% error (Unknown)","PA (Unknown)","PA error (Unknown)",
		#						"Catalog Name","Catalog Reference","Heiles IDCAT",
		#						"Heiles RA","Heiles Vmag","Heiles SpTy",
		#						"Blinov RA","Blinov Dec","Blinov G Mag","Low SNR Flag"
		#						]
		PreBuiltColumnNames = ["CalType"]

		for columnN in np.arange(len(PreBuiltColumnNames)):
			Simbadresult = addPrebuiltColumn2New(prebuiltTable,Simbadresult,i,columnName=PreBuiltColumnNames[columnN])





		if i == 0:
			fullTable = Simbadresult

		else:
			fullTable = vstack([fullTable,Simbadresult])

	print(fullTable)

	df = fullTable.to_pandas()

	df.to_csv(filebase+'_v2.csv')

def query_K22(hipname):
	hipname.replace(" ", "")

	Kervella22Keys = ['HIP','PlxH2','e_PlxH2','PlxG3','e_PlxG3','RUWE','snrPMaHG1','BinHG1','pmRAH2G2','e_pmRAH2G2','pmDEH2G2','e_pmDEH2G2','PMaRAH2G2',
					'e_PMaRAH2G2','pmDEH2G2','e_pmDEH2G2','PMaRAH2G2','e_PMaRAH2G2','PMaDEH2G2','e_PMaDEH2G2','snrPMaH2G2','BinH2G2','pmRAH2EG3a',
					'e_PMaRAH2EG3a','PMaDEH2EG3a','e_PMaDEH2EG3a','snrPMaH2EG3a','BinH2EG3a','pmRAH2EG3b','e_pmRAH2EG3b','pmDEH2EG3b','e_pmDEH2EG3b',
					'PMaRAH2EG3b','e_PMaRAH2EG3b','PMaDEH2EG3b','e_PMaDEH2EG3b','snrPMaH2EG3b','BinH2EG3b','dVt','e_dVt','dVtPA','e_dVtPA','M1','e_M1',
					'r_M1']

	#try:
	k22 = Vizier(catalog="J/A+A/657/A7/tablea1",columns=Kervella22Keys).query_constraints(HIP=hipname)[0]
	#print(k22)
	time.sleep(1)
	return k22

def query_K19(hipname):
	hipname.replace(" ", "")

	Kervella19Keys = ['HIP','RUWE','dVt','e_dVt','dVtPA','e_dVtPA']

	#try:
	k19 = Vizier(catalog="J/A+A/623/A72/hipgpma",columns=Kervella19Keys).query_constraints(HIP=hipname)[0]
	print(k19)
	time.sleep(1)
	return k19


def query_JSDC(starname,hdname):
	column_names = ['Name','Vmag','LDD','e_LDD','LDDCHI','UDDB','UDDV','UDDR','UDDI','UDDJ','UDDH','UDDK']
	hdname.replace(" ", "")

	try:

		jsdc = Vizier(catalog="II/346/jsdc_v2",columns=column_names).query_constraints(Name=starname)[0]
		time.sleep(1)
	except:
		jsdc = Vizier(catalog="II/346/jsdc_v2",columns=column_names).query_constraints(Name=hdname)[0]
		time.sleep(1)
	#except:
	#	jsdc = Table(names=('Name','Vmag','LDD','e_LDD','LDDCHI','UDDB','UDDV','UDDR','UDDI','UDDJ','UDDH','UDDK'),
	#					dtype=(str,float,float,float,float,float,float,float,float,float,float,float))
	

	return jsdc

def query_JMMC(hdname):
	column_names = ['ID1','ID2','UDdiam','LDdiam','e_LDdiam','Band','BibCode']

	#hdname.replace(" ", "")

	#try:

	#print(hdname)

	jmmc = Vizier(catalog="II/345/jmdc",columns=column_names).query_constraints(ID1=hdname)[0]
	time.sleep(1)

	#print(jmmc)

	jmmc = Table(jmmc,masked=True,copy=False)

	#Reverse chronological order

	jmmc_flip = jmmc[::-1]

	jmmc_flip = Table(jmmc_flip,masked=True,copy=False)

	counter = 0

	for row in np.arange(len(jmmc_flip)):
		singlerow = Table(jmmc_flip[row])
		#print(singlerow["UDdiam"][0])
		if singlerow["UDdiam"][0] > 0 and singlerow["LDdiam"][0] > 0:
			print('Both UD and LD meas available')
			#jmmc_flip[row]["e_UDdiam"] = jmmc_flip[row]["e_LDdiam"]
			#jmmc_flip[row]["UDband"] = jmmc_flip[row]["Band"]
			#jmmc_flip[row]["UDDmeas_ref"] = jmmc_flip[row]["BibCode"]

			singlerow.add_column(singlerow["e_LDdiam"][0],name="e_UDdiam")
			singlerow.add_column(singlerow["Band"][0],name="UDband")
			singlerow.add_column(singlerow["BibCode"][0],name="UDDmeas_ref")

			#print(singlerow)


			return singlerow
		else:
			counter+=1
			#print(counter)
	#print('Final ',counter)
	#print(len(jmmc_flip))
	if counter == len(jmmc_flip):
		for row in np.arange(len(jmmc_flip)):
			singlerow = Table(jmmc_flip[row])
			#print(singlerow)
			#print(singlerow["UDdiam"][0])
			#print(singlerow["LDdiam"][0])
			if singlerow["UDdiam"].mask == True and singlerow["LDdiam"][0] > 0:
				print('LDmeas Only')
				singlerow.add_column(np.ma.masked,name="e_UDdiam")
				singlerow.add_column(np.ma.masked,name="UDband")
				singlerow.add_column(np.ma.masked,name="UDDmeas_ref")

				return singlerow
			elif singlerow["UDdiam"][0] > 0 and singlerow["LDdiam"].mask == True:
				singlerow.add_column(singlerow["e_LDdiam"][0],name="e_UDdiam")
				singlerow.add_column(singlerow["Band"][0],name="UDband")
				singlerow.add_column(singlerow["BibCode"][0],name="UDDmeas_ref")
				print('UDmeas Only')
				#print(singlerow)
				
				return singlerow



def HD_id(idstring,startid='HD'):
	HDstart = idstring.find(startid+' ')
	#print(HDstart)

	if HDstart != -1:
		extracted = idstring[HDstart+len(startid)+1:]
		endpoint = extracted.find('|')
		if endpoint != -1:
			final = extracted[:endpoint]
		else:
			final = extracted
	elif HDstart == -1:
		final = str(np.ma.masked)
		#print('Masked',final)

	return final

def queryGaiaDR3_byID(designation):
	example_base = """SELECT TOP 1 gaia_source.designation,gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,
						gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,
						gaia_source.ruwe,gaia_source.phot_g_mean_mag,
						gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_variable_flag,gaia_source.non_single_star,
						gaia_source.has_xp_continuous,gaia_source.has_xp_sampled,gaia_source.has_rvs,gaia_source.has_epoch_photometry,
						gaia_source.has_epoch_rv,gaia_source.has_mcmc_gspphot,gaia_source.has_mcmc_msc,gaia_source.teff_gspphot,
						gaia_source.logg_gspphot,gaia_source.mh_gspphot,gaia_source.distance_gspphot,gaia_source.azero_gspphot,
						gaia_source.ag_gspphot,gaia_source.ebpminrp_gspphot
						FROM gaiadr3.gaia_source 
						WHERE 
						gaiadr3.gaia_source.designation = '{DESI}'"""

	example = example_base.format(DESI=designation)

	job = Gaia.launch_job_async(example)

	return job.get_results()

def queryGaiaDR2_byID(designation):
	example_base = """SELECT TOP 1 gaia_source.designation,gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,
						gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,
						gaia_source.phot_g_mean_mag,
						gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error
						FROM gaiadr2.gaia_source 
						WHERE 
						gaiadr2.gaia_source.designation = '{DESI}'"""

	example = example_base.format(DESI=designation)

	job = Gaia.launch_job_async(example)

	return job.get_results()

def EXCAM_PointingCheck(ra,dec):
	#coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.degree), frame='fk5')
	coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')

	#degRA = coord.ra.degree
	#degDec = coord.dec.degree

	#EXCAM_restrict_arcsec = 12.0*u.arcsec
	#EXCAM_restrict_degree = EXCAM_restrict_arcsec.to(u.degree).value
	#print(EXCAM_restrict_degree)

	#example_base = """SELECT TOP 50 gaia_source.designation,gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.pmra,gaia_source.pmdec,gaia_source.ruwe,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.non_single_star,gaia_source.has_xp_continuous,gaia_source.has_xp_sampled,gaia_source.has_rvs,gaia_source.has_epoch_photometry,gaia_source.has_epoch_rv,gaia_source.has_mcmc_gspphot,gaia_source.has_mcmc_msc,gaia_source.teff_gspphot,gaia_source.logg_gspphot,gaia_source.mh_gspphot,gaia_source.distance_gspphot,gaia_source.azero_gspphot,gaia_source.ag_gspphot,gaia_source.ebpminrp_gspphot,
	#					'17:26:35--48:03:42' AS target_id,{RA} AS target_ra,{DEC} AS target_dec,
	#					DISTANCE(
	#					POINT('ICRS', gaiadr3.gaia_source.ra, gaiadr3.gaia_source.dec), 
	#					POINT('ICRS', {RA}, {DEC}) ) AS 'target_separation (deg)'
	#					FROM gaiadr3.gaia_source 
	#					WHERE 
	#					CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),
	#					CIRCLE('ICRS',{RA},{DEC},{excam_restrict})
	#					)=1 
	#					ORDER BY target_id, 'target_separation (deg)'"""


	#example = example_base.format(RA=str(degRA),DEC=str(degDec),excam_restrict=str(EXCAM_restrict_degree))

	#example = """SELECT TOP 2000 gaia_source.designation,gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.pmra,gaia_source.pmdec,gaia_source.ruwe,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.non_single_star,gaia_source.has_xp_continuous,gaia_source.has_xp_sampled,gaia_source.has_rvs,gaia_source.has_epoch_photometry,gaia_source.has_epoch_rv,gaia_source.has_mcmc_gspphot,gaia_source.has_mcmc_msc,gaia_source.teff_gspphot,gaia_source.logg_gspphot,gaia_source.mh_gspphot,gaia_source.distance_gspphot,gaia_source.azero_gspphot,gaia_source.ag_gspphot,gaia_source.ebpminrp_gspphot,'2.2990654282325518-59.14897576410493' AS target_id,2.2990654282325518 AS target_ra,59.14897576410493 AS target_dec,
	#				DISTANCE(
	#				POINT('ICRS', gaiadr3.gaia_source.ra, gaiadr3.gaia_source.dec), 
	#				POINT('ICRS', 2.2990654282325518, 59.14897576410493) 
	#				) AS "target_separation (deg)"

	#				FROM gaiadr3.gaia_source 
	#				WHERE 
	#				CONTAINS(
	#				POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),
	#				CIRCLE('ICRS',2.2990654282325518,59.14897576410493,0.0033333333333333335)
	#				)=1 ORDER BY target_id, 'target_separation (deg)'"""

	#print(example)

	job = Gaia.cone_search_async(coordinate=coord, radius=u.Quantity(12.0,u.arcsec))

	#job = Gaia.launch_job_async(example)

	return job.get_results()

def EXCAM_Standards(filebase):
	#referencelist = LoadExistingDatabase(filename='PolStandards_Trialv2.csv')
	#referencelist = pandas.read_csv('UnpolStandards_Trial.csv')
	referencelist = pandas.read_csv(filebase+'_v2.csv')


	Simbad_RA = referencelist['ra'].values #Simbad RA in degrees (J2000)

	Simbad_Dec = referencelist['dec'].values #Simbad Dec in degrees (J2000)

	names = referencelist['main_id'].values #main identifiers

	GaiaDR3 = referencelist['Gaia_DR3'].values #Gaia DR3 identifiers

	DR3_RA = referencelist['ra_2016'].values #Gaia DR3 RA in degrees (J2016)

	DR3_Dec = referencelist['dec_2016'].values #Gaia DR3 Dec in degrees (J2016)

	Target_G16 = referencelist['Gmag_2016'].values #

	Target_G = referencelist['G'].values

	Target_V = referencelist['V'].values

	referencelist.insert(len(referencelist.columns),"EXCAM_PASS",[False for _ in range(len(names))])

	for i in np.arange(len(names)):
		print('Checking ',names[i])
		if np.isnan(DR3_RA[i]) == False and np.isnan(DR3_Dec[i]) == False:
			tableresults = EXCAM_PointingCheck(DR3_RA[i],DR3_Dec[i])
			
		elif np.isnan(DR3_RA[i]) == True or np.isnan(DR3_Dec[i]) == True:
			print('WARNING: USING SIMBAD COORDINATES')
			tableresults = EXCAM_PointingCheck(Simbad_RA[i],Simbad_Dec[i])
		QueryGmags = tableresults['phot_g_mean_mag']
		print(QueryGmags)
		if len(QueryGmags)>1:
			if np.isnan(Target_G16[i]) == False:
				print('Catalog G: ',Target_G16[i])
				failure = np.where(QueryGmags<Target_G16[i]+3)[0] #exclude first row since it is the target
				if len(failure) > 1:
					print('This fails EXCAM Acquisition!')
					EXCAM_test = False
				else:
					EXCAM_test = True
			elif np.isnan(Target_G16[i]) == True and np.isnan(Target_G[i]) == False:
				print('Catalog G: ',Target_G[i])
				failure = np.where(QueryGmags<Target_G[i]+3)[0] #exclude first row since it is the target
				if len(failure) > 1:
					print('This fails EXCAM Acquisition!')
					EXCAM_test = False
				else:
					EXCAM_test = True
			elif np.isnan(Target_G16[i]) == True and np.isnan(Target_G[i]) == True:
				print('Catalog V: ',Target_V[i])
				failure = np.where(QueryGmags<Target_V[i]+3)[0] #exclude first row since it is the target
				if len(failure) > 1:
					print('This fails EXCAM Acquisition!')
					EXCAM_test = False
				else:
					EXCAM_test = True
		else:
			EXCAM_test = True
		referencelist.loc[i,"EXCAM_PASS"] = EXCAM_test


		
		print('Completed ',names[i])
	print(referencelist["EXCAM_PASS"])



	#df = fullTable.to_pandas()

	#referencelist.to_csv('PolStandards_Trialv3_EXCAM.csv')
	#referencelist.to_csv('UnpolStandards_Trialv3_EXCAM.csv')
	referencelist.to_csv(filebase+'_EXCAMv2.csv')

def addPrebuiltColumn2New(prebuiltTable,newTable,rowIndex,columnName="P% (U)"):
	newTable[columnName] = prebuiltTable[columnName].values[rowIndex]

	return newTable










