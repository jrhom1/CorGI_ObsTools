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
	mag855 = queryTap_byVmagAsym(10.9,12.0)

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

	ax1.plot(ra_rad855,dec_rad855,marker='o',linestyle='None',color='indianred',markersize=10,label='V>10.9')

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
	mag855 = queryTap_byVmag(8.55,cone=0.05)

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

	ax1.plot(ra_rad855,dec_rad855,marker='o',linestyle='None',color='indianred',markersize=10,label='V~8.55')

	ax1.set_title('Commissioning Flats, Heliocentric True Ecliptic, J2027, eq=2027',fontsize=18)

	ax1.set_xlabel('l [deg]',fontsize=18)

	ax1.set_ylabel('b [deg]',fontsize=18)

	ax1.legend(bbox_to_anchor=(1.2,0),loc='lower right',fontsize=16)

	fig1.tight_layout()

	if plotsave == True:

		fig1.savefig('CommissioningFlatStandards.jpg')

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

	magList = np.array([5.5,8.5,11.5,13.5]) #Array of V magnitudes that we want to query
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
		mag = queryTap_byVmag(magList[i],cone=0.05)
    	#The following lines fill any missing values in your table query with the numbers in parentheses
    	#mag["plx_value"] = mag["plx_value"].filled(0.00001)
    	#mag["pmra"] = mag["pmra"].filled(0.0)
    	#mag["pmdec"] = mag["pmdec"].filled(0.0)
    	#mag["rvz_radvel"] = mag["rvz_radvel"].filled(0.0)
		magDist = 1000./mag["plx_value"].value #Calculate the system distance in pc
		ics = SkyCoord(ra=mag['ra'],dec=mag['dec'],unit=(u.deg,u.deg),frame='icrs',
		distance=magDist*u.pc,pm_ra_cosdec=mag["pmra"],pm_dec=mag["pmdec"],
		radial_velocity=mag["rvz_radvel"],obstime=astroTime('J2000',format='jyear_str')) #Set up a SkyCoord object. Requires the values from Simbad Query.

		mag['CalType'] = ['ImagCorr_'+str(magList[i]) for j in np.arange(50)]

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










