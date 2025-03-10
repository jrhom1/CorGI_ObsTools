import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from matplotlib.patches import Patch,Rectangle
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import matplotlib as mpl

mpl.rc('font',family='Times')
plt.rcParams['xtick.labelsize']=15
plt.rcParams['ytick.labelsize']=15

def offaxisPSF(save=False):
	hlc_profile = 'HLC_prof.csv'
	#hlc_profile = fits.getdata('')
	separation = np.loadtxt(hlc_profile,usecols=0,delimiter=',')
	contrast = np.loadtxt(hlc_profile,usecols=1,delimiter=',')

	fig1,ax1 = plt.subplots(1,1,figsize=(8,6))

	ax1.plot(separation,-2.5*np.log10(contrast),label='Unocculted PSF Profile')
	ax1.plot(separation,np.ones_like(contrast)*17.5,linestyle='--',label=r'10$^{-7}$ Contrast')
	ax1.plot(separation,np.ones_like(contrast)*20.0,linestyle='-.',label=r'10$^{-8}$ Contrast')
	ax1.set_xscale('log')
	#ax1.set_yscale('log')
	ax1.set_xlabel('Separation [arcsec]',fontsize=18)
	ax1.set_ylabel(r'$\Delta V$ (relative to source)',fontsize=18)
	ax1.yaxis.set_inverted(True)
	ax1.set_xlim([7e-3,1e2])
	#ax1.set_title('Unocculted PSF Profile',fontsize=18)
	ax1.legend(loc='best',fontsize=16)
	plt.show()
	if save == True:
		fig1.savefig('HLC_PSFprof.pdf')

def offaxisPSF_Demo(primarymag=0,comp_sep=0,comp_mag=0,save=False):
	'''
	offaxisPSF_Demo: Creates a plot that shows how much flux the PSF profile of an off-axis companion can introduce into the dark hole of a given primary star magnitude.
	Inputs:
		primarymag: Default: 0. Primary component apparent V magnitude.
		comp_sep: Default: 0. Off-axis component separation in arcseconds.
		comp_mag: Default: 0. Off-axis component apparent V magnitude.
		save: Default: False. Boolean for whether or not to save the figure.

	'''

	hlc_profile = 'HLC_prof.csv'
	#hlc_profile = fits.getdata('')
	separation = np.loadtxt(hlc_profile,usecols=0,delimiter=',')
	contrast = np.loadtxt(hlc_profile,usecols=1,delimiter=',')

	lefthalfseps = separation[::-1]*-1.0
	lefthalfcontrast = contrast[::-1]

	#print(lefthalfseps)

	fullseps = np.append(lefthalfseps,separation)
	fullcontrast = np.append(lefthalfcontrast,contrast)

	c_v = 'dodgerblue'
	c_bbvis = 'cadetblue'
	c_band3 = 'goldenrod'
	c_band4 = 'orange'

	legend_elements1 = [Patch(color='dodgerblue', alpha=0.5,label='HLC Band 1'),
    					#Patch(color='goldenrod',alpha=0.5, label='SPC-SPEC Band 3'),
    					Patch(color='orange', alpha=0.5,label='SPC-WFOV Band 4'),
    					#Line2D([0], [0], marker=None, linestyle='--',color='black', label='Pol Measured',
                        #  markerfacecolor='black', markersize=10)
    					]

	fig1,ax1 = plt.subplots(1,1,figsize=(8,6))

	ax1.plot(fullseps+comp_sep,-2.5*np.log10(fullcontrast)+comp_mag,label='Unocculted PSF Profile',color='navy')
	ax1.plot(fullseps,np.ones_like(fullcontrast)*(primarymag+17.5),linestyle='--',label=r'10$^{-7}$ Contrast',color='k')
	ax1.plot(fullseps,np.ones_like(fullcontrast)*(primarymag+20.0),linestyle='-.',label=r'10$^{-8}$ Contrast',color='k')
	if comp_sep>0:

		ax1.text(-5,(primarymag+17.5),r'$10^{-7}$ Dark Hole',color='k',horizontalalignment='left',va='bottom',fontsize=14)
		ax1.text(-5,(primarymag+20.0),r'$10^{-8}$ Dark Hole',color='k',horizontalalignment='left',va='bottom',fontsize=14)
	elif comp_sep<0:
		ax1.text(5,(primarymag+17.5),r'$10^{-7}$ Dark Hole',color='k',horizontalalignment='right',va='bottom',fontsize=14)
		ax1.text(5,(primarymag+20.0),r'$10^{-8}$ Dark Hole',color='k',horizontalalignment='right',va='bottom',fontsize=14)

	ax1.axvspan(0.110,0.450,alpha=0.5,color='dodgerblue')
	ax1.axvspan(-0.450,-0.110,alpha=0.5,color='dodgerblue')
	#ax1.axvspan(0.180,0.550,alpha=0.5,color='goldenrod')
	#ax1.axvspan(-0.550,-0.180,alpha=0.5,color='goldenrod')
	ax1.axvspan(0.450,1.40,alpha=0.5,color='orange')
	ax1.axvspan(-1.40,-0.450,alpha=0.5,color='orange')
	ax1.set_xscale('symlog',linthresh=1)
	ax1.set_xticks([-10,-5,-1,-0.5,-0.1, 0.1, 0.5, 1, 5,10])
	ax1.tick_params(axis='both', which='major', labelsize=12)

	#ax1.set_yscale('log')
	ax1.set_xlabel('Separation [arcsec]',fontsize=18)
	ax1.set_ylabel(r'$V$',fontsize=18)
	ax1.yaxis.set_inverted(True)
	#ax1.set_xlim([7e-3,1e2])
	ax1.set_xlim([-10,10])
	ax1.set_ylim([25+primarymag,primarymag])
	ax1.set_title(r'Ref. Star $V$='+str(primarymag)+r', Comp. $V$='+str(comp_mag)+', Sep = '+str(comp_sep)+' arcsec',fontsize=18)
	#ax1.legend(handles=legend_elements1,fontsize=14,bbox_to_anchor=(1.0,0.8))
	ax1.legend(handles=legend_elements1,fontsize=14,loc=1)
	ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2g'))
	fig1.tight_layout()
	if save==True:
		fig1.savefig('OffAxis_HostV'+str(primarymag)+'_CompV'+str(comp_mag)+'_Sep'+str(comp_sep)+'.jpg')
	plt.show()






