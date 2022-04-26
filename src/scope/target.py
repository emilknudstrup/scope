#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
from astropy.time import Time
import astropy.units as u
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import numpy as np
from astropy import constants as const

class Target(object):
	'''Target class.
	
	Class to represent orbital/system parameters.
	
	Attributes
	----------
	per : float
		Period in days.

	T0 : float
		Mid-transit time. Format Julian Date. Scale utc.

	duration : float
		Duration of transit in hours.

	RA : float
		Right ascension in degrees.

	Dec : float
		Declination in degrees.

	Vmag : float
		V magnitude.
	
	ID : str
		Name of planet.

	starID : str
		Name of star.

	full_transit : bool
		Whether to observe a full transit, or if a partial transit is okay.

	Methods
	-------
	make_target()
	'''
	
	def __init__(self,per,T0,duration,Vmag,RA,Dec,
			  ID='Planet name',starID='Star name',
			  b=0.0,vsini=0.0):
		'''Initialize target.

		Parameters
		----------
			per : float
				Period in days.

			T0 : float
				Mid-transit time. Format Julian Date. Scale utc.

			duration : float
				Duration of transit in hours.

			RA : float
				Right ascension in degrees.

			Dec : float
				Declination in degrees.

			Vmag : float
				V magnitude.
			
			ID : str
				Name of planet.

			starID : str
				Name of star.

			full_transit : bool
				Whether to observe a full transit, or if a partial transit is okay.

		'''
		self.per = per
		self.T0 = T0
		self.duration = duration
		self.Vmag = Vmag
		self.RA = RA
		self.Dec = Dec
		self.ID = 'Planet name'
		self.starID = 'Star name'
		self.b = b
		self.vsini = vsini


	def makeTarget(self):
		'''Add units to parameters.
		

		'''
		return self.per*u.day, Time(self.T0,format='jd'), self.duration*u.hour, self.RA*u.deg, self.Dec*u.deg

class GetTarget(object):
	
	ESSENTIALS = ['ra','dec','pl_orbper','pl_tranmid','pl_ratror','st_vsin','sy_vmag','pl_orbsmax','pl_radj','st_rad','pl_orbincl']#,'pl_trandur'
	
	def __init__(self,ADDITIONAL=['pl_trueobliq','pl_massj','st_teff','pl_trandur','pl_orbeccen','pl_imppar'],tab=None):
		self.ADDITIONAL = ADDITIONAL
		self.tab = None
		self.plDict = {}
		self.targets = []
	
	def getTargets(self,table='ps'):
		self.table = table

		tab = NasaExoplanetArchive.query_criteria_async(table=self.table)
		
		skip = []
		for pl in tab['pl_name']:
			if pl in skip:
				pass
			else:
				skip.append(pl)
				entries = np.where(pl == tab['pl_name'])[0]
				subDict = {}
				values = [np.nan]*len(self.ESSENTIALS)
				for ii, ess in enumerate(self.ESSENTIALS):
					for entry in entries:
						val = tab[ess][entry]
						
						if np.isfinite(val):
							subDict[ess] = val
							values[ii] = val
							try:
								unc1 = tab[ess+'err1'][entry]
								if np.isfinite(unc1): subDict[ess+'err1'] = unc1
								unc2 = tab[ess+'err2'][entry]
								if np.isfinite(unc2): subDict[ess+'err2'] = unc2
							except KeyError:
								pass
							continue
				subDict['hostname'] = tab['hostname'][entry]
				for ii, add in enumerate(self.ADDITIONAL):
					for entry in entries:
						val = tab[add][entry]
						
						if np.isfinite(val):
							subDict[add] = val
							values[ii] = val
							try:
								unc1 = tab[add+'err1'][entry]
								if np.isfinite(unc1): subDict[add+'err1'] = unc1
								unc2 = tab[add+'err2'][entry]
								if np.isfinite(unc2): subDict[add+'err2'] = unc2
							except KeyError:
								pass
							continue
				if all(np.isfinite(values)): 
					aAU = subDict['pl_orbsmax']
					Rp = subDict['pl_radj']
					Rs = subDict['st_rad']
					inc = subDict['pl_orbincl']
					per = subDict['pl_orbper']
					aR = aAU*const.au/(Rs*const.R_sun)
					aR = aR.value
					b = subDict['pl_imppar']
					if np.isnan(b):
						b = aR*np.cos(inc*np.pi/180.)
						subDict['pl_imppar'] = b
					#rp = Rp*const.R_jup/Rs*const.R_sun
					rp = subDict['pl_ratror']
					nom = np.sqrt((1+rp)**2-b**2)
					den = aR*np.sin(inc*np.pi/180.)
					dur = subDict['pl_trandur']
					if np.isnan(dur):
						dur = 24*per/np.pi*np.arcsin(nom/den)
						subDict['pl_trandur'] = dur

					
					
					self.plDict[pl] = subDict


	def getTarget(self,name):
		self.name = name
		
		tab = NasaExoplanetArchive.query_object(self.name) 
		keys = tab.keys()
		arr = tab.as_array()
		rows = arr.shape[0]
		subDict = {}
		idxs = np.argsort(Time(tab['rowupdate'],scale='utc',format='iso').jd)
		for key in keys:
			for idx in idxs:

				val = arr[key][idx]
				try:
					if np.isfinite(val): subDict[key] = val
				except TypeError:
					subDict[key] = val
					pass
				try:
					unc1 = tab[key+'err1'][idx]
					if np.isfinite(unc1): subDict[key+'err1'] = unc1
					unc2 = tab[key+'err2'][idx]
					if np.isfinite(unc2): subDict[key+'err2'] = unc2
				except KeyError:
					pass
		self.plDict[name] = subDict
		
	
	def createTargetlist(self):
		pls = self.plDict.keys()
		for pl in pls:
			per = self.plDict[pl]['pl_orbper']
			Vmag = self.plDict[pl]['sy_vmag']
			T0 = self.plDict[pl]['pl_tranmid']
			duration = self.plDict[pl]['pl_trandur']
			ra = self.plDict[pl]['ra']
			dec = self.plDict[pl]['dec']
			target = Target(per,T0,duration,Vmag,ra,dec)
			target.ID = pl
			target.starID = self.plDict[pl]['hostname']
			target.b = self.plDict[pl]['pl_imppar']
			target.vsini = self.plDict[pl]['st_vsin']
			self.targets.append(target)
		
#tab = GetTarget()
#targets = tab.getTarget('K2-261 b')

#cols = tab.fieldnames
#use = tab['pl_name'] == 'Kepler-11 c'
