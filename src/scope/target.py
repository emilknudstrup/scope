#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
from astropy.time import Time
import astropy.units as u
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import numpy as np
from astropy import constants as const
import pickle

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
			  b=np.nan,vsini=np.nan,Teff=np.nan,rp=np.nan,
			  RMamp=np.nan,lam=np.nan,psi=np.nan,a=np.nan):
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
		self.Teff = Teff
		self.vsini = vsini
		self.rp = rp
		self.a = a
		self.RMamp = RMamp
		self.lam = lam
		self.psi = psi


	def makeTarget(self):
		'''Add units to parameters.
		

		'''
		return self.per*u.day, Time(self.T0,format='jd'), self.duration*u.hour, self.RA*u.deg, self.Dec*u.deg


class GetTarget(object):
	
	ESSENTIALS = ['ra','dec','pl_orbper','pl_tranmid','pl_ratror','st_vsin','sy_vmag','pl_orbsmax','pl_radj','st_rad','pl_orbincl']
	
	def __init__(self,ADDITIONAL=['pl_trueobliq','st_teff','pl_trandur','pl_orbeccen','pl_imppar']):
		self.ADDITIONAL = ADDITIONAL
		#self.tab = None
		self.plDict = {}
		self.targets = []
	
	def getRMTargets(self,table='pscomppars',allkeys=True):
		table = table

		tab = NasaExoplanetArchive.query_criteria_async(table=table)
		
		skip = []
		for pl in tab['pl_name']:
			if table != 'pscomppars':
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
								try:
									unc1 = tab[add+'err1'][entry]
									if np.isfinite(unc1): subDict[add+'err1'] = unc1
									unc2 = tab[add+'err2'][entry]
									if np.isfinite(unc2): subDict[add+'err2'] = unc2
								except KeyError:
									pass
								continue
			else:
				subDict = {}
				entry = np.where(pl == tab['pl_name'])[0]
				subDict['hostname'] = tab['hostname'][entry][0]
				values = [np.nan]*len(self.ESSENTIALS)
				for ii, ess in enumerate(self.ESSENTIALS):
					val = tab[ess][entry].data
					try:
						subDict[ess] = val[0]
						values[ii] = val[0]
						try:
							unc1 = tab[ess+'err1'][entry].data[0]
							if np.isfinite(unc1): subDict[ess+'err1'] = unc1
							unc2 = tab[ess+'err2'][entry].data[0]
							if np.isfinite(unc2): subDict[ess+'err2'] = unc2
						except KeyError:
							pass
					except (TypeError,IndexError):
						continue
				for ii, add in enumerate(self.ADDITIONAL):
					val = tab[add][entry].data
					try:
						subDict[add] = val[0]
						try:
							unc1 = tab[add+'err1'][entry].data[0]
							if np.isfinite(unc1): subDict[add+'err1'] = unc1
							unc2 = tab[add+'err2'][entry].data[0]
							if np.isfinite(unc2): subDict[add+'err2'] = unc2
						except KeyError:
							pass
					except (TypeError,IndexError):
						continue

				if all(np.isfinite(values)): 
					aAU = subDict['pl_orbsmax']
					Rp = subDict['pl_radj']
					Rs = subDict['st_rad']
					inc = subDict['pl_orbincl']
					per = subDict['pl_orbper']
					#aR = aAU*const.au/(Rs*const.R_sun)
					#aR = aR.value
					#subDict['pl_ratdor'] = aR
					aR = subDict['pl_ratdor'] 
					try:
						b = subDict['pl_imppar']
					except KeyError:
						b = np.nan
					if np.isnan(b):
						if inc > 90: inc = 90-(inc-90)
						b = aR*np.cos(inc*np.pi/180.)
						subDict['pl_imppar'] = b
					#rp = Rp*const.R_jup/Rs*const.R_sun
					subDict['pl_RMamp'] = 0.7*subDict['st_vsin']*np.sqrt(1-b**2)*np.power(subDict['pl_ratror'],2)
					try:
						dur = subDict['pl_trandur']
					except KeyError:
						dur = np.nan
					if np.isnan(dur):
						rp = subDict['pl_ratror'] 
						nom = np.sqrt((1+rp)**2-b**2)
						den = aR*np.sin(inc*np.pi/180.)
						dur = 24*per/np.pi*np.arcsin(nom/den)
						subDict['pl_trandur'] = dur

					if allkeys:
						keys = tab.fieldnames
						for key in keys:
							#if (key not in self.ESSENTIALS) & (key not in self.ADDITIONAL):
							try:
								val = subDict[key]
							except KeyError:
								subDict[key] = np.nan
								val = tab[key][entry].data
								try:
									subDict[key] = val[0]
								except (TypeError,IndexError):
									pass
						
					
					self.plDict[pl] = subDict

	def getTargets(self,table='pscomppars',requirements={'sy_vmag' : [0,16]}):
		self.table = table
		self.requirements = requirements
		tab = NasaExoplanetArchive.query_criteria_async(table=self.table)
		self.keys = tab.fieldnames
		for pl in tab['pl_name']:
			subDict = {}
			idx = np.where(pl == tab['pl_name'])[0]
			for key in self.keys:
				val = tab[key][idx]
				subDict[key] = val
				# val = tab[key][idx]
				# try:
				# 	if np.isfinite(val): subDict[key] = val
				# except TypeError:
				# 	subDict[key] = val
				# 	pass
				# try:
				# 	unc1 = tab[key+'err1'][idx]
				# 	if np.isfinite(unc1): subDict[key+'err1'] = unc1
				# 	unc2 = tab[key+'err2'][idx]
				# 	if np.isfinite(unc2): subDict[key+'err2'] = unc2
				# except KeyError:
				# 	pass

			req = True
			for key in self.requirements.keys(): 
				low = requirements[key][0]
				high = requirements[key][1]
				try:
					val = subDict[key]
					if (low > val) or (val > high):
						req = False
				except KeyError:
					req = False

			if req:	self.plDict[pl] = subDict



	def byName(self,name):
		self.name = name
		print('Query to NASA Exoplanet Archive for {}.'.format(name))
		tab = NasaExoplanetArchive.query_object(self.name) 
		self.keys = tab.keys()
		arr = tab.as_array()
		rows = arr.shape[0]
		subDict = {}
		idxs = np.argsort(Time(tab['rowupdate'],scale='utc',format='iso').jd)
		if len(idxs):
			for key in self.keys:
				for idx in idxs:
					val = arr[key][idx]
					try:
						if np.isfinite(val): subDict[key] = val
					except TypeError:
						subDict[key] = val
						pass
					# try:
					# 	unc1 = tab[key+'err1'][idx]
					# 	if np.isfinite(unc1): subDict[key+'err1'] = unc1
					# 	unc2 = tab[key+'err2'][idx]
					# 	if np.isfinite(unc2): subDict[key+'err2'] = unc2
					# except KeyError:
					# 	pass
				if (key in self.ESSENTIALS):
					try:
						subDict[key]
					except KeyError:
						#subDict[key] = np.nan
						print('Essential parameter {} appears to be nan.'.format(key))
			self.plDict[name] = subDict
		
	
	def createTransitTargetlist(self):
		pls = self.plDict.keys()
		for pl in pls:
			print('Appending target {}'.format(pl))
			per = self.plDict[pl]['pl_orbper']
			Vmag = self.plDict[pl]['sy_vmag']
			T0 = self.plDict[pl]['pl_tranmid']
			duration = self.plDict[pl]['pl_trandur']
			ra = self.plDict[pl]['ra']
			dec = self.plDict[pl]['dec']
			target = Target(per,T0,duration,Vmag,ra,dec)
			target.ID = pl
			try:
				target.starID = self.plDict[pl]['hostname']
			except KeyError:
				print('No name for host.')
			try:
				target.b = self.plDict[pl]['pl_imppar']
			except KeyError:
				print('Impact parameter missing.')
			try:
				target.vsini = self.plDict[pl]['st_vsin']
			except KeyError:
				print('vsini missing.')
			try:
				target.lam = self.plDict[pl]['pl_projobliq']
			except KeyError:
				print('lambda missing.')
			try:
				target.psi = self.plDict[pl]['pl_trueprojobliq']
			except KeyError:
				print('psi missing.')
			try:
				target.Teff = self.plDict[pl]['st_teff']
			except KeyError:
				print('Teff missing.')
			try:
				target.rp = self.plDict[pl]['pl_ratror']
			except KeyError:
				print('Rp/Rs missing.')
			try:
				target.a = self.plDict[pl]['pl_ratdor']
			except KeyError:
				print('a/Rs missing.')

			self.targets.append(target)

	def createTargetlist(self):
		pls = self.plDict.keys()
		for pl in pls:
			Vmag = self.plDict[pl]['sy_vmag']
			ra = self.plDict[pl]['ra']
			dec = self.plDict[pl]['dec']
			target = Target(0.0,0.0,0.0,Vmag,ra,dec)
			target.ID = pl
			target.starID = self.plDict[pl]['hostname']
			self.targets.append(target)
	

	def writeTargetDict(self,path):
		with open(path,'wb') as file:
			pickle.dump(self.plDict,file)

	def readTargetDict(self,path):
		with open(path,'rb') as file:
			self.plDict = pickle.load(file)

	def printPlanet(self,plname=None):
		
		if plname == None: plname = list(self.plDict.keys())[0]


		#print('\nFound planet {} in {}.\n'.format(plname,cname))
		#plname = self.plDict[plname]['pl_name']

		print(plname+':')
		disc = self.plDict[plname]['discoverymethod']
		starname = self.plDict[plname]['hostname']
		vsini = self.plDict[plname]['st_vsin']
		Vmag = self.plDict[plname]['sy_vmag']
		per = self.plDict[plname]['pl_orbper']
		Rp = self.plDict[plname]['pl_ratror']
		RJ = self.plDict[plname]['pl_radj']
		teff = int(self.plDict[plname]['st_teff'])
		ra, dec = self.plDict[plname]['ra'], self.plDict[plname]['dec']
		T0 = self.plDict[plname]['pl_tranmid'] 
		dur = self.plDict[plname]['pl_trandur']
		per = self.plDict[plname]['pl_orbper']
		aAU = self.plDict[plname]['pl_orbsmax']
		inc = self.plDict[plname]['pl_orbincl']
		b = self.plDict[plname]['pl_imppar']
		Rs = self.plDict[plname]['st_rad']
		aR = self.plDict[plname]['pl_ratdor']

		#rmA = RM(Rp*pc.RJ/(Rs*pc.RSun),b,vsini,90.)*1000

		print('Rp/Rs = {:0.2f}'.format(Rp))
		print('a/Rs = {:0.2f}'.format(aR))
		print('P = {:0.1f} d'.format(per))
		print('b = {:0.1f}'.format(b))
		print('e = {:0.2f}'.format(self.plDict[plname]['pl_orbeccen']))


		#print('\nRp = {:0.2f} RE'.format(RJ*pc.RJ/pc.REarth))
		print('Rp = {:0.2f} RJ'.format(RJ))

		#print('\nRM-amplitude = {:0.1f} m/s'.format(rmA))
		try:
			print('lambda = {}+/-{} deg'.format(self.plDict[plname]['pl_projobliq'],self.plDict[plname]['pl_projobliqerr1']))
		except KeyError:
			print('lambda = nan')

		print('\n'+starname+':')
		print('Teff = {:d} K'.format(teff))
		print('Rs = {:0.1f} RSun'.format(Rs))
		print('vsini = {:0.1f} km/s'.format(vsini))
		print('Vmag = {:0.1f}'.format(Vmag))
		print('SpT, '.format(self.plDict[plname]['st_spectype']))

		# target = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
		# print('\nR.A. {:d}:{:d}:{:0.2f}, Dec. {:d}:{:d}:{:0.2f}'.format(\
		#   int(target.ra.hms.h),int(target.ra.hms.m),target.ra.hms.s,\
		#   int(target.dec.dms.d),abs(int(target.dec.dms.m)),abs(target.dec.dms.s)))
		# print()

	def prepare(self,syname):
		pls = ['b','c','d','e','f','g','h','i','j']

		self.byName(syname)

		npls = self.plDict[syname]['sy_pnum']

		#if plname == None: plname = list(self.plDict.keys())[0]
		self.system = {
			'Teff' : self.plDict[syname]['st_teff'],#err1,err2
			'logg' : self.plDict[syname]['st_logg'],
			'FeH' : self.plDict[syname]['st_met'],
			'vsini' : self.plDict[syname]['st_vsin'],

		}

		self.planets = {
			'pls' : []
		}

		fnpars = ['P','T0','rp','aR','inc','K']
		NASApars = ['pl_orbper','pl_tranmid','pl_ratror','pl_ratdor','pl_orbincl','pl_rvamp']

		for ii in range(npls):
			pl = pls[ii]
			plname = syname+' ' + pl
			self.byName(plname)
			self.planets['pls'].append(pl)

			self.planets[pl] = {}
			for jj, fnpar in enumerate(fnpars):
				NASApar = NASApars[jj]
				try:
					self.planets[pl][fnpar] = self.plDict[plname][NASApar]
				except KeyError:
					print('{} not defined setting to np.nan.'.format(fnpar))
					self.planets[pl][fnpar] = np.nan
			# self.planets[pl]['P'] = self.plDict[plname]['pl_orbper']
			# self.planets[pl]['T0'] = self.plDict[plname]['pl_tranmid']
			# self.planets[pl]['rp'] = self.plDict[plname]['pl_ratror']
			# self.planets[pl]['aR'] = self.plDict[plname]['pl_ratdor']
			# self.planets[pl]['inc'] = self.plDict[plname]['pl_orbincl']
			# self.planets[pl]['K'] = self.plDict[plname]['pl_rvamp']
			self.planets[pl]['law'] = 'quadratic'
			self.planets[pl]['cs'] = [0.3,0.2]

			try:
				e = self.plDict[plname]['pl_orbeccen']
				w = self.plDict[plname]['pl_orblper']
			except KeyError:
				e = 0.0
				w = 90.
	 
			self.planets[pl]['ecc'] = e
			self.planets[pl]['w'] = w


		# self.planets = {
		# 	'pls': [self.plDict[plname]['pl_letter']],
		# 	self.plDict[plname]['pl_letter']  : {
		# 	'P'  : self.plDict[plname]['pl_orbper'],
		# 	'T0' : self.plDict[plname]['pl_tranmid'],
		# 	'rp' : self.plDict[plname]['pl_ratror'],
		# 	'aR' : self.plDict[plname]['pl_ratdor'],
		# 	'inc': self.plDict[plname]['pl_orbincl'],
		# 	'ecc': e,
		# 	'w'  : w,
		# 	'law': 'quadratic',
		# 	'cs' : [0.3,0.2]
		# 	}

		# }

	def prepareFromtracit(self,df,npls=1):
		pls = ['b','c','d','e','f','g','h','i','j']

		self.planets = {
			'pls' : []
		}

		fnpars = ['P','T0','rp','aR','inc','K']
		tpars = ['P','T0','Rp_Rs','a_Rs','inc','K']

		for ii in range(npls):
			pl = pls[ii]

			self.planets['pls'].append(pl)
			self.planets[pl] = {}
			for jj, fnpar in enumerate(fnpars):
				tpar = tpars[jj] + '_' + pl
				try:
					self.planets[pl][fnpar] = float(df[tpar][4])#self.plDict[plname][NASApar]
				except KeyError:
					print('{} not defined setting to np.nan.'.format(fnpar))
					self.planets[pl][fnpar] = np.nan
			try:
				e = self.planets[pl]['e_'+pl]
				w = self.planets[pl]['w_'+pl]
			except KeyError:
				e = 0.0
				w = 90.
			self.planets[pl]['ecc'] = e
			self.planets[pl]['w'] = w

			self.planets[pl]['law'] = 'quadratic'
			self.planets[pl]['cs'] = [0.3,0.2]

	# 		try:
	# 			e = self.plDict[plname]['pl_orbeccen']
	# 			w = self.plDict[plname]['pl_orblper']
	# 		except KeyError:
	# 			e = 0.0
	# 			w = 90.
	 
	# 		self.planets[pl]['ecc'] = e
	# 		self.planets[pl]['w'] = w

	def prepare4tracit(self,syname=None,par=None,n_phot=1,n_spec=1):
		if not par:
			import tracit
			par = tracit.par_struct(n_planets=len(planets),n_phot=n_phot,n_spec=n_spec)

		try:
			planets = self.planets['pls']
		except AttributeError:
			if not syname: print('You need to give the name of the system,\nwhen `prepare` has not been initialized.')
			self.prepare(syname)
			planets = self.planets['pls']

		#par = {}
		par['Planets'] = planets

		#pars = ['P','T0','inc','ecc','K','']

		for pl in par['Planets']:
			par['P_'+pl]['Value'] = self.planets[pl]['P']
			par['P_'+pl]['Prior_vals'] = [self.planets[pl]['P'],0.5,0,self.planets[pl]['P']+10]
			par['T0_'+pl]['Value'] = self.planets[pl]['T0']
			par['T0_'+pl]['Prior_vals'] = [self.planets[pl]['T0'],0.5,self.planets[pl]['T0']-10,self.planets[pl]['T0']+10]
			par['a_Rs_'+pl]['Value'] = self.planets[pl]['aR']
			par['a_Rs_'+pl]['Prior_vals'] = [self.planets[pl]['aR'],0.5,0,self.planets[pl]['aR']+10]
			par['Rp_Rs_'+pl]['Value'] = self.planets[pl]['rp']
			par['Rp_Rs_'+pl]['Prior_vals'] = [self.planets[pl]['rp'],0.01,0,1.0]
			par['inc_'+pl]['Value'] = self.planets[pl]['inc']
			par['inc_'+pl]['Prior_vals'] = [self.planets[pl]['inc'],0.01,0,90]
			par['e_'+pl]['Value'] = self.planets[pl]['ecc']
			par['e_'+pl]['Prior_vals'] = [self.planets[pl]['ecc'],0.01,0,1]
			par['w_'+pl]['Value'] = self.planets[pl]['w']
			par['w_'+pl]['Prior_vals'] = [self.planets[pl]['w'],1,-180,180]
			par['K_'+pl]['Value'] = self.planets[pl]['K']
			par['K_'+pl]['Prior_vals'] = [self.planets[pl]['K'],1,0,self.planets[pl]['K']+200]


		#par['ECs'] = []
		#par['LinCombs'] = []
		#par['FPs'] = []



#

#tab = GetTarget()
#targets = tab.getTarget('K2-261 b')

#cols = tab.fieldnames
#use = tab['pl_name'] == 'Kepler-11 c'
