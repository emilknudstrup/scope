#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 21:30:22 2022

@author: emil
"""
from astropy.coordinates import SkyCoord
from astroplan import (FixedTarget, AtNightConstraint, AltitudeConstraint, Observer,
					is_always_observable, is_observable, is_event_observable, months_observable,
					EclipsingSystem, MoonSeparationConstraint)
import astropy.units as u
from astropy.time import Time
import pandas as pd
import os
import numpy as np
from .plots import transit_plots, showtime,  visPlot

def get_twilights(site,time,twilight = 'astronomical'):
	if twilight == 'astronomical':
		twi_ev = site.twilight_evening_astronomical(Time(time), which='previous')
		twi_mo = site.twilight_morning_astronomical(Time(time), which='next')
	elif twilight == 'nautical':
		twi_ev = site.twilight_evening_nautical(Time(time), which='previous')
		twi_mo = site.twilight_morning_nautical(Time(time), which='next')
	elif twilight == 'civil':
		twi_ev = site.twilight_evening_civil(Time(time), which='previous')
		twi_mo = site.twilight_morning_civil(Time(time), which='next')
	return twi_ev, twi_mo

def transitTimes(observer,target,obs_time,ntransits=20,alt_const=30,moon_sep=30,full_transit=True,night=True): 
	'''
	Returns the times for the observations.

	Parameters:
		observer : class ``astroplan.Observer``
			Location of telescope.

		target : class ``target.Target``
			Target parameters.
		


	Returns:
		times : list of arrays
			Times for observations (start to finish).
		obs_target : class ``astroplan.FixedTarget``
			Target.

	'''
	constraints = [AltitudeConstraint(min=alt_const*u.deg),
				MoonSeparationConstraint(min=moon_sep*u.deg),
				AtNightConstraint.twilight_civil()]

	name = target.ID
	per, T0, dur, ra, dec = target.makeTarget()
	obs_target = FixedTarget(SkyCoord(ra=ra,dec=dec),name=name) 
	transiting_system = EclipsingSystem(primary_eclipse_time=T0, 
										orbital_period=per,duration=dur,name=name)
	if night:
		if full_transit:
			times = transiting_system.next_primary_ingress_egress_time(obs_time, n_eclipses=ntransits)
			obs = is_event_observable(constraints,observer,obs_target,times_ingress_egress=times)
		else:
			ing_egr = transiting_system.next_primary_ingress_egress_time(obs_time, n_eclipses=ntransits)
			times = transiting_system.next_primary_eclipse_time(obs_time, n_eclipses=ntransits)
			obs = is_event_observable(constraints,observer,obs_target,times=times)
		return times[obs[0][:]], obs_target
	else:
		times = transiting_system.next_primary_ingress_egress_time(obs_time, n_eclipses=ntransits)
	return times, obs_target


def getTransits(targets,telescope,start,end,path,plDict=None,
	alt_const=30,moon_sep=30,full_transit=True,night=True,
	limits={}):
	'''Pars

	'''

	if not path.endswith('/'): path += '/'

	if plDict: ndict = {}

	obs_time = Time(start)
	end_time = Time(end)
	ext = (end_time - obs_time).value
	for target in targets:
		if (target.Vmag < telescope.Vmag):
		#& (target.b > b_lim) & (target.vsini > vsini_lim):
			go = True
			for key in limits.keys():
				val = target.__dict__[key]
				if np.isnan(val): 
					go = False
					continue
				if (val <= limits[key][0]) or (val >= limits[key][1]):
					go = False
					continue
			if not go: continue
			
			try:
				nn = int(ext/target.per) + 1
			except OverflowError:
				print('Period for {} is {}.'.format(target.ID,target.per))
				print('Make sure you call scope.target.createTransitTargetlist().')
				continue
			times,fix_target = transitTimes(telescope.location,target,obs_time,
								ntransits=nn,alt_const=alt_const,moon_sep=moon_sep,
								full_transit=full_transit,night=night)
			
			plnames = target.ID.split()
			starname = target.starID
			dname = ''
			for plname in plnames: dname += plname
			# dname += '_'
			# starnames = starname.split()
			# for stname in starnames: dname += stname
			fullpath = path + dname	
			if len(times):		
				try:
					os.mkdir(fullpath)
				except FileExistsError:
					print(fullpath+' exists.')

			ntransits = len(times.value)
			if ntransits:
				if full_transit:
					print('Found {} transits for {} on:'.format(ntransits,target.ID))
				else:
					print('Found {} (partial) transits for {} on:'.format(ntransits,target.ID))
			
			atLeast1 = False
			for time in times.value:
				if full_transit:
					time = Time(time,scale='utc',format='jd')
					twi_ev = telescope.location.twilight_evening_civil(time[0], which='previous')
					fname = fullpath + '/' + fullpath.split('/')[-1] + '_' + twi_ev.isot.split('T')[0] + '.png'
					mid = Time(time[0].jd + (time[1].jd-time[0].jd)*0.5,scale='utc',format='jd')
					print('The night {}, midtransit at {}'.format(twi_ev.isot.split('T')[0],mid.isot))
				else:
					time = Time(time,scale='utc',format='iso')
					twi_ev = telescope.location.twilight_evening_civil(time, which='previous')
					fname = fullpath + '/' + fullpath.split('/')[-1] + '_' + twi_ev.isot.split('T')[0] + '.png'
					print('The night {}, midtransit at {}'.format(twi_ev.isot.split('T')[0],time.isot))
				try:
					transit_plots(time,fix_target,telescope.location,fname,target,full_transit=full_transit,night=night)
					atLeast1 = True
				except ValueError:
					pass
			
			if atLeast1: 
				if plDict:
					for key in plDict[target.ID].keys():
						try:
							ndict[key].append(plDict[target.ID][key])
						except KeyError:
							ndict[key] = [plDict[target.ID][key]]



	if plDict:
		newdf = pd.DataFrame(ndict)
		newdf.to_csv(path+'targets.csv',index=False)
		with open(path+'obs.log','w') as file:
			file.write('Telescope: {}.\n'.format(telescope.name))
			file.write('Observations from {} to {}.\n'.format(obs_time,end_time))
			file.write('Vmag < {}\n'.format(telescope.Vmag))
			for key in limits.keys():
				file.write('{} <= {} <= {}\n'.format(limits[key][0],key,limits[key][1]))	

def getVisPlot(targets,telescope,time,
	path=False,moon=True,
	time_scale='utc',time_format='isot',
	legend_outside=False,interact=True):

	time = Time(time,scale=time_scale,format=time_format)
	obs_targets = []
	for ii, target in enumerate(targets):
		per, T0, dur, ra, dec = target.makeTarget()
		name = target.ID
		obs_targets.append(FixedTarget(SkyCoord(ra=ra,dec=dec),name=name))
	visPlot(obs_targets,telescope,time,moon=moon,path=path,legend_outside=False,interact=True)


def checkVisibility(targets,telescope,start,end,path=None,alt_const=30,moon_sep=30,weeks=2,
	obslog=False,hours_from_astronomical=1,return_targets=False,plDict=None,keys=None):

	months = {
		'01' : 'January',
		'02' : 'February',
		'03' : 'March',
		'04' : 'April',
		'05' : 'May',
		'06' : 'June',
		'07' : 'July',
		'08' : 'August',
		'09' : 'September',
		'10' : 'October',
		'11' : 'November',
		'12' : 'December',
	}

	obs_time = Time(start)
	end_time = Time(end)	
	trange = [obs_time,end_time]

	constraints = [AltitudeConstraint(min=alt_const*u.deg),
				MoonSeparationConstraint(min=moon_sep*u.deg),
				AtNightConstraint.twilight_astronomical()]

	visible_targets = []

	if plDict: ndict = {}
	if path:
		if not path.endswith('/'): path += '/'

	if obslog: 
		line = 'Observations from the {}.\n'.format(telescope.name)
		line += 'From {} to {}.\n###\n'.format(obs_time.isot.split('T')[0],end_time.isot.split('T')[0])
	for target in targets:
		if (target.Vmag < telescope.Vmag):
			obs_target = FixedTarget(SkyCoord(ra=target.RA,dec=target.Dec,unit='deg'))
			vis = is_observable(constraints,telescope.location,obs_target,time_range=trange)[0]
			if vis:
				visible_targets.append(target.ID)
				print('{} is observable from {} in the interval {} to {}.'.format(target.starID,telescope.name,obs_time.isot.split('T')[0],end_time.isot.split('T')[0]))
				if obslog: line += '\n{} is observable.\nBest observed in:\n'.format(target.starID)
				new_time = obs_time
				days = [new_time]
				print('Best observed in:')
				obs_months = []
				while new_time < end_time:
					twi_ev = telescope.location.twilight_evening_astronomical(new_time, which='next') + hours_from_astronomical*u.hour
					twi_mo = telescope.location.twilight_morning_astronomical(new_time, which='next') - hours_from_astronomical*u.hour
					times = Time(np.arange(twi_ev.jd,twi_mo.jd,0.02),scale='utc',format='jd')
					vis = is_observable(constraints,telescope.location,obs_target,times=times)[0]
					if vis:
						if obslog: 
							month = months[twi_ev.isot.split('T')[0].split('-')[1]]
							if month not in obs_months: obs_months.append(month)
						print(twi_ev.isot.split('T')[0])
					days.append(new_time)
					new_time += weeks*7*u.day
				if obslog: 
					for month in months: line +=  month + ', '
					line += '\n###'
				if path:
					showtime(telescope.location,obs_target,days,path=path+target.starID)

				if plDict:
					if keys == None: keys = plDict[target.ID].keys()
					for key in keys:
						val = plDict[target.ID][key].data[0]
						#if not len(val): val = ' '
						try:
							ndict[key].append(val)
						except KeyError:
							ndict[key] = [val]

	if obslog:
		if not path:
			print('You need to specify a path to save the obslog.')
		else:
			with open(path+'obslog','w') as file:
				file.write(line)

	if plDict:
		newdf = pd.DataFrame(ndict)
		if not path:
			print('You need to specify a path to save the targets.')
		else:
			newdf.to_csv(path+'targets.csv',index=False)
	

	if return_targets: return visible_targets




				#vis = is_observable(constraints,telescope.location,obs_target,times=days)#[0]
				#print(vis,len(days))
				#best_months = months_observable(constraints,telescope.location,[obs_target])
				#print(best_months)
			# 		#print(new_time)
			# 		#trange = [new_time,new_time + weeks*7*u.day] 
			# 		#np = new_time + weeks*7*u.day
			# 		#print(new_time.isot,np.isot)
			# 		#twi_ev = telescope.location.twilight_evening_civil(new_time, which='next')
			# 		#twi_mo = telescope.location.twilight_morning_civil(new_time, which='next')
			# 		#times = Time(np.arange(twi_ev.jd,twi_mo.jd,0.02),scale='utc',format='jd')
			# 		#mid = twi_ev + (twi_mo-twi_ev)*0.5
			# 		#new_vis = is_observable(constraints,telescope.location,obs_target,times=mid,time_range=[twi_ev,twi_mo])[0]
			# 		new_vis = is_observable(constraints,telescope.location,obs_target,times=mid,time_range=trange)[0]
			# 		#new_vis = is_observable(constraints,telescope.location,obs_target,time_range=trange)[0]
			# 		#new_vis = is_observable(constraints,telescope.location,obs_target,times=np.array([twi_ev]))#[0]
			# 		#print(new_vis)
			# 		if new_vis:
			# 			print('{} is observable from {} in:'.format(target.starID,telescope.name))
			# 			print(new_time.isot.split('T')[0])

			#vis = is_observable(constraints,telescope.location,obs_target,time_range=trange)[0]
			# if vis:
			# 	new_time = obs_time# + weeks*7*u.day
			# 	print('{} is observable from {} in:'.format(target.starID,telescope.name))
			# 	best_months = months_observable(constraints,telescope.location,[obs_target],times=trange)
			# 	print(best_months)
			# # 	times = []
			# 	while new_time < end_time:
			# 		#print(new_time)
			# 		#trange = [new_time,new_time + weeks*7*u.day] 
			# 		#np = new_time + weeks*7*u.day
			# 		#print(new_time.isot,np.isot)
			# 		twi_ev = telescope.location.twilight_evening_civil(new_time, which='next')
			# 		twi_mo = telescope.location.twilight_morning_civil(new_time, which='next')

			# 		#new_vis = is_observable(constraints,telescope.location,obs_target,time_range=[twi_ev,twi_mo])[0]
			# 		new_vis = is_observable(constraints,telescope.location,obs_target,times=np.array([twi_ev]))#[0]
			# 		print(new_vis)
			# 		# if new_vis:
			# 		# 	print(new_time.isot.split('T')[0])
			# 		new_time += weeks*7*u.day

			# #print(vis)
