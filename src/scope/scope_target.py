#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 21:30:22 2022

@author: emil
"""
from astropy.coordinates import SkyCoord
from astroplan import (FixedTarget, AtNightConstraint, AltitudeConstraint, Observer,
					is_event_observable, EclipsingSystem, MoonSeparationConstraint)
from .plots import transit_plots
import astropy.units as u
from astropy.time import Time
import os

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
			times = transiting_system.next_primary_eclipse_time(obs_time, n_eclipses=ntransits)
			obs = is_event_observable(constraints,observer,obs_target,times=times)
		return times[obs[0][:]], obs_target
	else:
		times = transiting_system.next_primary_ingress_egress_time(obs_time, n_eclipses=ntransits)
	return times, obs_target


def getTransits(targets,telescope,path,start,end,bLim=0.0,vsinilim=0.0,alt_const=30,moon_sep=30,full_transit=True,night=True):
	'''Pars

	'''

	if not path.endswith('/'): path += '/'

	for target in targets:
		if (target.Vmag < telescope.Vmag) & (target.b > bLim) & (target.vsini > vsinilim):
			obs_time = Time(start)
			end_time = Time(end)
			ext = (end_time - obs_time).value
			nn = int(ext/target.per) + 1
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
			try:
				os.mkdir(fullpath)
			except FileExistsError:
				print(fullpath+' exists.')
			for time in times.value:
				time = Time(time,scale='utc',format='jd')
				fname = fullpath + '/' + fullpath.split('/')[-1] + '_' + time[0].isot.split('T')[0] + '.png'
				try:
					print(fname, 'asdasd')
					transit_plots(time,fix_target,telescope.location,fname,target,full_transit=full_transit,night=night)
				except ValueError:
					print(fname,time)




