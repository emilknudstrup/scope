#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 22:18:20 2022

@author: emil
"""
import matplotlib.pyplot as plt
from astropy.time import Time
import numpy as np
from astroplan.plots import plot_altitude, plot_sky
import astropy.units as u
from matplotlib.patches import Wedge, Circle

FONT = 12

def dual_half_circle(center, radius, angle=0, ax=None, colors=('lightyellow','k'),phase=0.0,ill=1.0,
					 **kwargs):
	"""
	Add two half circles to the axes *ax* (or the current axes) with the 
	specified facecolors *colors* rotated at *angle* (in degrees).
	"""
	if ax is None: ax = plt.gca()
	theta1, theta2 = angle, angle + 180
	C1, C2 = colors[0], colors[1]
	full, new = False, False
	if ill > 0.97: 
		C1, C2 = colors[0], colors[0]
		full = True
	elif ill < 0.07: 
		C1, C2 = colors[1], colors[1]
		new = True
	
	c = Circle(center,radius,fc='none',ec='k')
	ax.add_artist(c)
	w1 = Wedge(center, radius, theta1, theta2, fc=C1, ec='none', alpha=1.0, **kwargs)
	w2 = Wedge(center, radius, theta2, theta1, fc=C2, ec='none', alpha=1.0, **kwargs)
	for wedge in [w1, w2]: ax.add_artist(wedge)
	#if not (ill > 0.48) and (ill < 0.52):
	#print(ill)
	if not full or new:
		Cp = C1
		if ill < 0.5: 
			Cp = C2
			ill = 1 - ill
		elif ill < 0.515 and ill > 0.485: Cp = 'none'
	
	
		width, height = 1.0, 1.0*ill
		
		theta = np.deg2rad(np.arange(0.0, 360.0, 1.0))
		x = 0.5 * width * np.cos(theta)
		y = 0.5 * height * np.sin(theta)
	
		rtheta = np.radians(angle)
		R = np.array([
			[np.cos(rtheta), -np.sin(rtheta)],
			[np.sin(rtheta),  np.cos(rtheta)],
			])
		x, y = np.dot(R, np.array([x, y]))
		x += center[0]
		y += center[1]
	
		ax.fill(x, y,facecolor=Cp, edgecolor='none', linewidth=1, zorder=10)
	
	return [w1, w2, c]

def transit_plots(time,fix_target,observatory,path,target,
			   full_transit=True,night=True,plot_moon=True,moon_phase=True):
	'''


	'''
	duration = target.duration*u.hour

	if full_transit:
		mid = time[0] + duration/2
		ing, eg = time[0].isot.split('T')[-1][:8], time[1].isot.split('T')[-1][:8]
		day = (mid - 0.5).datetime
	else:
		mid = time
		ing, eg = (mid - duration/2).isot.split('T')[-1][:8], (mid + duration/2).isot.split('T')[-1][:8]
		day = (mid - 0.8).datetime
	ot = mid + np.linspace(-1*duration/2,duration/2,100)

	fig = plt.figure(figsize=(10,10))
	midT = mid.isot.split('T')[-1][:8]
	date = mid.isot.split('T')[0]
	obs_name = observatory.name
	top = '\mathrm{Transit} \ \mathrm{from} \ \ ' + ing + '\ \ \mathrm{to} \ \ ' + eg
	mid1 = '\mathrm{Midtransit} \ \mathrm{time} \ \ ' + midT
	#mid2 = date + '\ \mathrm{UTC}'
	mid2 = date.split('-')[-1] + '-' + date.split('-')[1] + '-' + date.split('-')[0] + '\ \mathrm{UTC}'
	onames = obs_name.split()
	obs_name = ' '
	for oname in onames: obs_name += '\mathrm{' + oname + '} \ '
	bottom = '\mathrm{Observed} \ \mathrm{from} \ ' + obs_name

	ax = fig.add_subplot(111)
	if moon_phase:
		axm = fig.add_subplot(9,9,1, autoscale_on=False, aspect='equal', xlim=[-0.05,1.05], ylim=[-0.05,1.05])

	ax.text(0.5,1.135,r'${}$'.format(top),
	ha='center',va='center',transform=ax.transAxes,fontsize=FONT)
	ax.text(0.5,1.10,r'${}, \ {}$'.format(mid1,mid2),
	ha='center',va='center',transform=ax.transAxes,fontsize=FONT)
	ax.text(0.5,1.06,r'${}$'.format(bottom),
	ha='center',va='center',transform=ax.transAxes,fontsize=FONT)

	transit = {'linewidth' : 3.5, 'zorder' : 10}
	plot_altitude(fix_target, observatory, ot, ax=ax, style_kwargs=transit)

	full_back = {'linewidth' : 5.0, 'color' : 'k', 'zorder' : 6}
	plot_altitude(fix_target, observatory, mid, ax=ax, style_kwargs=full_back)
	full = {'linewidth' : 3.5, 'color' : 'C1', 'zorder' : 8}
	kk = plot_altitude(fix_target, observatory, mid, ax=ax, style_kwargs=full,
	           brightness_shading=True, airmass_yaxis=False, min_altitude=3.0)


	ax.axhline(y=30.,linestyle='--',color='C3',alpha=0.6)

	newxlabels = []
	xlabels = ax.get_xticklabels()
	xs = []
	for xlabel in xlabels:
		text = r'$\rm {}$'.format(xlabel.get_text())
		xx, yy = xlabel.get_position()
		newxlabels.append(text)
		xs.append(xx)
	ax.set_xticklabels(newxlabels)
	ax.set_xticks(xs)

	xl = kk.get_xlabel()
	nl = [ll for ll in xl.split()]

	xmin, xmax = ax.get_xlim()
	if plot_moon:
		mt = Time(np.linspace(xmin,xmax,200),format='plot_date')
		moon = observatory.moon_altaz(mt)
		ax.plot(mt.plot_date,moon.alt,'--',color='k')

	if night:
		low = (observatory.sun_set_time(Time(day), which='next') - 0.5*u.hour).datetime
		high = (observatory.sun_rise_time(Time(low), which='next') + 0.5*u.hour).datetime
		ax.set_xlim(low,high)
		if low.hour < 5:
			dd = int(nl[2].split('-')[-1])
			new_string = nl[2][:-2]
			if dd < 9: new_string = '{}0{:d}'.format(new_string,dd+1)
			else: new_string = '{}{:d}'.format(new_string,dd+1)
			nl[2] = new_string
	else:
		xmin = Time(xmin,format='plot_date')
		xmax = Time(xmax,format='plot_date')

		srise = observatory.twilight_morning_astronomical(xmin, which='next')
		xlows = np.linspace(xmin.plot_date,srise.plot_date,10)
		sset = observatory.twilight_evening_astronomical(xmin, which='next')
		xhighs = np.linspace(sset.plot_date,xmax.plot_date,10)
		if xlows[-1] < xhighs[0]:
			ax.fill_between(xlows,0,y2=90,color='C7',alpha=0.5)
			ax.fill_between(xhighs,0,y2=90,color='C7',alpha=0.5)

	xl = '$'
	for ll in nl: xl += '\mathrm{' + ll + '} \ '
	xl += '$'
	ax.set_xlabel(r'{}'.format(xl),fontsize=FONT)

	ax.set_ylabel(r'$\rm Altitude \ (^\circ)$',fontsize=FONT)

	airmass_ticks = np.array([1.0, 1.25, 1.5, 1.75, 2.0, 3.0])
	altitude_ticks = 90 - np.degrees(np.arccos(1/airmass_ticks))

	ax2 = ax.twinx()
	ax2.set_yticks(altitude_ticks)
	airmass_labels = [r'${:0.2f}$'.format(ll) for ll in airmass_ticks]
	ax2.set_yticklabels(airmass_labels)
	ax2.set_ylim(ax.get_ylim())
	ax2.set_ylabel(r'$\rm Airmass$',fontsize=FONT)



	if moon_phase:
		ph = observatory.moon_phase(mid).value
		kk = observatory.moon_illumination(mid)
		dd = int(mid.isot.split('T')[0].split('-')[-1])
		angle = 90
		if dd > 14: angle += 180		
		
		dual_half_circle((0.5, 0.5), radius=0.5, angle=angle, ax=axm, phase=ph*180./np.pi,ill=kk)
		plt.setp(axm.get_xticklabels(),visible=False)
		plt.setp(axm.get_yticklabels(),visible=False)
		axm.set_xticks([])
		axm.set_yticks([])

	#ax.set_xlim(xmin, xmax)
	print(path)
	plt.savefig(path,format='png')
	plt.close()