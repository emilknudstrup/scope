#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 22:18:20 2022

@author: emil
"""
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)
from astropy.time import Time
import numpy as np
from astroplan.plots import plot_altitude, plot_sky
import astropy.units as u
from matplotlib.patches import Wedge, Circle
import matplotlib.ticker as mticker
from astroquery.simbad import Simbad
import datetime
import subprocess
import os
FONT = 12

from matplotlib.legend_handler import HandlerBase
class AnyObjectHandler(HandlerBase):
	def create_artists(self, legend, handles,
		x0, y0, width, height, fontsize, trans):
		l1 = plt.Line2D([x0,y0+width], [0.5*height,0.5*height],
			linestyle='-', color='k',lw=3.0,zorder=-1)
		l2 = plt.Line2D([x0,y0+width], [0.5*height,0.5*height],
			linestyle=handles[1], color=handles[0],lw=1.5,zorder=-1)          

		return [l1, l2] 

import matplotlib.legend_handler
from matplotlib.container import ErrorbarContainer
class re_order_errorbarHandler(matplotlib.legend_handler.HandlerErrorbar):
	"""
	Sub-class the standard error-bar handler 
	"""
	def create_artists(self, *args, **kwargs):
		#  call the parent class function
		a_list = matplotlib.legend_handler.HandlerErrorbar.create_artists(self, *args, **kwargs)
		# re-order the artist list, only the first artist is added to the
		# legend artist list, this is the one that corresponds to the markers
		a_list = a_list[-1:] + a_list[:-1]
		return a_list

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

def transit_plots(
				time,
				fix_target,
				observatory,
				path,
				target,
				full_transit=True,
				night=True,
				plot_moon=True,
				moon_phase=True,
				verbose=False
				):
	'''


	'''
	duration = target.duration*u.hour
	if full_transit:
		mid = time[0] + duration/2
		ing, eg = time[0].isot.split('T')[-1][:8], time[1].isot.split('T')[-1][:8]
		#day = (mid - 0.5).datetime
	else:
		mid = time
		ing, eg = (mid - duration/2).isot.split('T')[-1][:8], (mid + duration/2).isot.split('T')[-1][:8]
		day = (mid - 0.8).datetime
	ot = mid + np.linspace(-1*duration/2,duration/2,100)

	fig = plt.figure(figsize=(10,10))
	midT = mid.isot.split('T')[-1][:8]
	date = mid.isot.split('T')[0]
	#print(mid.isot)
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


	ax.axhline(y=30.,ls='--',color='C3',alpha=0.6)

	newxlabels = []
	xlabels = ax.get_xticklabels()
	xs = []
	for xlabel in xlabels:
		#text = r'$\rm {}$'.format(xlabel.get_text())
		text = r'$\rm {}$'.format(xlabel.get_text()[1:-1])
		xx, yy = xlabel.get_position()
		newxlabels.append(text)
		xs.append(xx)
	ticks_loc = ax.get_xticks().tolist()
	ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
	ax.set_xticklabels(newxlabels)
	#ax.set_xticks(xs)

	xl = kk.get_xlabel()

	xmin, xmax = ax.get_xlim()
	if plot_moon:
		mt = Time(np.linspace(xmin,xmax,200),format='plot_date')
		moon = observatory.moon_altaz(mt)
		ax.plot(mt.plot_date,moon.alt,'--',color='k')

	if night:
		#low = (observatory.sun_set_time(Time(day), which='next') - 0.5*u.hour).datetime
		low = (observatory.sun_set_time(Time(mid), which='previous') - 0.5*u.hour).datetime
		high = (observatory.sun_rise_time(Time(mid), which='next') + 0.5*u.hour).datetime
		ax.set_xlim(low,high)
		nf = observatory.sun_set_time(Time(mid), which='previous') - 0.5*u.hour
		xl = xl[:10] + nf.isot.split('T')[0] + xl[20:]
		# if low.hour < 5:
		# 	dd = int(nl[2].split('-')[-1])
		# 	new_string = nl[2][:-2]
		# 	if dd < 9: new_string = '{}0{:d}'.format(new_string,dd+1)
		# 	else: new_string = '{}{:d}'.format(new_string,dd+1)
		# 	nl[2] = new_string
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

	nl = [ll for ll in xl.split()]
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
	if verbose: print(path)
	plt.savefig(path,format='png')
	plt.close()

def visPlot(targets,observatory,time,moon=True,path=False,legend_outside=False,interact=True):

	site_name = observatory.name
	sun_set_next = observatory.sun_set_time(time, which='next')
	sun_set_previous = observatory.sun_set_time(time, which='previous')
	sun_rise_next = observatory.sun_rise_time(time, which='next')
	sun_rise_previous = observatory.sun_rise_time(time, which='previous')


	if sun_rise_next < sun_set_next:
		time = (sun_set_previous - 1*u.hour)
	else:
		time = (sun_set_next - 1*u.hour)
	ent = observatory.twilight_evening_nautical(Time(time), which='next')
	eat = observatory.twilight_evening_astronomical(Time(time), which='next')
	mnt = observatory.twilight_morning_nautical(Time(time), which='next')
	mat = observatory.twilight_morning_astronomical(Time(time), which='next')

	sun_rise = observatory.sun_rise_time(Time(time), which='next')
	sun_set = observatory.sun_set_time(Time(time), which='next')
	low = (sun_set - 1*u.hour).datetime
	high = (sun_rise + 1*u.hour).datetime

	time = (sun_set - 1*u.hour) + np.arange(0,25,0.125)*u.hour
	fig = plt.figure(figsize=(9,6))

	site_name = '\ '.join(site_name.split())
	title = r'$\rm Observed\ from\ {}.$'.format(site_name)

	# if len(targets) > 6:
	# 	legend_outside = True 
	#fig.suptitle(title)
	#ax = fig.add_subplot(121)


	# if legend_outside:
	# 	# fig = plt.figure(figsize=(12,8))

	# 	# site_name = observatory.name
	# 	# site_name = '\ '.join(site_name.split())		
	# 	# title = r'$\rm Observed\ from\ {}.$'.format(site_name)
	# 	# fig = plt.figure(figsize=(12,8))
	# 	# fig.suptitle(title,fontsize=FONT)
	# 	# ax = fig.add_subplot(111)

	# 	fig, axes = plt.subplots(ncols=1, nrows=3)
	# 	gs = axes[0].get_gridspec()
	# 	for ax in axes[:2]: ax.remove()

	# 	ax = fig.add_subplot(gs[:2])

	# 	axo = axes[-1]#plt.subplot(gs[2])
	# 	axo.spines['top'].set_color('none')
	# 	axo.spines['bottom'].set_color('none')
	# 	axo.spines['left'].set_color('none')
	# 	axo.spines['right'].set_color('none')
	# 	axo.patch.set_visible(False)
	# 	axo.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
	# 	fig.suptitle(title,y=1.008)
	# else:
	fig.suptitle(title,y=0.985,x=0.4)
	ax = fig.add_subplot(111)



	ax.text(ent.plot_date-0.02,90.3,'NAU\n{}'.format(ent.isot.split('T')[-1][:5]),fontsize=0.65*FONT)
	ax.text(eat.plot_date,90.3,'AST\n{}'.format(eat.isot.split('T')[-1][:5]),fontsize=0.65*FONT)
	ax.text(mnt.plot_date,90.3,'NAU\n{}'.format(mnt.isot.split('T')[-1][:5]),fontsize=0.65*FONT)
	ax.text(mat.plot_date-0.022,90.3,'AST\n{}'.format(mat.isot.split('T')[-1][:5]),fontsize=0.65*FONT)


	handles = []
	labels = []
	colors = ['C{}'.format(ii) for ii in range(0,10) if ii != 7]
	ncs = len(colors)
	lstyles = ['-','--','-.']
	lines = []
	brightness_shading = True
	#axsky = None
	for ii, obs_target in enumerate(targets):
		#per, T0, dur, ra, dec = target.makeTarget()
		name = obs_target.name
		#obs_target = FixedTarget(SkyCoord(ra=ra,dec=dec),name=name) 
		if not name: name = r'$\rm Target\ ' + str(ii) + '$'
		else: 
			name = name.replace(' ','\,')
			name = r'$\rm ' + name + '$'
		if ii > 0: brightness_shading = False
		color = colors[ii - ncs*(ii//ncs)]
		ls = lstyles[ii//ncs]
		handle = (color,ls)
		#title = r'$\rm {}\ observed\ from\ {}.$'.format(name,site_name)

		bstyle = {'linewidth' : 3.0, 'color' : 'k', 'zorder' : 4}
		plot_altitude(obs_target, observatory, time, ax=ax, style_kwargs=bstyle)
		style = {'linewidth' : 2.0, 'color' : color, 'zorder' : 8, 'linestyle' : ls}
		kk = plot_altitude(obs_target, observatory, time, ax=ax, style_kwargs=style,
		airmass_yaxis=False, min_altitude=7.0, brightness_shading=brightness_shading)
		lines.append((kk.get_lines()[ii*2],kk.get_lines()[ii*2+1],name))
		handles.append(handle)
		labels.append(name)

	if moon:
		mm = observatory.moon_altaz(time)
		mstyle = {'linewidth' : 1.5, 'color' : 'k', 'zorder' : 2,
		'linestyle' : '--', 'alpha' : 0.6}
		plot_altitude(mm, observatory, time, ax=ax, style_kwargs=mstyle)  

	ax.axhline(y=30.,linestyle='--',color=colors[3],alpha=0.6,zorder=1,linewidth=1.5)

	xl = kk.get_xlabel()
	nl = [ll for ll in xl.split()]
	xl = '$'
	for ll in nl: xl += '\mathrm{' + ll + '} \ '
	xl += '$'
	ax.set_xlabel(r'{}'.format(xl),fontsize=FONT)
	ax.set_xlim(low,high)
	newxlabels = []
	safe = low
	for nn in range(10):
		label = safe.isoformat().split('T')[1]
		label = int((int(label[:2]) + 2)%24)
		newxlabels.append(r'${:d}:00$'.format(label))
		safe += datetime.timedelta(hours = 2)
		if safe > high:
			break

	ax.set_xticklabels(newxlabels)

	ax.set_ylabel(r'$\rm Altitude \ (^\circ)$',fontsize=FONT)

	set_from_airmass = False
	if set_from_airmass:
		airmass_ticks = np.array([1.0, 1.05, 1.10, 1.25, 1.5, 1.75, 2.0, 3.0, 5.5])
		altitude_ticks = 90 - 180.*np.arccos(1/airmass_ticks)/np.pi
	else:
		altitude_ticks = np.linspace(30.,90.,10)
		altitude_ticks = np.append(altitude_ticks,19.47122063449069)
		altitude_ticks = np.append(altitude_ticks,9.594068226860458)
		airmass_ticks = 1./np.cos(np.pi*(90 - altitude_ticks)/180.)



	ax2 = ax.twinx()
	ax2.set_yticks(altitude_ticks)
	airmass_labels = [r'${:0.2f}$'.format(ll) for ll in airmass_ticks]
	ax2.set_yticklabels(airmass_labels)
	ax2.set_ylim(ax.get_ylim())
	ax2.set_ylabel(r'$\rm Airmass$')
	if legend_outside:
		#leg = ax.legend(handles, labels,
		leg = axo.legend(handles, labels,
			handler_map={tuple: AnyObjectHandler()},
			fancybox=True,shadow=True,
			fontsize=0.8*FONT,ncol=3,
			loc='upper left',bbox_to_anchor=(0.1, 0.7))
		if interact:
			lined = dict()
			for legline, origlines in zip(leg.get_lines(), lines):
				legline.set_picker(5)  # 5 pts tolerance
				lined[legline] = origlines

				def onpick(event):
					# on the pick event, find the orig line corresponding to the
					# legend proxy line, and toggle the visibility
					legline = event.artist
					bline, fline, text = lined[legline]
					for txt in ax.texts: txt.set_visible(False)
					alf = legline.get_alpha()
					if alf == 0.2:
						legline.set_alpha(1.0)

						plt.setp(bline,lw=3.0,zorder=25)
						plt.setp(fline,lw=2.0,zorder=26)
					else:
						legline.set_alpha(0.2)
						plt.setp(bline,lw=4.0,zorder=50)
						plt.setp(fline,lw=3.0,zorder=51)
						ax.text(0.9,0.9,text, 
						horizontalalignment='center',
						verticalalignment='center', 
						transform=ax.transAxes,
						bbox=dict(facecolor='none', alpha=1.0))

					fig.canvas.draw()

			fig.canvas.mpl_connect('pick_event', onpick)


	else:
		# plt.legend(handles, labels,
		# 		handler_map={tuple: AnyObjectHandler()},
		# 		fancybox=True,#shadow=True,
		# 		fontsize=0.8*FONT,ncol=3,
		# 		#loc='upper left',bbox_to_anchor=(0.15, 0.95))
		# 		loc='upper left',bbox_to_anchor=(0.1, 1.1))
		ncol = 1
		if len(targets)//12: ncol = len(targets)//12

		f = 0.8
		if len(targets) > 9: f = 0.7
		elif len(targets) > 15: f = 0.6
		elif len(targets) > 21: f = 0.5
		elif len(targets) > 27: f = 0.4

		xl = 1.28
		fa = 0.8
		if ncol == 2: 
			xl = 1.4
			fa = 1
		leg = ax.legend(handles, labels,loc="lower right",
					handler_map={tuple: AnyObjectHandler()},
					fancybox=True,shadow=True,
					fontsize=f*FONT,#framealpha=fa,
					bbox_to_anchor=(xl,0.15),
					ncol=ncol) 
		leg.set_zorder(100)
		if interact:
			lined = dict()
			for legline, origlines in zip(leg.get_lines(), lines):
				legline.set_picker(5)  # 5 pts tolerance
				lined[legline] = origlines

				def onpick(event):
					# on the pick event, find the orig line corresponding to the
					# legend proxy line, and toggle the visibility
					legline = event.artist
					bline, fline, text = lined[legline]
					for txt in ax.texts: txt.set_visible(False)
					alf = legline.get_alpha()
					if alf == 0.2:
						legline.set_alpha(1.0)

						plt.setp(bline,lw=3.0,zorder=25)
						plt.setp(fline,lw=2.0,zorder=26)
					else:
						legline.set_alpha(0.2)
						plt.setp(bline,lw=4.0,zorder=50)
						plt.setp(fline,lw=3.0,zorder=51)
						ax.text(0.9,0.9,text, 
						horizontalalignment='center',
						verticalalignment='center', 
						transform=ax.transAxes,zorder=30,
						bbox=dict(facecolor='white', alpha=1.0))

					fig.canvas.draw()

			fig.canvas.mpl_connect('pick_event', onpick)

	moon_phase = 1
	if moon_phase:
		axm = fig.add_subplot(9,9,1, autoscale_on=False, aspect='equal', xlim=[-0.05,1.05], ylim=[-0.05,1.05])
		mid = time[len(time)//2]
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

	plt.tight_layout()
	if path:
		pformat = path[-3:]
		plt.savefig(path,format=pformat,dpi=1000)
		plt.close()
	else:
		plt.show()

def showtime(observatory,target,times,moon=True,path=None,legend_outside=False,font=12):
	'''
	'''

	time = times[0]
	site_name = observatory.name
	fnames = []




	
	for ii, time in enumerate(times):

		fig = plt.figure()
		ax = fig.add_subplot(111)
		sun_set_next = observatory.sun_set_time(Time(time), which='next')
		sun_set_previous = observatory.sun_set_time(Time(time), which='previous')
		sun_rise_next = observatory.sun_rise_time(Time(time), which='next')
		sun_rise_previous = observatory.sun_rise_time(Time(time), which='previous')


		if sun_rise_next < sun_set_next:
			time = (sun_set_previous - 1*u.hour)
		else:
			time = (sun_set_next - 1*u.hour)
		ent = observatory.twilight_evening_nautical(Time(time), which='next')
		eat = observatory.twilight_evening_astronomical(Time(time), which='next')
		mnt = observatory.twilight_morning_nautical(Time(time), which='next')
		mat = observatory.twilight_morning_astronomical(Time(time), which='next')

		sun_rise = observatory.sun_rise_time(Time(time), which='next')
		sun_set = observatory.sun_set_time(Time(time), which='next')
		low = (sun_set - 1*u.hour).datetime
		high = (sun_rise + 1*u.hour).datetime	

		date = time.isot.split('T')[0]
		time = (sun_set - 1*u.hour) + np.arange(0,25,0.125)*u.hour

		site_name = '\, '.join(site_name.split())
		#site_name = ' '.join(site_name.split())
		title = r'$\rm Observed\ from\ {}.$'.format(site_name)
		plt.title(title)
		

		ax.text(ent.plot_date-0.025,85,'NAU\n{}'.format(ent.isot.split('T')[-1][:5]),fontsize=0.6*FONT)
		ax.text(eat.plot_date,85,'AST\n{}'.format(eat.isot.split('T')[-1][:5]),fontsize=0.6*FONT)
		ax.text(mnt.plot_date,85,'NAU\n{}'.format(mnt.isot.split('T')[-1][:5]),fontsize=0.6*FONT)
		ax.text(mat.plot_date-0.025,85,'AST\n{}'.format(mat.isot.split('T')[-1][:5]),fontsize=0.6*FONT)



		brightness_shading = True

		#name = target.ID
		#obs_target = FixedTarget(SkyCoord(ra=target.ra,dec=target.dec),name=name) 

		#title = r'$\rm {}\ observed\ from\ {}.$'.format(name,site_name)
		bstyle = {'linewidth' : 3.0, 'color' : 'k', 'zorder' : 4}
		plot_altitude(target, observatory, time, ax=ax, style_kwargs=bstyle)
		#plot_altitude(target, observatory, time, style_kwargs=bstyle)
		style = {'linewidth' : 2.0, 'color' : 'C0', 'zorder' : 8, 'linestyle' : '-'}
		#kk = plot_altitude(target, observatory, time, style_kwargs=style,
		kk = plot_altitude(target, observatory, time, ax=ax, style_kwargs=style,
						airmass_yaxis=False, min_altitude=7.0, 
						brightness_shading=brightness_shading)

		ax.axhline(y=30.,linestyle='--',color='C3',alpha=0.6,zorder=1,linewidth=1.5)
		if moon:
			mm = observatory.moon_altaz(time)
			mstyle = {
					'linewidth' : 1.5, 
					'color' : 'k', 
					'zorder' : 2,
					'linestyle' : '--', 
					'alpha' : 0.6
					}
			plot_altitude(mm, observatory, time, ax=ax, style_kwargs=mstyle)  

		xl = kk.get_xlabel()
		nl = [ll for ll in xl.split()]
		xl = '$'
		for ll in nl: xl += '\mathrm{' + ll + '} \ '
		xl += '$'
		ax.set_xlabel(r'{}'.format(xl),fontsize=font)
		
		newxlabels = []
		safe = low
		for nn in range(10):
			label = safe.isoformat().split('T')[1]
			label = int((int(label[:2]) + 2)%24)
			newxlabels.append(r'${:d}:00$'.format(label))
			safe += datetime.timedelta(hours = 2)
			if safe > high:
				break

		ax.set_xticklabels(newxlabels)

		ax.set_ylabel(r'$\rm Altitude \ (^\circ)$',fontsize=font)

		set_from_airmass = False
		if set_from_airmass:
			airmass_ticks = np.array([1.0, 1.05, 1.10, 1.25, 1.5, 1.75, 2.0, 3.0, 5.5])
			altitude_ticks = 90 - 180.*np.arccos(1/airmass_ticks)/np.pi
		else:
			altitude_ticks = np.linspace(30.,90.,10)
			altitude_ticks = np.append(altitude_ticks,19.47122063449069)
			altitude_ticks = np.append(altitude_ticks,9.594068226860458)
			airmass_ticks = 1./np.cos(np.pi*(90 - altitude_ticks)/180.)



		ax2 = ax.twinx()
		ax2.set_yticks(altitude_ticks)
		airmass_labels = [r'${:0.2f}$'.format(ll) for ll in airmass_ticks]
		ax2.set_yticklabels(airmass_labels)
		ax2.set_ylim(ax.get_ylim())
		ax2.set_ylabel(r'$\rm Airmass$')

		ax.set_xlim(low,high)

		pname = ''.join(path.split())
		try:
			os.mkdir(pname)
		except FileExistsError:
			print(pname+' exists.')
		fname = pname + '/' + pname.split('/')[-1] + '_{}.png'.format(date)
		plt.tight_layout()
		plt.savefig(fname,dpi=1000)
		#plt.gca()
		plt.close()


def schedulePlot(targets,site,time,
				fig = None,
				moon_phase = True,
				phase = False,
				exposures=[],
				**kwargs
				):

	if fig == None:
		site_name = site.name
		site_name = '\ '.join(site_name.split())		
		title = r'$\rm Observed\ from\ {}.$'.format(site_name)
		fig = plt.figure(figsize=(12,8))
		fig.suptitle(title,fontsize=FONT)
		ax1 = fig.add_subplot(211)
		if phase:
			ax2 = fig.add_subplot(413,sharex=ax1)
			ax3 = fig.add_subplot(414,sharex=ax1)
		else:
			ax3 = fig.add_subplot(212,sharex=ax1)

		if moon_phase:
			axm = fig.add_subplot(9,9,1, autoscale_on=False, aspect='equal', xlim=[-0.05,1.05], ylim=[-0.05,1.05])
		
	else:
		axes = fig.axes
		ax1, ax2, ax3 = axes[:3]
		if moon_phase:
			axm = axes[-1]
		
	colors = ['C{}'.format(ii) for ii in range(0,10) if ii != 7]
	ncs = len(colors)
	lstyles = ['-','--','-.']		
		
	
	prev_sun_set = site.sun_set_time(time, which='previous')
	next_sun_set = site.sun_set_time(time, which='next')

	next_sun_set_time = Time(str(next_sun_set.datetime),format='iso',scale='utc')
	prev_sun_set_time = Time(str(prev_sun_set.datetime),format='iso',scale='utc')

	tnext = abs(Time(time).jd-next_sun_set_time.jd)
	tprev = abs(Time(time).jd-prev_sun_set_time.jd)

	if tnext < tprev:
		sun_set = next_sun_set
	else:
		sun_set = prev_sun_set

	sun_set_time = Time(str(sun_set.datetime),format='iso',scale='utc')
	sun_rise_time = Time(str(sun_set.datetime),format='iso',scale='utc')
	sun_rise = site.sun_rise_time(sun_rise_time, which='next')

	time = Time(sun_set.jd + 0.5*(sun_rise.jd-sun_set.jd),format='jd',scale='utc')
	
	low = (sun_set - 1*u.hour).datetime
	high = (sun_rise + 1*u.hour).datetime
	handles = []
	labels = []
	brightness_shading = True
	
	for ii, obs_target in enumerate(targets):
		# if phase:
		# 	per = target.per
		# 	t0 = target.#df['T0 (BJD)'][idx]
		# 	ID = df['Planet'][idx]

		#Vmag = df['Vmag'][idx]
		#K = df['K (m/s)'][idx]
		#exp = int(df['Exp (s)'][idx])
		if not len(exposures):
			exp = 900
		else:
			exp = exposures[ii]
		#simbad = df['Name'][idx]	
		color = colors[ii - ncs*(ii//ncs)]
		linestyle = lstyles[ii//ncs]

		# try:
		# 	result_table = Simbad.query_object(simbad)			
		# 	ra, dec = result_table['RA'][0], result_table['DEC'][0]
		# 	cc = SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg))	
		# except RecursionError:
		# 	ra = df['R.A. (deg)'][idx]
		# 	dec = df['Dec (deg)'][idx]
		# 	cc = SkyCoord(ra=ra,dec=dec,unit=(u.deg, u.deg))	
		#obs_target = FixedTarget(SkyCoord(ra=cc.ra,dec=cc.dec),name=ID)
		
		if ii > 0: brightness_shading = False		
	
		#tw_ev, tw_mo = get_twilights(site,time,twilight='astronomical')
		
		tw_mo = site.twilight_morning_astronomical(Time(time), which='next')
		tw_ev = site.twilight_evening_astronomical(Time(time), which='previous')
		
		twi_mo_exp = tw_mo + exp*u.second/2
		twi_ev_exp = tw_ev - exp*u.second/2
		
		night_time = Time(np.linspace(twi_ev_exp.jd,twi_mo_exp.jd,400),format='jd',scale='utc')
		is_up = site.target_is_up(night_time,obs_target,horizon=15*u.deg)
		
		if any(is_up):
	
			
			night_altaz = site.altaz(night_time,obs_target)
			alt_const =  (15 < night_altaz.alt.value) & (85 > night_altaz.alt.value)
	
			if any(alt_const):
				night_time = night_time[alt_const]
				night_altaz = night_altaz[alt_const]
				
				max_alt = np.argmax(night_altaz.alt.value)
				
				st = night_time[max_alt] - exp*u.second/2
				et = night_time[max_alt] + exp*u.second/2
				if et > tw_mo:
					et = twi_mo_exp
					st = twi_mo_exp - exp*u.second
				elif st < tw_ev:
					st = twi_ev_exp
					et = twi_ev_exp + exp*u.second				
	
				obs_window = Time(np.linspace(st.jd,et.jd,50),scale='utc',format='jd')
	
	
			
				bstyle = {'linewidth' : 1.0, 'color' : 'k', 'zorder' : 4}
				plot_altitude(obs_target, site, time, ax=ax1, style_kwargs=bstyle, min_altitude=3.0,
							airmass_yaxis=False, brightness_shading=brightness_shading)
				
				obs_bkw = {'linewidth' : 3.0, 'color' : 'k', 'zorder' : 8, 'linestyle' : '-'}
				kk = plot_altitude(obs_target, site, obs_window, ax=ax1, style_kwargs=obs_bkw)
				
				
				
				obs_kw = {'linewidth' : 2.0, 'zorder' : 10, 'color' : color, 'linestyle' : linestyle}
				plot_altitude(obs_target, site, obs_window, ax=ax1, style_kwargs=obs_kw, min_altitude=3.)
				#lines.append((kk.get_lines()[ii*2],kk.get_lines()[ii*2+1],ID))
	
				handle = (color,linestyle)
				label = r'$\rm {}$'.format(obs_target.name)	
				handles.append(handle)
				labels.append(label)

				
				### ALTERNATIVELY USE LOW AND HIGH HERE
				x1, x2  = Time(sun_set - 1*u.hour,format='jd',scale='utc'),  Time(sun_rise + 1*u.hour,format='jd',scale='utc')
				xx = Time(np.linspace(x1.jd,x2.jd,200),format='jd',scale='utc')
				### Phase axis
				if phase:
					#ax2.set_ylim(-0.55,0.55)
					ax2.set_ylim(-0.05,1.05)
					xp = time2phase(xx.jd,per,t0) + 0.5
					ax2.plot(xx.plot_date,xp,'k-',lw=1.0)
					exp_idxs = (st.plot_date < xx.plot_date) & (xx.plot_date < et.plot_date)
					ax2.plot(xx.plot_date[exp_idxs],xp[exp_idxs],'-',color='k',lw=3.0,zorder=4)
					ax2.plot(xx.plot_date[exp_idxs],xp[exp_idxs],'-',color=color,lw=2.0,linestyle=linestyle,zorder=6)
	
				
				### Azimuth axis
				aaz = site.altaz(xx,obs_target)
				ax3.set_ylim(-5,365)
				phi = 180
				shift = True
				if shift:
					y_aaz = np.array([])
					for val in aaz.az.value: 
						
						if val > phi: val -= phi
						else: val += phi
						y_aaz = np.append(y_aaz,val)
				else:
					y_aaz = aaz.az.value
				ax3.plot(xx.plot_date,y_aaz,'k-',lw=1.0,zorder=-1)
				
				
				aaz_exp = site.altaz(obs_window,obs_target)
				if shift:
					y_aaz_exp = np.array([])
					for val in aaz_exp.az.value: 
						if val > phi: val -= phi
						else: val += phi
						y_aaz_exp = np.append(y_aaz_exp,val)
				else:
					y_aaz_exp = aaz_exp.az.value
				ax3.plot(obs_window.plot_date,y_aaz_exp,'k-',lw=3.0,zorder=6)
				ax3.plot(obs_window.plot_date,y_aaz_exp,'-',color=color,lw=2.0,linestyle=linestyle,zorder=6)		

			ax1.set_xlim(low,high)
			xmin, xmax = ax1.get_xlim()		
			
			#ax2.fill_between([xmin,xmax],0.2,0.3,color='C7',zorder=-2,alpha=0.4)
			#ax2.fill_between([xmin,xmax],-0.3,-0.2,color='C7',zorder=-2,alpha=0.4)
			if phase:
				ax2.fill_between([xmin,xmax],0.7,0.8,color='C7',zorder=-2,alpha=0.4)
				ax2.fill_between([xmin,xmax],0.2,0.3,color='C7',zorder=-2,alpha=0.4)
				ax2.set_ylabel(r'$\rm Phase$',fontsize=FONT*0.7)
	
				ax2.yaxis.set_major_locator(MultipleLocator(0.25))
				ax2.yaxis.set_minor_locator(MultipleLocator(0.125))
				plt.setp(ax2.get_xticklabels(),visible=False)

			#xlab = ax1.get_xlabel()
			#nl = [ll for ll in xlab.split()]
			#print(time,nl)
			#print(time.isot.split('T')[0],nl)
			#xl = '$'
			#for ll in nl: xl += '\mathrm{' + ll + '} \ '
			#xl += '$'
			ax1.set_xlabel('')
			#ax3.set_xlabel(r'${}$'.format(xl))
			ax3.set_xlabel(r'$\rm Time \ from \ {} \ [UTC]$'.format(time.isot.split('T')[0]))
			
			ax1.set_ylabel(r'$\rm Altitude \ (^\circ)$',fontsize=FONT*0.7)
			
			mt = Time(np.linspace(xmin,xmax,200),format='plot_date',scale='utc')
			moon = site.moon_altaz(mt)
	
			ax1.axhline(20.0,linestyle='--',color='C3',zorder=1)
			ax1.axhline(35.0,linestyle='--',color='C7',zorder=1)
			ax1.plot(mt.plot_date,moon.alt,'--',color='k',zorder=1)
	
			
			ax3cp = ax3.twinx()
			ax3cp.set_ylabel(r'$\rm Azimuth \ (^\circ)$',fontsize=FONT*0.7)
			
			ax1cp = ax1.twinx()
			altitude_ticks = np.linspace(30.,90.,6)
			altitude_ticks = np.append(altitude_ticks,19.47122063449069)
			altitude_ticks = np.append(altitude_ticks,9.594068226860458)
			airmass_ticks = 1./np.cos(np.pi*(90 - altitude_ticks)/180.)
			airmass_labels = [r'${:0.2f}$'.format(ll) for ll in airmass_ticks]
			ax1cp.set_yticklabels(airmass_labels)
			ax1cp.set_yticks(altitude_ticks)
			ax1cp.set_ylabel(r'$\rm Airmass$',fontsize=FONT*0.7)
			ax3.set_yticks([0,90,180,270,360])
			ax3cp.set_yticks([0,90,180,270,360])		
			if shift:
				ax3_ticklabels = [r'$\rm {:d}$'.format(lab) for lab in [180,270,0,90,180]]
				ax3cp.set_yticklabels(ax3_ticklabels)
				ax3.set_yticklabels([r'$\rm S$',r'$\rm W$',r'$\rm N$',r'$\rm E$',r'$\rm S$'])
			else:
				ax3.set_yticklabels([r'$\rm N$',r'$\rm E$',r'$\rm S$',r'$\rm W$',r'$\rm N$'])
			newxlabels = []
			safe = low
			for nn in range(10):
				label = safe.isoformat().split('T')[1]
				label = int(label[:2]) + 1
				newxlabels.append(r'${:d}:00$'.format(label))
				safe += datetime.timedelta(hours = 2)
				if safe > high: break
			ax3.set_xticklabels(newxlabels)	
			# if shift:
			# 	ax3.annotate(r'$\rm Windy$',xy=(-0.019,0.5),xycoords='axes fraction',
			# 				xytext=(-0.1,0.5),textcoords='axes fraction', rotation=0.0, fontsize=9, 
			# 				arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
			# 				bbox=dict(boxstyle="round", fc="white", ec="black"))
			# 	ax3.annotate('',xy=(-0.015,0.72),xycoords='axes fraction',
			# 				xytext=(-0.05,0.54),textcoords='axes fraction', rotation=0.0, fontsize=9, 
			# 				arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
			# 	ax3.annotate('',xy=(-0.015,0.28),xycoords='axes fraction',
			# 				xytext=(-0.05,0.49),textcoords='axes fraction', rotation=0.0, fontsize=9, 
			# 				arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
		
			plt.setp(ax1.get_xticklabels(),visible=False)
			plt.subplots_adjust(hspace=0.05)

			plt.legend(handles, labels,
					handler_map={tuple: AnyObjectHandler()},
					fancybox=True,shadow=True,
					fontsize=0.8*FONT,ncol=6,
					loc='upper left',bbox_to_anchor=(.06, 1.15))
		
			#midtime = night_time[int(len(night_time)/2)]
			if moon_phase:
				ph = site.moon_phase(time).value
				kk = site.moon_illumination(time)
				dd = int(time.isot.split('T')[0].split('-')[-1])
				angle = 90
				if dd > 14: angle += 180		
				
				dual_half_circle((0.5, 0.5), radius=0.5, angle=angle, ax=axm, phase=ph*180./np.pi,ill=kk)
				plt.setp(axm.get_xticklabels(),visible=False)
				plt.setp(axm.get_yticklabels(),visible=False)
				axm.set_xticks([])
				axm.set_yticks([])
	plt.show()

def skyPlot(targets,observatory,time,moon,path,warn_below_horizon=False):
	site_name = observatory.name
	sun_set_next = observatory.sun_set_time(Time(time), which='next')
	sun_set_previous = observatory.sun_set_time(Time(time), which='previous')
	sun_rise_next = observatory.sun_rise_time(Time(time), which='next')
	sun_rise_previous = observatory.sun_rise_time(Time(time), which='previous')

	if sun_rise_next < sun_set_next:
		time = (sun_set_previous - 1*u.hour)
	else:
		time = (sun_set_next - 1*u.hour)
	ent = observatory.twilight_evening_nautical(Time(time), which='next')
	eat = observatory.twilight_evening_astronomical(Time(time), which='next')
	mnt = observatory.twilight_morning_nautical(Time(time), which='next')
	mat = observatory.twilight_morning_astronomical(Time(time), which='next')

	sun_rise = observatory.sun_rise_time(Time(time), which='next')
	sun_set = observatory.sun_set_time(Time(time), which='next')
	low = (sun_set - 1*u.hour).datetime
	high = (sun_rise + 1*u.hour).datetime

	#time = (sun_set - 1*u.hour) + np.arange(0,25,0.125)*u.hour
	#time = (sun_set - 1*u.hour) + np.arange(0,25,0.5)*u.hour
	time = Time(np.linspace(sun_set.jd,sun_rise.jd,15),format='jd',scale='utc')#Time()
	#print(sun_rise.isot)
	#print(time.isot)


	site_name = '\ '.join(site_name.split())
	title = r'$\rm Observed\ from\ {}.\ Sunset\ UTC\ {}.$'.format(site_name,time[0].isot[:-4])
	fig = plt.figure()
	axsky = fig.add_subplot(projection='polar', facecolor="lightgoldenrodyellow")
	fig.suptitle(title,y=1.00)

	handles = []
	labels = []
	marks = []
	colors = ['C{}'.format(ii) for ii in range(0,10) if ii != 7]
	ncs = len(colors)
	mstyles = ['o','s','^']
	#for ii, target in enumerate(targets):
		#per, T0, dur, ra, dec = target.make_target()
		#name = target.ID
		#obs_target = FixedTarget(SkyCoord(ra=ra,dec=dec),name=name) 
	for ii, obs_target in enumerate(targets):
		name = obs_target.name
		if not name: name = r'$\rm Target\ ' + str(ii) + '$'
		else: 
			name = name.replace(' ','\,')
			name = r'$\rm ' + name + '$'
		color = colors[ii - ncs*(ii//ncs)]
		ms = mstyles[ii//ncs]
		style = {'marker' : ms,'facecolor' : color, 'edgecolor' : 'k', 'zorder' : 8, 'label' : name}

		kk = plot_sky(obs_target, observatory, time, ax=axsky,style_kwargs=style,warn_below_horizon=warn_below_horizon)
		marks.append(((kk.get_children()[ii],color,name)))
		nstyle = {'marker' : ms,'facecolor' : color, 'edgecolor' : 'k', 'zorder' : 9, 's' : 60, 'label' : None}
		plot_sky(obs_target, observatory, time[0], ax=axsky,style_kwargs=nstyle)


	my_handler_map = {ErrorbarContainer: re_order_errorbarHandler(numpoints=len(targets))}


	#angle = np.deg2rad(67.5)
	ncol = 1
	if len(targets)//12:
		ncol = len(targets)//12
	leg = axsky.legend(handler_map=my_handler_map,#loc="lower center",
				fancybox=True,#shadow=True,
				fontsize=0.8*FONT,
				#bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2),
				bbox_to_anchor=(1.37,1.0),
				ncol=ncol) 
	


	theta = np.linspace(0,2*np.pi,100)
	r = np.ones(len(theta))*70
	axsky.plot(theta,r,color='r')

	interact = 1
	if interact:
		lined = dict()
		#for legline, origlines in zip(leg.get_lines(), lines):
		for legline, origlines in zip(leg.legendHandles, marks):
			legline.set_picker(0.1)  # 5 pts tolerance
			#legline.set_pickradius(15)  # 5 pts tolerance
			lined[legline] = origlines

			def onpick(event):
				# on the pick event, find the orig line corresponding to the
				# legend proxy line, and toggle the visibility
				legline = event.artist
				bline, fline, text = lined[legline]
				#bline, text = lined[legline]
				for txt in axsky.texts: txt.set_visible(False)
				alf = legline.get_alpha()
				if alf == 0.2:
					legline.set_alpha(1.0)

					#plt.setp(bline,color='k',zorder=25)
					#plt.setp(bline,lw=3.0,zorder=25)
					#plt.setp(fline,lw=2.0,zorder=26)
					plt.setp(bline,facecolor=fline,zorder=25)
				else:
					legline.set_alpha(0.2)
					plt.setp(bline,color='k',zorder=50)
					#plt.setp(bline,lw=4.0,zorder=50)
					#plt.setp(fline,lw=3.0,zorder=51)
					axsky.text(0.9,0.9,text, 
					horizontalalignment='center',
					verticalalignment='center', 
					transform=axsky.transAxes,zorder=30,
					bbox=dict(facecolor='white', alpha=1.0))

				fig.canvas.draw()

		fig.canvas.mpl_connect('pick_event', onpick)


	plt.tight_layout()

	plt.show()


# def vis_plot(observatory,targets,time,moon,path,legend_outside=False,interact=True):
#   '''
#   '''
#   from matplotlib.legend_handler import HandlerBase
#   class AnyObjectHandler(HandlerBase):
#       def create_artists(self, legend, handles,
#                          x0, y0, width, height, fontsize, trans):
#           l1 = plt.Line2D([x0,y0+width], [0.5*height,0.5*height],
#                              linestyle='-', color='k',lw=3.0,zorder=-1)
#           l2 = plt.Line2D([x0,y0+width], [0.5*height,0.5*height],
#                              linestyle=handles[1], color=handles[0],lw=1.5,zorder=-1)          

#           return [l1, l2] 
  
#   site_name = observatory.name
#   sun_set_next = observatory.sun_set_time(Time(time), which='next')
#   sun_set_previous = observatory.sun_set_time(Time(time), which='previous')
#   sun_rise_next = observatory.sun_rise_time(Time(time), which='next')
#   sun_rise_previous = observatory.sun_rise_time(Time(time), which='previous')
  

#   if sun_rise_next < sun_set_next:
#     time = (sun_set_previous - 1*u.hour)
#   else:
#     time = (sun_set_next - 1*u.hour)
#   ent = observatory.twilight_evening_nautical(Time(time), which='next')
#   eat = observatory.twilight_evening_astronomical(Time(time), which='next')
#   mnt = observatory.twilight_morning_nautical(Time(time), which='next')
#   mat = observatory.twilight_morning_astronomical(Time(time), which='next')

#   sun_rise = observatory.sun_rise_time(Time(time), which='next')
#   sun_set = observatory.sun_set_time(Time(time), which='next')
#   low = (sun_set - 1*u.hour).datetime
#   high = (sun_rise + 1*u.hour).datetime
  
#   time = (sun_set - 1*u.hour) + np.arange(0,25,0.125)*u.hour
#   fig = plt.figure()
  
#   site_name = '\ '.join(site_name.split())
#   title = r'$\rm Observed\ from\ {}.$'.format(site_name)
  
#   if len(targets) > 6:
#     legend_outside = True 
#   #fig.suptitle(title)
#   #ax = fig.add_subplot(121)
  

#   if legend_outside:
#     fig, axes = plt.subplots(ncols=1, nrows=3)
#     gs = axes[0].get_gridspec()
#     for ax in axes[:2]: ax.remove()
	
#     ax = fig.add_subplot(gs[:2])

#     #from matplotlib.gridspec import GridSpec
#     #gs = GridSpec(3, 1)
#     #ax = plt.subplot(gs[:2])
	

#     axo = axes[-1]#plt.subplot(gs[2])
#     #axo = fig.add_subplot(122)
#     axo.spines['top'].set_color('none')
#     axo.spines['bottom'].set_color('none')
#     axo.spines['left'].set_color('none')
#     axo.spines['right'].set_color('none')
#     axo.patch.set_visible(False)
#     axo.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
#     fig.suptitle(title)
#   else:
#     fig.suptitle(title)
#     ax = fig.add_subplot(111)



#   ax.text(ent.plot_date-0.025,85,'NAU\n{}'.format(ent.isot.split('T')[-1][:5]),fontsize=0.6*FONT)
#   ax.text(eat.plot_date,85,'AST\n{}'.format(eat.isot.split('T')[-1][:5]),fontsize=0.6*FONT)
#   ax.text(mnt.plot_date,85,'NAU\n{}'.format(mnt.isot.split('T')[-1][:5]),fontsize=0.6*FONT)
#   ax.text(mat.plot_date-0.025,85,'AST\n{}'.format(mat.isot.split('T')[-1][:5]),fontsize=0.6*FONT)


#   handles = []
#   labels = []
#   colors = ['C{}'.format(ii) for ii in range(0,10) if ii != 7]
#   ncs = len(colors)
#   lstyles = ['-','--','-.']
#   lines = []
#   brightness_shading = True
#   #axsky = None
#   for ii, target in enumerate(targets):
#     per, T0, dur, ra, dec = target.make_target()
#     name = target.ID
#     obs_target = FixedTarget(SkyCoord(ra=ra,dec=dec),name=name) 
#     if not name: name = r'$\rm Target\ ' + str(ii) + '$'
#     else: 
#       name = name.replace(' ','\,')
#       name = r'$\rm ' + name + '$'
#     if ii > 0: brightness_shading = False
#     color = colors[ii - ncs*(ii//ncs)]
#     ls = lstyles[ii//ncs]
#     handle = (color,ls)
#     #title = r'$\rm {}\ observed\ from\ {}.$'.format(name,site_name)
	
#     bstyle = {'linewidth' : 3.0, 'color' : 'k', 'zorder' : 4}
#     plot_altitude(obs_target, observatory, time, ax=ax, style_kwargs=bstyle)
#     style = {'linewidth' : 2.0, 'color' : color, 'zorder' : 8, 'linestyle' : ls}
#     kk = plot_altitude(obs_target, observatory, time, ax=ax, style_kwargs=style,
#                  airmass_yaxis=False, min_altitude=7.0, brightness_shading=brightness_shading)
#     lines.append((kk.get_lines()[ii*2],kk.get_lines()[ii*2+1],name))
#     handles.append(handle)
#     labels.append(name)

#   if moon:
#     mm = observatory.moon_altaz(time)
#     mstyle = {'linewidth' : 1.5, 'color' : 'k', 'zorder' : 2,
#               'linestyle' : '--', 'alpha' : 0.6}
#     plot_altitude(mm, observatory, time, ax=ax, style_kwargs=mstyle)  
  
#   ax.axhline(y=30.,linestyle='--',color=colors[3],alpha=0.6,zorder=1,linewidth=1.5)

#   xl = kk.get_xlabel()
#   nl = [ll for ll in xl.split()]
#   xl = '$'
#   for ll in nl: xl += '\mathrm{' + ll + '} \ '
#   xl += '$'
#   ax.set_xlabel(r'{}'.format(xl),fontsize=font)
#   ax.set_xlim(low,high)
#   newxlabels = []
#   safe = low
#   for nn in range(10):
#     label = safe.isoformat().split('T')[1]
#     label = int((int(label[:2]) + 2)%24)
#     newxlabels.append(r'${:d}:00$'.format(label))
#     safe += datetime.timedelta(hours = 2)
#     if safe > high:
#       break
	
#   ax.set_xticklabels(newxlabels)

#   ax.set_ylabel(r'$\rm Altitude \ (^\circ)$',fontsize=font)
  
#   set_from_airmass = False
#   if set_from_airmass:
#     airmass_ticks = np.array([1.0, 1.05, 1.10, 1.25, 1.5, 1.75, 2.0, 3.0, 5.5])
#     altitude_ticks = 90 - 180.*np.arccos(1/airmass_ticks)/np.pi
#   else:
#     altitude_ticks = np.linspace(30.,90.,10)
#     altitude_ticks = np.append(altitude_ticks,19.47122063449069)
#     altitude_ticks = np.append(altitude_ticks,9.594068226860458)
#     airmass_ticks = 1./np.cos(np.pi*(90 - altitude_ticks)/180.)



#   ax2 = ax.twinx()
#   ax2.set_yticks(altitude_ticks)
#   airmass_labels = [r'${:0.2f}$'.format(ll) for ll in airmass_ticks]
#   ax2.set_yticklabels(airmass_labels)
#   ax2.set_ylim(ax.get_ylim())
#   ax2.set_ylabel(r'$\rm Airmass$')
#   if legend_outside:
#     leg = axo.legend(handles, labels,
#              handler_map={tuple: AnyObjectHandler()},
#              fancybox=True,shadow=True,
#              fontsize=0.8*FONT,ncol=3,
#              loc='upper left',bbox_to_anchor=(0.1, 0.7))
#     if interact:
#       lined = dict()
#       for legline, origlines in zip(leg.get_lines(), lines):
#         legline.set_picker(5)  # 5 pts tolerance
#         lined[legline] = origlines
	  
#       def onpick(event):
#         # on the pick event, find the orig line corresponding to the
#         # legend proxy line, and toggle the visibility
#         legline = event.artist
#         bline, fline, text = lined[legline]
#         for txt in ax.texts: txt.set_visible(False)
#         alf = legline.get_alpha()
#         if alf == 0.2:
#           legline.set_alpha(1.0)

#           plt.setp(bline,lw=3.0,zorder=25)
#           plt.setp(fline,lw=2.0,zorder=26)
#         else:
#           legline.set_alpha(0.2)
#           plt.setp(bline,lw=4.0,zorder=50)
#           plt.setp(fline,lw=3.0,zorder=51)
#           ax.text(0.9,0.9,text, 
#             horizontalalignment='center',
#             verticalalignment='center', 
#             transform=ax.transAxes,
#             bbox=dict(facecolor='white', alpha=1.0))
	   
#         fig.canvas.draw()

#       fig.canvas.mpl_connect('pick_event', onpick)
#       # import mpld3
#       # from mpld3 import plugins
#       # handles, labels = axo.get_legend_handles_labels() # return lines and labels
#       # interactive_legend = plugins.InteractiveLegendPlugin(
#       #   zip(handles,axo.collections),labels,
#       #   alpha_unsel=0.5,alpha_over=1.5,start_visible=True)     
#       # plugins.connect(fig, interactive_legend)
#       # mpld3.enable_notebook()

#   else:
#       plt.legend(handles, labels,
#            handler_map={tuple: AnyObjectHandler()},
#            fancybox=True,shadow=True,
#            fontsize=0.8*FONT,ncol=3,
#            loc='upper left',bbox_to_anchor=(.1, 1.1))
#   #plt.tight_layout()
#   #plt.legend(handles=hands, title='title', bbox_to_anchor=(1.05, 1), loc='upper left')#, prop=fontP)

#   if path:
#     pformat = path[-3:]
#     plt.savefig(path,format=pformat,dpi=1000)
#     plt.close()
#   else:
#     plt.show()