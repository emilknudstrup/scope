#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: emil
"""
import sys
import os
from scope.target import *
from scope.telescope import *
from scope.scope_target import *
from astropy.time import Time
import astropy.units as u

if __name__ == '__main__':
	import argparse 
	
	### Command line arguments to parse to script  
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
			description='Script to run scope target')

	

	function = parser.add_argument_group('Function')

	function.add_argument('function', type=str,
						help='Function to call',
						choices=['vis','sky','transit'])

	parser.add_argument('-t', '--targets', type=str, default='', 
						help='Target to plot')
	## input site
	parser.add_argument('-s', '--site', type=str, default='NOT', 
						help='Telescope/observing site')
	## input date in YYYY-MM-DD
	parser.add_argument('-d', '--date', type=str,default=Time.now().isot,
						help='Time for observation to be carried out must be formatetd |yyyy-mm-ddThh:mm:ss| UTC.')
	parser.add_argument('-moon', '--plot_moon', type=bool, default=True,
						help='Plot position of the Moon.')
	parser.add_argument('-sky', '--plot_sky', type=bool, default=False,
						help='Plot the position of the target on the sky.')

	## number of days from 'date' to search for transits
	parser.add_argument('-nd', '--ndays', type=int, default=180,
						help='Number of days from ´date´ to search for transits.')

	## periods for targets to search for transits
	parser.add_argument('-p', '--period', type=str, default=None,
						help='Update periods for transiting targets.')
	## T0s for targets to search for transits
	parser.add_argument('-t0', '--midtransit', type=str, default=None,
						help='Update T0s for transiting targets.')				
	## durations for targets to search for transits
	parser.add_argument('-dur', '--duration', type=str, default=None,
						help='Update T0s for transiting targets.')				

	## output file name
	parser.add_argument('-o', '--output', type=str, default=os.getcwd(),
						help='Path to store output')

	args, unknown = parser.parse_known_args(sys.argv[1:])

	target = GetTarget()
	observer = telescopeSites(args.site)


	inputTargets = args.targets.split(',')
	inputPeriods= args.period
	if inputPeriods:
		inputPeriods = inputPeriods.split(',')
	else:
		inputPeriods = []
	inputMidtransits = args.midtransit
	if inputMidtransits:
		inputMidtransits = inputMidtransits.split(',')
	else:
		inputMidtransits = []
	inputDurations = args.duration
	if inputDurations:
		inputDurations = inputDurations.split(',')
	else:
		inputDurations = []

	
	for ii, t in enumerate(inputTargets):
		target.byName(t)
		if len(inputPeriods):
			target[t]['pl_orbper'] = float(inputPeriods[ii])
		if len(inputMidtransits):
			target[t]['pl_tranmid'] = float(inputMidtransits[ii])
		if len(inputDurations):
			target[t]['pl_trandur'] = float(inputDurations[ii])


	target.createTargetlist()
	if args.function == 'vis':
		getVisPlot(target.targets, observer, args.date)
		if args.plot_sky:
			getSkyPlot(target.targets, observer, args.date, args.plot_moon)	
	elif args.function == 'sky':
		getSkyPlot(target.targets, observer, args.date, args.plot_moon)
	elif args.function == 'transit':
		start = Time(args.date)
		end = start + args.ndays*u.day
		start = start.isot
		end = end.isot
		target.createTransitTargetlist()
		telescope = Telescope(observer,20,observer.name)
		getTransits(target.targets, telescope, start, end,args.output)