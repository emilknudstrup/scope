#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: emil
"""
import sys
from scope.target import *
from scope.telescope import *
from scope.scope_target import *
from astropy.time import Time

if __name__ == '__main__':
	import argparse 
	
	### Command line arguments to parse to script  
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
			description='Script to run scope target')

	

	function = parser.add_argument_group('Function')

	function.add_argument('function', type=str,
						help='Function to call',
						choices=['vis','sky'])

	parser.add_argument('-t', '--targets', type=str, default='', 
						help='Target to plot')
	## input site
	parser.add_argument('-s', '--site', type=str, default='NOT', 
						help='Telescope/ observing site')
	## input date in YYYY-MM-DD
	parser.add_argument('-d', '--date', type=str,default=Time.now().isot,
						help='Time for observation to be carried out must be formatetd |yyyy-mm-ddThh:mm:ss| UTC.')
	parser.add_argument('-moon', '--plot_moon', type=bool, default=True,
						help='Plot position of the Moon.')
	parser.add_argument('-sky', '--plot_sky', type=bool, default=False,
						help='Plot the position of the target on the sky.')


	args, unknown = parser.parse_known_args(sys.argv[1:])

	target = GetTarget()
	telescope = telescopeSites(args.site)


	inputTargets = args.targets.split(',')
	for t in inputTargets:
		target.byName(t)

	target.createTargetlist()
	if args.function == 'vis':
		getVisPlot(target.targets, telescope, args.date)
		if args.plot_sky:
			getSkyPlot(target.targets, telescope, args.date, args.plot_moon)	
	elif args.function == 'sky':
		getSkyPlot(target.targets, telescope, args.date, args.plot_moon)