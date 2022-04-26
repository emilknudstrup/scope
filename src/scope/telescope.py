#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Default telescopes
------------------

NOT
TNG
SONG
VLT

"""
from astroplan import Observer
import astropy.units as u
from astropy.coordinates import EarthLocation


class Telescope:
	'''Telescope class.
	
	Class to represent telescope parameters.
	
	Attributes
	----------
	location : class ``astroplan.Observer``
		Location of telescope.
	
	Vmag : float
		Limiting V magnitude.
	
	Methods
	-------
	print_sites()
	'''
	
	def __init__(self,location,Vmag=12.):
		'''Initialize
		Constructs all the necessary attributes for the telescope object.

		Parameters
		----------
			location : class ``astroplan.Observer``
				Location of telescope.
			Vmag : float
				Limiting V magnitude.
		'''
		self.location = location
		self.Vmag = Vmag
		

	def print_sites(self):
		'''Print sites.
		
		Prints all the default sites from ``astropy.coordinates.EarthLocation``.
		
		'''
		print(EarthLocation.get_site_names())

NOT = Telescope(Observer.at_site('Roque de los Muchachos'),10.5)
TNG = Telescope(Observer.at_site('Roque de los Muchachos'),11)
VLT = Telescope(Observer.at_site('Paranal Observatory'),12)
SONG = Telescope(Observer(longitude = -16.509, latitude = 28.2917,
						  elevation = 2395*u.m,timezone = 'Europe/London',
						  name = 'Teide Observatory'),8)
