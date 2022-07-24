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
	
	def __init__(self,location,Vmag=12.,name='Name'):
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
		self.name = name
		

	def print_sites(self):
		'''Print sites.
		
		Prints all the default sites from ``astropy.coordinates.EarthLocation``.
		
		'''
		print(EarthLocation.get_site_names())

NOT = Telescope(Observer.at_site('Roque de los Muchachos'),10.5,'NOT')
TNG = Telescope(Observer.at_site('Roque de los Muchachos'),11,'TNG')
VLT = Telescope(Observer.at_site('Paranal Observatory'),12,'VLT')
SONG = Telescope(Observer(longitude = -16.509, latitude = 28.2917,
						  elevation = 2395*u.m,timezone = 'Europe/London',
						  name = 'Teide Observatory'),8,'SONG')
SONG_AUS = Telescope(Observer(longitude = 151.8554, latitude = -27.7977,
						  elevation = 682*u.m,timezone = 'Australia/Queensland',
						  name = 'Mount Kent Observatory'),5,'SONG_AUS')
FUT = Telescope(Observer(longitude = 151.8554, latitude = -27.7977,
						  elevation = 682*u.m,timezone = 'Australia/Queensland',
						  name = 'Mount Kent Observatory'),15,'FUT')
def telescopeSites(telescope):
	telescopes = {
		'SONG' : Observer(longitude = -16.509,#343.5033*u.deg,
			latitude = 28.2917*u.deg,
			elevation = 2395*u.m,
			timezone = 'Europe/London',
			name = 'Teide Observatory'), 
		'MUSCAT2' : Observer(longitude = -16.509,#343.5033*u.deg,
			latitude = 28.2917*u.deg,
			elevation = 2395*u.m,
			timezone = 'Europe/London',
			name = 'Teide Observatory'),
		'ORM' : Observer.at_site('Roque de los Muchachos'),
		'NOT' : Observer.at_site('Roque de los Muchachos'),
		'TNG' : Observer.at_site('Roque de los Muchachos'),
		'HARPS-N' : Observer.at_site('Roque de los Muchachos'),
		'DCT' : Observer.at_site('Discovery Channel Telescope'),
		'VLT' : Observer.at_site('Paranal Observatory'),
		'Paranal' : Observer.at_site('Paranal Observatory'),
		'LaSilla' : Observer.at_site('La Silla Observatory'),
		'HARPS' : Observer.at_site('La Silla Observatory'),
		'Keck' : Observer.at_site('Keck Observatory'),
		'Subaru' : Observer.at_site('Subaru'),
		'OHP' : Observer(longitude = 43.9308*u.deg,
			latitude = 5.7133*u.deg,
			elevation = 650*u.m,
			timezone = 'Europe/Paris',
			name = 'Haute-Provence Observatory'),
		'CASLEO' : Observer(longitude = -31.7986*u.deg,
			latitude = -69.2956*u.deg,
			elevation = 2483*u.m,
			timezone = 'Etc/GMT-3',
			name = 'Leoncito Astronomical Complex'),
		'BTA' : Observer(longitude = 43.6468*u.deg, 
			latitude = 41.4405*u.deg,
			elevation = 2070*u.m,
			timezone = 'Etc/GMT+5',
			name = 'Special Astrophysical Observatory'),
		'TLS' : Observer(longitude = 11.711167*u.deg,
			latitude = 50.98011*u.deg,
			elevation = 341*u.m,
			timezone = 'Europe/Berlin',
			name = 'Thueringer Landessternwarte Tautenbug')
	}
	return telescopes[telescope]