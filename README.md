# scope
Finding transits of exoplanets observable from given locations


### Installation
`cd /path/to`

`git clone https://github.com/emilknudstrup/scope.git`

`cd /path/to/scope`

`python -m pip install .`

`pip install -r requirements.txt`


### Example for visibility and/or sky plot
```python
import scope

## Initialize target
tar = scope.target.GetTarget()
## Search by name (NASA Exoplanet Archive)
tar.byName('K2-290')
## ...or from SIMBAD
sim = 0
if sim:
	tar.fromSIMBAD('TYC 6193-663-1')
## ...or manually
man = 0
if man:
	ra, dec = '15 39 25.86253', '-20 11 55.77049'
	tar.byHand(ra,dec,name='K2-290')

## Create target list
tar.createTargetlist()

## Choose a telescope
tel = scope.telescope.VLT

## Visibility and sky plot
scope.scope_target.getVisPlot(tar.targets,tel.location,'2022-07-23T08:00:00')
scope.scope_target.getSkyPlot(tar.targets,tel.location,'2022-07-23T08:00:00')

```

### Example to search for transit for a single target
```python
import os
import scope

## Initialize the telescope 
tel = scope.telescope.TNG
#tel.Vmag = 12 ## Increase/decrease Vmag limit. Default TNG is 11

## Initialize the target and query NASA Exoplanet Archive for the planet
tar = scope.target.GetTarget()
tar.byName('HD 118203 b')
## This is only needed if we know (/think) we have better transit parameters than the default from NASA Exoplanet Archive
## or if we see (printed) some (essential) parameters are missing
betterPars = 1
if betterPars:
	tar.plDict['HD 118203 b']['pl_orbper'] = 6.134985
	tar.plDict['HD 118203 b']['pl_orbincl'] = 88.88
	tar.plDict['HD 118203 b']['pl_tranmid'] = 2458712.66147
	tar.plDict['HD 118203 b']['pl_radj'] = 1.136
	tar.plDict['HD 118203 b']['pl_ratror'] = 0.05552
	tar.plDict['HD 118203 b']['pl_trandur'] = 0.23543*24
	tar.plDict['HD 118203 b']['pl_imppar'] = 0.1111

## Create the transit target list
tar.createTransitTargetlist()

## Create path to output
path = os.getcwd()
path += '/hd118203'
try:
	os.mkdir(path)
except FileExistsError:
	print('{} exists.'.format(path))

## Start and end dates
start = '2023-04-01 12:00'
end = '2023-10-01 12:00'

## Search for transits
scope.scope_target.getTransits(tar.targets,tel,start,end,path)	

```