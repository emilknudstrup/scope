# scope
Finding transits of exoplanets observable from given locations


### Installation
`cd /path/to`

`git clone https://github.com/emilknudstrup/scope.git`

`cd /path/to/scope`

`python -m pip install .`

`pip install -r requirements.txt`

### Example: checking for transits of all planets at given site in given time interval
```python
import pandas as pd
import os

import scope
## Initialize target
tar = scope.target.GetTarget()
## Get, prepare, and list transit targets
tar.getRMTargets()
tar.createTransitTargetlist()

## Choose a telescope/site
tel = scope.telescope.VLT

## Create path to store output files
path = os.getcwd()
path += '/vis'
try:
	os.mkdir(path)
except FileExistsError:
	pass

## Start and end dates
start = '2023-04-01 12:00'
end = '2023-10-01 12:00'

## Search for transits
plim = 0 #1 for a range in period
if plim:
	scope.scope_target.getTransits(tar.targets,tel,start,end,path,plDict=tar.plDict,limits={'per' : [8,1e3]})
else:
	scope.scope_target.getTransits(tar.targets,tel,start,end,path,plDict=tar.plDict)

## Store output with some selected, key parameters 
tardf = pd.read_csv(path+'/targets.csv')
keys = ['pl_name','pl_orbper','pl_orbeccen','pl_imppar','st_vsin']

ndict = {}
for key in keys: ndict[key] = list(tardf[key])
newdf = pd.DataFrame(ndict)
newdf.to_csv(path+'/keypars_targets.csv',index=False)

```

### Example: visibility and/or sky plot
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

### Example: search for transits for a single target
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

### Example: create target list from a CSV file

scope can generate a target list and plots from a CSV file if it's in the following format:

| name     | host          | RA           | Dec           |
| -------- | --------------| ------------ | ------------- |
| TOI-2025 | TYC 4595-797-1| 18 51 10.839 | 82 14 43.5636 |
| ...      | ...           | ...          | ...           |
| TOI-2158 | HD 348661     | 18 27 14.461 | 82 14 43.5636 |

This is how you would load it in with scope
```python
import scope

## Observing date
date = '2022-08-22T12:00:00'

## Choose a telescope
tel = scope.telescope.NOT

## Initialize target
tar = scope.target.GetTarget()
## Read in a target list from a CSV file
tar.targetFromcsv('spreadsheet.csv',skiprows=0)#The skiprows is to remove a "header"
## Create target lsit
tar.createTargetlist()

## Visibility and sky plot
scope.scope_target.getVisPlot(tar.targets,tel.location,date)
scope.scope_target.getSkyPlot(tar.targets,tel.location,date)

```


