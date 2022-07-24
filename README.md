# scope
Finding transits of exoplanets observable from given locations


### Installation
`cd /path/to`

`git clone https://github.com/emilknudstrup/scope.git`

`cd /path/to/scope`

`python -m pip install .`

`pip install -r requirements.txt`


### Example
```python
import scope

## Initialize target
tar = scope.target.GetTarget()
## Search by name
tar.byName('K2-290')
## Create target lsit
tar.createTargetlist()

## Choose a telescope
tel = scope.telescope.VLT

## Visibility and sky plot
scope.scope_target.getVisPlot(tar.targets,tel.location,'2022-07-23T08:00:00')
scope.scope_target.getSkyPlot(tar.targets,tel.location,'2022-07-23T08:00:00')

```