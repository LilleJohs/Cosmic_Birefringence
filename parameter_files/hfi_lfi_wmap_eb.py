# The ordering of the frequency bands is a bit weird, but the reason is explained in the
# params variable in run_cosmic_birefringence_analysis.py

maps_param = {
	'030': {
		'data_split': ['A', 'B'],
		'map_root': '/path/to/map/npipe6v20{}/npipe6v20{}_030_map_K.fits',
		'n_side': 1024,
		'data_set': 'npipe'
	},
	'044': {
		'data_split': ['A', 'B'],
		'map_root': '/path/to/map/npipe6v20{}/npipe6v20{}_044_map_K.fits',
		'n_side': 1024,
		'data_set': 'npipe'
	},
	'K': {
		# The code was initially coded to only support channels with time or detectors splits.
		# K and Ka only have one map each, so I treat the as two detector split maps. They each have their own
		# miscalibration angle. The only problem is doing a frequency dependence search for beta. But
		# I have a solution for that in the parameter file used for cb.py
		'data_split': ['1', 'a1'],
		'map_root': '/path/to/map/wmap_iqumap_r9_9yr_K{}_v5.fits',
		'n_side': 512,
		'data_set': 'wmap',
	},
	'Q': {
		'data_split': ['1', '2'],
		'map_root': '/path/to/map/wmap_iqumap_r9_9yr_Q{}_v5.fits',
		'n_side': 512,
		'data_set': 'wmap',
	},
	'V': {
		'data_split': ['1', '2'],
		'map_root': '/path/to/map/wmap_iqumap_r9_9yr_V{}_v5.fits',
		'n_side': 512,
		'data_set': 'wmap'
	},
	'070': {
		'data_split': ['A', 'B'],
		'map_root': '/path/to/map/npipe6v20{}/npipe6v20{}_070_map_K.fits',
		'n_side': 1024,
		'data_set': 'npipe'
	},
	'W': {
		# W1 and W2 are treated as data split for W
		# and W3 and W4 are treated as data split for W_
		'data_split': ['1', '2'],
		'map_root': '/path/to/map/wmap_iqumap_r9_9yr_W{}_v5.fits',
		'n_side': 512,
		'data_set': 'wmap'
	},
	'W_': {
		'data_split': ['3', '4'],
		'map_root': '/path/to/map/wmap_iqumap_r9_9yr_W{}_v5.fits',
		'n_side': 512,
		'data_set': 'wmap'
	},
	'100': {
		'data_split': ['A', 'B'],
		'map_root': '/path/to/map/npipe6v20{}/npipe6v20{}_100_map_K.fits',
		'n_side': 2048,
		'data_set': 'npipe'
	},
	'143': {
		'data_split': ['A', 'B'],
		'map_root': '/path/to/map/npipe6v20{}/npipe6v20{}_143_map_K.fits',
		'n_side': 2048,
		'data_set': 'npipe'
	},
	'217': {
		'data_split': ['A', 'B'],
		'map_root': '/path/to/map/npipe6v20{}/npipe6v20{}_217_map_K.fits',
		'n_side': 2048,
		'data_set': 'npipe'
	},
	'353': {
		'data_split': ['A', 'B'],
		'map_root': '/path/to/map/npipe6v20{}/npipe6v20{}_353_map_K.fits',
		'n_side': 2048,
		'data_set': 'npipe'
	},
}