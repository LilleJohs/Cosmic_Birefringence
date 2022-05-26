import numpy as np
import healpy as hp
import camb
from astropy.io import fits
import tools_fast as tools

import sys
sys.path.append('parameter_files/')
from hfi_lfi_wmap_eb import maps_param

# This file gets the theoretical LCDM parameters from CAMB, then beam smooth them
# and stores them in a fileformat that cb.py can read.

l_max = 1600

#Pixel window functions
w_pix = {
	'2048': hp.sphtfunc.pixwin(2048, pol=True, lmax=l_max+1),
	'1024': hp.sphtfunc.pixwin(1024, pol=True, lmax=l_max+1),
	'512': hp.sphtfunc.pixwin(512, pol=True, lmax=l_max+1)
}

#NPIPE
beam_function_folder = 'path/to/beams/window_functions/npipe_aux/beam_window_functions/AB/'
wl_root_name = 'Wl_npipe6v20_{}{}x{}{}.fits'
bl_root_name = 'Bl_TEB_npipe6v20_{}{}x{}{}.fits'

#WMAP
beam_function_folder_wmap = 'path/to/beams/wmap_ampl_bl_{}{}_9yr_v5p1.txt'

# LCDM spectra. Cosmological parameters from Planck 2018
cp = camb.set_params(tau=0.0544, ns=0.9649, H0=67.36, ombh2=0.02237, omch2=0.12, As=2.1e-09, lmax=l_max)
camb_results = camb.get_results(cp)
all_cls_th = camb_results.get_cmb_power_spectra(lmax=l_max, raw_cl=True)['total']

# Number of bands
nob = len(maps_param.keys())

# Number of map combinations (excluding auto spectra)
size = tools.size(nob)
print(size)

# First index is multipole
# Second and third indcies are frequency bands
# Fourth and fifth indices for the time/detector split for the given band.
# Sixth index is what polarization spectra. 0 = EE, 1 = BB
c_l_th = np.zeros((l_max+1, nob, nob, 2, 2, 2))	

# Loop through all combinations of maps
i = 0
for freq1, freq1_param in maps_param.items():
	print(i)
	j = 0
	for freq2, freq2_param in maps_param.items():
		for im in range(2):
			for jm in range(2):
				if freq1_param['data_set'] == 'npipe' and freq2_param['data_set'] == 'npipe':
					# If both maps are NPIPE, there is a window function for this combination available
					wl_data = wl_data = fits.open(beam_function_folder + wl_root_name.format(freq1, freq1_param['data_split'][im], freq2, freq2_param['data_split'][jm]))
					
					for l in range(l_max+1):
						c_l_th[l, i, j, im, jm, 0] = all_cls_th[l, 1] * w_pix[str(freq1_param['n_side'])][1][l] * w_pix[str(freq2_param['n_side'])][1][l] * wl_data['EE'].data['EE_2_EE'][0, l] * 2.7255**2 #EE
						c_l_th[l, i, j, im, jm, 1] = all_cls_th[l, 2] * w_pix[str(freq1_param['n_side'])][1][l] * w_pix[str(freq2_param['n_side'])][1][l] * wl_data['BB'].data['BB_2_BB'][0, l] * 2.7255**2 #BB
				
				
				else:
					# At least one map is not NPIPE. So we use beam transfer functions instead of window functions
					if freq1_param['data_set'] == 'npipe':
						data_1 = fits.open(beam_function_folder + bl_root_name.format(freq1.replace("_", ""), freq1_param['data_split'][im], freq1, freq1_param['data_split'][im]))
						bl_data_1 = np.array([data_1[1].data['E'][:], data_1[1].data['B'][:]])
					elif freq1_param['data_set'] == 'wmap':
						data_1 = np.loadtxt(beam_function_folder_wmap.format(freq1.replace("_", ""), freq1_param['data_split'][im]))
						
						if len(data_1[:, 1]) < l_max:
							bl_data_1 = np.zeros((2, l_max+1))
							bl_data_1[:, :len(data_1[:, 1])] = np.array([data_1[:, 1], data_1[:, 1]])
						else:
							bl_data_1 = np.array([data_1[:, 1], data_1[:, 1]])

					if freq2_param['data_set'] == 'npipe':
						data_2 = fits.open(beam_function_folder + bl_root_name.format(freq2.replace("_", ""), freq2_param['data_split'][jm], freq2, freq2_param['data_split'][jm]))
						bl_data_2 = np.array([data_2[1].data['E'][:], data_2[1].data['B'][:]])
					elif freq2_param['data_set'] == 'wmap':
						data_2 = np.loadtxt(beam_function_folder_wmap.format(freq2.replace("_", ""), freq2_param['data_split'][jm]))
						
						if len(data_2[:, 1]) < l_max:
							# The beam transfer function for this WMAP map does not go all the way up to lmax
							# so we set the beam transfer functions to zero for the multipoles that dont have one
							bl_data_2 = np.zeros((2, l_max+1))
							bl_data_2[:, :len(data_2[:, 1])] = np.array([data_2[:, 1], data_2[:, 1]])
						else:
							bl_data_2 = np.array([data_2[:, 1], data_2[:, 1]])

					for l in range(l_max+1):
						#EE
						c_l_th[l, i, j, im, jm, 0] = all_cls_th[l, 1] * w_pix[str(freq1_param['n_side'])][1][l] * w_pix[str(freq2_param['n_side'])][1][l] * bl_data_1[0, l] * bl_data_2[0, l] * 2.7255**2
						#BB
						c_l_th[l, i, j, im, jm, 1] = all_cls_th[l, 2] * w_pix[str(freq1_param['n_side'])][1][l] * w_pix[str(freq2_param['n_side'])][1][l] * bl_data_1[1, l] * bl_data_2[1, l] * 2.7255**2
		j += 1
	i += 1
print(c_l_th.shape)
np.save('pre_proc/beam_corrected_lcdm_spectra_{}.npy'.format('_'.join(maps_param.keys())), c_l_th)