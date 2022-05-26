import numpy as np
import tools_fast as tools
import os.path
import sys
sys.path.append('parameter_files/')
from hfi_lfi_wmap_eb import maps_param

# If you only work with HFI
# from hfi_eb import maps_param


# This file takes the observed power spectra from PolSpice and puts them into a format
# that cb.py can read. It does not alter the spectra

l_max = 1520
nob = len(maps_param.keys())

size = tools.size(nob)

main_data_set = 'cl_wmap_hfi_lfi'

# If you only work with HFI
# main_data_set = 'cl_hfi'

each_folder = 'mask_percent_{}'

output_folder = main_data_set

print('Data set:', main_data_set)

galactic_mask_percent_list = ['0', '30']

for h, mask in enumerate(galactic_mask_percent_list):
	cl = np.zeros((nob, nob, 2, 2, 3, l_max+1))
	i = 0
	for freq1, freq1_param in maps_param.items():
		j = 0
		for freq2, freq2_param in maps_param.items():
			for im in range(2):
				for jm in range(2):
					filename = "generate_observed_power_sepctra/{}/{}/{}{}_{}{}.dat".format(
						main_data_set,
						each_folder.format(mask),
						freq1.replace('_', ''),
						freq1_param['data_split'][im],
						freq2.replace('_', ''),
						freq2_param['data_split'][jm])
					reverse = False
					if os.path.isfile(filename) == False:
						reverse = True
						filename = "generate_observed_power_sepctra/{}/{}/{}{}_{}{}.dat".format(
							main_data_set,
							each_folder.format(mask),
							freq2.replace('_', ''),
							freq2_param['data_split'][jm],
							freq1.replace('_', ''),
							freq1_param['data_split'][im])
					print(filename)
					curr_file = np.loadtxt(filename, dtype = float)

					# WMAP uses mK
					if freq1_param['data_set'] == 'wmap':
						curr_file[:, 1:] *= 1e-3
					if freq2_param['data_set'] == 'wmap':
						curr_file[:, 1:] *= 1e-3
					cl[i, j, im, jm, 0, :] = curr_file[:l_max+1, 2] #EE
					cl[i, j, im, jm, 1, :] = curr_file[:l_max+1, 3] #BB

					if reverse:
						cl[i, j, im, jm, 2, :] = curr_file[:l_max+1, 9] #BE
					else:
						cl[i, j, im, jm, 2, :] = curr_file[:l_max+1, 6] #EB

			j += 1
		i += 1

	output_file = 'pre_proc/{}/cl_{}_freq_{}.npy'.format(output_folder, each_folder.format(mask), '_'.join(maps_param.keys()))
	np.save(output_file, cl)
	print(output_file)

	
	
	
	

