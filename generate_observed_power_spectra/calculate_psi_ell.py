import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

# This files generate the psi_ell angle, which is the multipole
# averaged misalignment between local magnetic fields and filamentary
# dust structures. Psi_ell is used to model the dust EB power spectrum

l_max = 1490
l_min = 51

main_data_set = 'cl_wmap_hfi_lfi'
# If you only work with HFI
# main_data_set = 'cl_hfi'

each_folder = 'mask_percent_{}'

def get_psi(cl):
	# sigma = 1.5 was used for the paper, but you can make psi_ell more or less smooth by changing sigma
	sigma = 1.5

	# First we smooth the binne, symmetric power spectra
	smooth_tb = gaussian_filter1d(cl[1, :], sigma = sigma)
	smooth_te = gaussian_filter1d(cl[0, :], sigma = sigma)

	# Then we find psi_ell
	psi_ell = 0.5 * np.arctan(smooth_tb / smooth_te)

	return psi_ell

def get_binned_cl(mask):
	# Bin the power spectra

	data_set = main_data_set
	filename = "{}/{}/353{}_353{}.dat".format(data_set, each_folder.format(mask), 'A', 'B')

	curr_file = np.loadtxt(filename, dtype = float)

	cl = np.zeros((2, l_max+1))
	cl[0, :] = (curr_file[:l_max+1, 4] + curr_file[:l_max+1, 7])/2  #Symmetric TE
	cl[1, :] = (curr_file[:l_max+1, 5] + curr_file[:l_max+1, 8])/2  #Symmetric TB

	bin_size = 20
	nr_bin = int((l_max+1 - l_min)/bin_size)
	binned_cl = np.zeros((2, nr_bin))
	
	for i in range(2):
		for j in range(nr_bin):
			binned_cl[i, j] = cl[i, l_min + bin_size*j : l_min+bin_size*(j+1)].mean()
	return binned_cl

l_array = np.arange(l_min, l_max, 20)
psi_ell_list = np.zeros((72, 2))
masks = ['0', '30']

for i, mask in enumerate(masks):
	binned_cl = get_binned_cl(mask)
	psi_ell_list[:, i] = get_psi(binned_cl)

np.save('psi_l_sigma15.npy', psi_ell_list)