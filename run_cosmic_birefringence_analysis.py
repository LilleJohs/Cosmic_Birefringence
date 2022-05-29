from cb import *
import copy
### If on owl
import os
os.environ["OMP_NUM_THREADS"] = "1"
###

folder = 'cl_wmap_hfi_lfi'
mask_percent = ['0', '30']

print('Doing', folder)
freqs = ['030', '044', 'K', 'Q', 'V', '070', 'W', 'W_', '100', '143', '217', '353']

# If you are running this for the first time, I suggest doing HFI only. It is much faster.
# Running WMAP + Planck will take you many days
# freqs = ['100', '143', '217', '353']

nob = len(freqs)
params = {
	'alpha_labels': freqs,
	'initial_beta': 0.0,
	'nob': nob,
	'initial_alphas': np.zeros(nob),

	#These bands have only one miscalibration angle
	'alpha_no_split_list': ['030', '044'],
	
	'l_min': 51,
	'l_max': 1491,
	'l_bin': 20,
	'data_set': folder,
	'flat_beta_prior': 5.0, # |beta| < 5 deg
	'flat_alpha_prior': 5.0, # |alpha_i| < 5
	'cl_folder': '',
	'index': -1,
	'with_ln_det_cov': True, # is ln(det(M_b)) in the likelihood?

	# Backend should be true if doing WMAP + Planck since it is slow and you have to be careful about
	# convergence. If you only do HFI, it should be pretty quick to do a lot of samples quickly and you can
	# turn off backend
	'backend': True,


  # Beta model
	# To sample the frequency dependence of beta, change this to true and uncomment the following lines
	'frequency_dependent_beta': False,
	#'beta_model': 'power', # sample n and beta_0 in the power law beta = beta_0(nu/nu_0)^n
	#'nu_0': 150, # in n
	#'flat_n_prior': 3, # Flat prior on |n|
	#'frequency_list_beta_power_law': np.array([30, 44, 23, 41, 61, 70, 94, 94, 100, 143, 217, 353]), # list of frequencies, make sure the order is correct
	#'frequency_beta_split_AB': { # K and Ka are grouped together as one freuqnecy band. So we need to tell cb.py that they are different frequencies
	#	2: np.array([23, 33])
	#},
	
  # EB-modeling parameters
	'psi_l': 'psi_l_sigma15.npy', # The psi_ell file that you need to generate
	'number_of_A': 4, # Number of A_ell bins

	# The x number of bands are not dust dominated
	'turn_off_for_lowest_indices': 6, # In the list of bands (freq variable above), the first 6 bands are synch dominated. So we dont model dust EB for them
	'upper_bound_A': 1.0, # 0 <= A_ell <= 1

	# If you dont want to model the EB of dust, uncomment these lines and comment the lines above
	#'psi_l': '',
	#'number_of_A': 0,
}

# Make sure we have the correct number of sampled parameters
# 30 and 44 GHz have one miscalibration angle
number_of_sampled_params = 2*params['nob'] - len(params['alpha_no_split_list']) + params['number_of_A']
if params['frequency_dependent_beta']:
	if params['beta_model'] == 'power':
		betas = 2
	else:
		betas = params['nob']
else:
	betas = 1

number_of_sampled_params += betas
saved = np.zeros((len(mask_percent), number_of_sampled_params, 3))
print(saved.shape)

iterations = [100000, 100000]
for i, mask in enumerate(mask_percent):
	print('Galactic mask,', mask)
	params['cl_folder'] = 'mask_percent_{}'.format(mask)
	params['index'] = i

	a = cosmic_birefringence(parameters=copy.deepcopy(params))

	fig_name = 'wmap_planck_dust_eb_modelled'

	cur_mean, cur_upper, cur_lower = a.run_sampler(it=iterations[i], burnin=5000, fig_name=fig_name, no_beta=False)
	saved[i, :, 0] = cur_mean
	saved[i, :, 1] = cur_upper
	saved[i, :, 2] = cur_lower

np.save('results.npy', saved)
print('*** JOB FINISHED ***')