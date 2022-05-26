import numpy as np
from numpy import cos,sin,pi,tan,sqrt
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import emcee
import corner
import os.path
import tqdm
import tools_fast as tools
import scipy
import time

class cosmic_birefringence:
	# This is the class that stores all the variables and formats it correctly.
	# tools_fast.py does the actual math

	def __init__(self, parameters):
		self.params = parameters
		nob = self.params['nob']
		l_max = self.params['l_max']
		l_min = self.params['l_min']

		# Make sure we only deal with full ell-bins
		assert(((l_max - l_min) / self.params['l_bin']).is_integer())

		# Degrees to radians.
		if self.params['frequency_dependent_beta']:
			if self.params['beta_model'] == 'power':
				# Assume beta = beta_0 * (nu/nu_0)^n for some given nu_0
				self.params['initial_beta'] = np.array([self.params['initial_beta'] * pi/180, 0])
			elif self.params['beta_model'] == 'individual':
				# Samples each beta individually for each frequency
				self.params['initial_beta'] = np.ones(nob) * self.params['initial_beta'] * pi/180
		else:
			self.params['initial_beta'] = self.params['initial_beta'] * pi/180

		# Initialize the initial_alphas
		self.params['initial_alphas'] = np.concatenate((self.params['initial_alphas'], self.params['initial_alphas']))
		if len(self.params['alpha_no_split_list']) >= 1:
			self.params['initial_alphas'] = self.params['initial_alphas'][:-len(self.params['alpha_no_split_list'])]
		print('Number of alphas:', len(self.params['initial_alphas']))
			
		# Get f_sky	
		if self.params['data_set'] == 'cl_wmap_hfi_lfi':
			index = self.params['index']
			if index != -1:
				print('Using hfi/lfi ps mask with CO. Index:', index)
				cur_f = np.load('generate_observed_power_spectra/masks/f.npy')[index]
				self.params['f'] = np.ones((nob, 2)) * cur_f
			else:
				print('WARNING: USING f=1')
				self.params['f'] = np.ones((4, 2))
		
		size = tools.size(self.params['nob'])
		self.params['size'] = size
		self.number_of_bins = int((l_max-l_min)/self.params['l_bin'])
		self.import_cls()

		# Load LCDM spectra
		lcdm_filename = 'pre_proc/beam_corrected_lcdm_spectra_{}.npy'.format('_'.join(self.params['alpha_labels']))
		print('Importing beam smoothed LCDM cl. Filename:', lcdm_filename)
		self.c_l_th = np.load(lcdm_filename)

		# Are we trying to model the EB of the foreground?
		if self.params['psi_l'] != '':
			print('Loading', self.params['psi_l'])
			self.psi_l = np.load(self.params['psi_l'])[:, self.params['index']]
		else:
			# No EB modeling. Make psi_l into a one element list so that Numba can compile
			self.psi_l = np.array([])
		
		# Get covariance matrix		
		cov_filename = 'pre_proc/cov_bin_{}_{}_{}_lmin{}_lmax{}.npy'.format(self.params['data_set'], self.params['cl_folder'], '_'.join(self.params['alpha_labels']), l_min, l_max)
		calculate_cov = os.path.isfile(cov_filename) == False or nob <= 7
		if calculate_cov:
			# Here we calculate the covariance matrix. With many bands, this will take some time
			start = time.time()
			print('Calculating covariance matrix')
			self.obs_cov = tools.get_covariance(size, l_max, nob, self.observed_cross_spectra, self.params['f'])
			print(time.time() - start, ' seconds to calculate the covariance matrix')
		else:
			# To make things easier, we still calculate the binned obs_cov from array of zeros
			# even though we will load the binned one later
			self.obs_cov = np.zeros((size*3, size*3, l_max+1))

		# Set and store the binned spectra
		n_bins = self.number_of_bins	
		self.c_l_o_bin = np.zeros((n_bins, size * 3))
		self.c_l_th_bin = np.zeros((n_bins, size * 2))
		self.obs_cov_bin =  np.zeros((n_bins, size * 3, size * 3))
		for l_b in range(n_bins):
			c_l_o_bin, c_l_th_bin, obs_cov_bin = self.get_bin(l_b)
			self.c_l_o_bin[l_b, :] = c_l_o_bin
			self.c_l_th_bin[l_b, :] = c_l_th_bin
			self.obs_cov_bin[l_b, :, :] = obs_cov_bin
		
		if calculate_cov and nob >= 8:
			# Dont save the binned covariance matrix if we dont have many bands
			np.save(cov_filename, self.obs_cov_bin)
		elif calculate_cov == False:
			# Load if it exists to save time
			print('Loading saved covariance matrix. Filename:', cov_filename)
			self.obs_cov_bin = np.load(cov_filename)

	def import_cls(self):
		# Load the observed power spectra
		cl_filename = 'pre_proc/{}/cl_{}_freq_{}.npy'.format(self.params['data_set'], self.params['cl_folder'], '_'.join(self.params['alpha_labels']))
		if os.path.isfile(cl_filename):
			print('Importing observed cl. Filename:', cl_filename)
			self.observed_cross_spectra = np.load(cl_filename)
		else:
			print('Cant find:', cl_filename)
			raise Exception('Cl file does not exist')			
	
	def get_bin(self, l_b):
		# Here we get the binned power spectra and covariance matrices
		l_bin = self.params['l_bin']
		nob = self.params['nob']
		size = self.params['size']
		obs_cov = self.obs_cov
		ocs = self.observed_cross_spectra
		c_l_th = self.c_l_th

		c_l_o_bin = np.zeros(size*3)
		c_l_th_bin = np.zeros(size*2)
		obs_cov_bin = np.zeros((size*3, size*3))
		
		# Loop through all bins
		for l_0 in range(l_bin):
			l = self.params['l_min'] + l_bin*l_b + l_0
			c_l_o_vec = np.zeros(size*3)
			c_l_th_vec = np.zeros(size*2)

			v_i = 0
			for i in range(nob):
				for j in range(nob):
					for im in range(2):
						for jm in range(2):
							if tools.do_not_count_this_eq(i, j, im, jm):
								# Do not add auto spectra
								continue
							c_l_th_vec[2*v_i:2*v_i+2] = c_l_th[l, i, j, im, jm, :2]
							c_l_o_vec[3*v_i:3*v_i+3] = ocs[i, j, im, jm, :3, l]

							v_i += 1
			obs_cov_bin += obs_cov[:, :, l]
			c_l_o_bin += c_l_o_vec
			c_l_th_bin += c_l_th_vec

		c_l_th_bin /= l_bin
		c_l_o_bin /= l_bin
		obs_cov_bin /= l_bin**2
		return c_l_o_bin, c_l_th_bin, obs_cov_bin

	def log_prob(self, theta, beta_i='none'):
		# The log likelihood used when sampling the parameters

		nob = self.params['nob']
		number_of_alphas = len(self.params['initial_alphas'])
		if self.params['frequency_dependent_beta']:
			if self.params['beta_model'] == 'power':
				# Power model for beta
				beta_0 = theta[0]
				n = theta[1]
				nu_0 = self.params['nu_0']

				#nu = np.array([float(i) for i in self.params['alpha_labels']])
				nu = self.params['frequency_list_beta_power_law']
				beta = beta_0 * np.power(nu/nu_0, n)

				# One beta for each split
				beta = np.tile(beta, 2)

				# All our maps are treated as half mission or detector split maps of the same channel.
				# But K and Ka of WMAP only has one map each. The only problem with this approach is
				# when we look at the frequency dependence of beta. So here we set different freq for
				# the K and Ka channels.
				if 'frequency_beta_split_AB' in self.params:
					for key in self.params['frequency_beta_split_AB']:
						beta_A = nob * 0 + key
						beta_B = nob * 1 + key
						nu_A = self.params['frequency_beta_split_AB'][key][0]
						nu_B = self.params['frequency_beta_split_AB'][key][1]
						beta[beta_A] = beta_0 * np.power(nu_A/nu_0, n)
						beta[beta_B] = beta_0 * np.power(nu_B/nu_0, n)
				
				alpha = theta[2 : number_of_alphas + 2]
				A = theta[number_of_alphas + 2 :]
				lp = self.flat_prior(theta[:2], alpha, A)
			else:
				# Individual beta for each band
				beta = theta[:nob]
				alpha = theta[nob : number_of_alphas+nob]
				A = theta[number_of_alphas+nob:]
				lp = self.flat_prior(beta, alpha, A)

				# One beta for each split
				beta = np.tile(beta, 2)
		else:
			if beta_i == 'none':
				# We sample beta
				beta = np.array([theta[0]])
				alpha = theta[1 : number_of_alphas+1]
				A = theta[number_of_alphas + 1 : ]
				lp = self.flat_prior(beta, alpha, A)
			else:
				# We do not sample beta. It is set by beta_i
				beta = np.array([beta_i])
				alpha = theta[:number_of_alphas]
				A = theta[number_of_alphas : ]
				lp = self.flat_prior(beta, alpha, A)
			
		#30 and 44GHz use time splits and not detector splits
		#So they should have the same miscalibration angle so we need to format alpha
		correct_format_alpha = self.get_correct_alpha_format(alpha)

		if not np.isfinite(lp):
			return -np.inf

		if self.params['psi_l'] == '':
			# No EB modeling
			prob = tools.likelihood_prob(self.c_l_o_bin, self.c_l_th_bin, self.obs_cov_bin, self.params['with_ln_det_cov'], bin_range=self.number_of_bins, nob=nob, alpha=correct_format_alpha, beta=beta)
		else:
			# We model the dust EB
			prob = tools.likelihood_prob_model_eb(self.c_l_o_bin, self.c_l_th_bin, self.obs_cov_bin, self.params['with_ln_det_cov'], bin_range=self.number_of_bins, nob=nob, alpha=correct_format_alpha, beta=beta, psi_l=self.psi_l, A=A, turn_off_for_lowest_indices = self.params['turn_off_for_lowest_indices'])

		return prob

	def run_beta_grid(self, plot=True):
		# If we set alpha=0, it is easiest to just grid beta instead of sampling it since it is just one parameter
		num = 400
		beta_x = np.linspace(0.0*pi/180, 0.6*pi/180, num)
		p = np.zeros(num)
		for x in tqdm.tqdm(range(num)):
			alpha = np.zeros(22)
			omega = np.concatenate(([beta_x[x]], alpha))
			
			prob = self.log_prob(omega)
			p[x] = prob

		p = np.exp(p-max(p))
		beta_x = beta_x*180/pi
		norm = scipy.integrate.simps(p, x=beta_x)
		mean = beta_x[np.argmax(p)]
		sigma = norm/(np.max(p)*sqrt(2*pi))

		print('Assume no miscalibration: Beta is mean:', mean, 'std:', sigma)

		if plot:
			plt.plot(beta_x, p, label='actual run')
			plt.legend()
			plt.savefig('tmp/__beta_grid.pdf')
	
		return mean, sigma

	def flat_prior(self, beta, alpha, A):
		# Flat priors for our sampled parameters. Most importantly, we use A_ell>=0
		p = 0
		ap = self.params['flat_alpha_prior']*pi/180

		for i in range(len(alpha)):
			if abs(alpha[i]) > ap:
				p = -np.inf

		for i in range(len(A)):
			if A[i] < 0 or A[i] > self.params['upper_bound_A']: p = -np.inf

		if len(beta) == 2 and self.params['beta_model'] == 'power':
			# Priors on n in beta = beta_0 * (nu / nu_0)**n
			if np.abs(beta[1]) > self.params['flat_n_prior']: p = -np.inf

		return p

	def return_likelihood(self, beta='none'):
		# This minimizes the likelihood function
		a_0 = 2*pi/180
		angles = len(self.params['initial_alphas'])
		if beta == 'none':
			angles += 1

		n_params = angles + self.params['number_of_A']
		params = np.zeros(n_params) + 0.05
		
		bnds = np.zeros((n_params, 2))
		bnds[:angles, 0] = -a_0
		bnds[:angles, 1] = a_0
		bnds[angles:, 1] = self.params['upper_bound_A']
		bnds = tuple(map(tuple, bnds))
		if beta=='none':
			# Fit beta as well
			func = lambda x: -2*self.log_prob(x)
		else:
			# Beta is not fit for put set by hand
			func = lambda x: -2*self.log_prob(x, beta_i=beta)

		time_start = time.time()
		min_params = scipy.optimize.minimize(func, params, method = 'L-BFGS-B', bounds=bnds, options={'ftol':1e-12}).x
		print('Time taken:', time.time() - time_start)
		return func(min_params), min_params

	def run_sampler(self, it=500, burnin=100, fig_name='default', no_beta = False):
		# This starts the sampling process
		print('Running:', fig_name)
		print(self.params)
		params = self.params
		nob = params['nob']

		# Get correct alpha labels
		alpha_labels = []
		for element in params['alpha_no_split_list']:
			alpha_labels.append(r"$\alpha_{{{}}}$".format(element))
		for j in range(2):
			for i in range(nob):
				if params['alpha_labels'][i] in params['alpha_no_split_list']:
					continue
				alpha_labels.append(r"$\alpha_{{{}{}}}$".format(params['alpha_labels'][i], 'A' if j==0 else 'B'))

		# Get correct beta labels
		if no_beta:
			labels = alpha_labels
			start_pos = params['initial_alphas']
		else:
			if params['frequency_dependent_beta']:
				if params['beta_model'] == 'power':
					labels = np.concatenate((np.array([r"$\beta_{0}$", r"$n$"]), alpha_labels))
				else:
					beta_labels = np.array([r"$\beta_{{{}}}$".format(params['alpha_labels'][i]) for i in range(nob)])
					labels = np.concatenate((beta_labels, alpha_labels))
				start_pos = np.concatenate((params['initial_beta'], params['initial_alphas']))
			else:
				labels = np.concatenate((np.array([r"$\beta$"]), alpha_labels))
				start_pos = np.concatenate((np.array([params['initial_beta']]), params['initial_alphas']))
		
		if len(self.psi_l) > 0:
			# We are sampling As
			A_labels = np.array([r"$A_{{{}}}$".format(i) for i in range(params['number_of_A'])])
			labels = np.concatenate((labels, A_labels))
			start_pos = np.concatenate((start_pos, np.ones(params['number_of_A'])*0.05))
		print('Labels:', labels)

		# Initial guesses for the MCMC
		start_stack = start_pos
		for i in range(len(start_pos)*2):
			start_stack = np.vstack((start_stack, start_pos + 0.1*pi/180*np.random.randn(len(start_pos))))
		nwalkers, ndim = start_stack.shape

		# Set initial values for beta and A_ell
		if params['frequency_dependent_beta'] and params['beta_model'] == 'power':
			start_stack[:, 1] = np.random.normal(loc = 0, scale = 0.5, size = nwalkers)
		if len(self.psi_l) > 0:
			start_stack[:, -params['number_of_A']:] = np.random.normal(loc = 0.05, scale = 0.015, size = (nwalkers, params['number_of_A']))

		# Do we need to worry about convergence? If the code is slow, we need to make sure we have
		# reached convergence.
		if self.params['backend']:
			# Much data means slower convergence time
			# When the code is slow we will regularly check when the chain converges

			filename = 'chains/{}.h5'.format(fig_name)
			backend = emcee.backends.HDFBackend(filename)
			backend.reset(nwalkers, ndim)

			# Initialize the sampler
			sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_prob, backend=backend)

			# We'll track how the average autocorrelation time estimate changes
			index=0
			autocorr = np.empty(it)

			# This will be useful to testing convergence
			old_tau = np.inf

			start_stack = sampler.get_last_sample()

			for sample in sampler.sample(start_stack, iterations=it, progress=True):
				# Only check convergence every 100 steps
				if sampler.iteration % 100:
						continue

				# Compute the autocorrelation time so far
				# Using tol=0 means that we'll always get an estimate even
				# if it isn't trustworthy
				tau = sampler.get_autocorr_time(tol=0)
				autocorr[index] = np.mean(tau)
				index += 1

				# Check convergence
				converged = np.all(tau * 100 < sampler.iteration)
				converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
				'''if converged:
						tau = sampler.get_autocorr_time()
						cal_burnin = int(2 * np.max(tau))
						thin = int(0.5 * np.min(tau))
						print('Finished!')
						print('Burning:', cal_burnin)
						print('Thin:', thin)
						flat_samples = sampler.get_chain(discard=cal_burnin, thin=thin, flat=True, )
						break'''
				old_tau = tau

				n = 100 * np.arange(1, index + 1)
				y = autocorr[:index]
				plt.figure()
				plt.plot(n, n / 100.0, "--k")
				plt.plot(n, n / 50.0, "--k")
				plt.plot(n, y)
				plt.xlim(0, n.max())
				plt.ylim(0, y.max() + 0.1 * (y.max() - y.min()))
				plt.xlabel("number of steps")
				plt.ylabel(r"mean $\hat{\tau}$")
				plt.savefig('chains/convergence_{}.pdf'.format(fig_name))
				plt.close()

		else:
			sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_prob)
			sampler.run_mcmc(start_stack, it, progress=True)
			flat_samples = sampler.get_chain(discard=burnin, thin=5, flat=True)
	
		# Rad to deg
		if self.params['number_of_A'] == 0:
			# Set everything to degree
			flat_samples *= 180/pi
		else:
			# Set everything to degree except A_\ell
			flat_samples[:, :-self.params['number_of_A']] *= 180/pi
		if self.params['frequency_dependent_beta'] and self.params['beta_model'] == 'power':
			# This is n in beta = beta_0 * (nu/nu_0)**n
			# Should not be converted to rad
			flat_samples[:, 1] *= pi/180

		corner.corner(
		    flat_samples, labels=labels, show_titles=True, fontsize=12, color='g', title_fmt='.3f'
		);
		number_of_beta = len(self.params['initial_beta']) if type(self.params['initial_beta']) is list else 1
		np.save('chains/full_chain_{}_A_{}_beta_{}_samples_{}_{}.npy'.format(self.params['number_of_A'], number_of_beta, self.params['cl_folder'], it, fig_name), flat_samples)
		plt.savefig('img/{}/{}_{}_it_{}.pdf'.format(self.params['data_set'], self.params['cl_folder'], fig_name, it))
		
		print('Full Analysis Beta:', np.mean(flat_samples[:, 0]), 'Std:', np.std(flat_samples[:, 0]))
		print('Autocorrelation length:', sampler.get_autocorr_time(tol=0))

		mean, upper, lower = np.zeros(ndim), np.zeros(ndim), np.zeros(ndim)
		for i in range(ndim):
			mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
			mean[i] = mcmc[1]
			upper[i] = mcmc[2] - mcmc[1]
			lower[i] = mcmc[1] - mcmc[0]

		return mean, upper, lower

	def get_correct_alpha_format(self, alpha):
		# This formats alpha into the format that can be read by tools_fast methods so that there
		# is one miscalibration angle per map.
		# For example, 30 and 44 channels have one miscalibration angle for the two hald mission maps
		# we therefore have to repeat the same misclibration angle for both maps.

		nob = self.params['nob']
		params = self.params
		correct_format_alpha = np.zeros(2*nob)
		number_of_no_alpha_split_freqs = len(params['alpha_no_split_list'])
		for i in range(nob):
			if params['alpha_labels'][i] in params['alpha_no_split_list']:
				correct_format_alpha[i] = alpha[i]
				correct_format_alpha[i+nob] = alpha[i]
			else:
				for j in range(2):
					index = j*nob + i
					if j == 0:
						correct_format_alpha[index] = alpha[index]
					elif j == 1:
						correct_format_alpha[index] = alpha[index - number_of_no_alpha_split_freqs]
		return correct_format_alpha

	# Plotting functions
	#--------------------

	def plot_eb_mean(self, beta, alpha, A, best_fit_params):
		# This function make the EB plot in the WMAP + Planck paper

		nob = self.params['nob']

		eb_mean, eb_std, chi_squared_all, chi_squared = tools.get_inverse_variance_mean_std(
			self.c_l_o_bin,
			self.c_l_th_bin,
			self.obs_cov_bin,
			bin_range=self.number_of_bins,
			nob=nob,
			alpha=np.zeros(22),
			beta=np.array([0]),
			psi_l=self.psi_l,
			A=np.zeros(4),
			turn_off_for_lowest_indices = self.params['turn_off_for_lowest_indices'],
			with_band = False
		)
		print('chi squared:', chi_squared_all, chi_squared)
		l_min = self.params['l_min']
		l_max = self.params['l_max']
		bin_size = self.params['l_bin']
		ell_bin = np.arange(l_min, l_max, bin_size) + bin_size/2

		plt.rc('text', usetex=True)
		fig, axs = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios':[2,1]})
		fig.tight_layout(h_pad=2)
		axs[0].errorbar(ell_bin, ell_bin * eb_mean*1e12, yerr = ell_bin * eb_std*1e12, fmt='.', color='black')

		# We make a histogram out of 2500 samples
		n = 2500
		eb_beta, eb_alpha = np.zeros((n, 72)), np.zeros((n, 72))
		for i in tqdm.tqdm(range(n)):
			eb_alpha[i, :], eb_beta[i, :], _, _  = tools.get_inverse_variance_mean_std(
				self.c_l_o_bin,
				self.c_l_th_bin,
				self.obs_cov_bin,
				bin_range=self.number_of_bins,
				nob=nob,
				alpha=self.get_correct_alpha_format(alpha[i*10, :]),
				beta=np.array([beta[i*10]]),
				psi_l=self.psi_l,
				A=A[i*10, :],
				turn_off_for_lowest_indices = self.params['turn_off_for_lowest_indices'],
				with_band = True
			)

		cmb_mean = np.mean(eb_beta, axis=0)
		cmb_lower = cmb_mean - np.std(eb_beta, axis=0)
		cmb_upper = cmb_mean + np.std(eb_beta, axis=0)
		mis_mean = np.mean(eb_alpha, axis=0)
		mis_lower = mis_mean - np.std(eb_alpha, axis=0)
		mis_upper = mis_mean + np.std(eb_alpha, axis=0)

		axs[0].fill_between(ell_bin, ell_bin *cmb_lower*1e12, ell_bin*cmb_upper*1e12, color='blue', alpha=0.4)
		axs[0].plot(ell_bin, ell_bin * cmb_mean *1e12, color='blue', label='Cosmic birefringence')

		axs[0].fill_between(ell_bin, ell_bin *mis_lower * 1e12, ell_bin*mis_upper*1e12, color='red', alpha=0.4)
		axs[0].plot(ell_bin, ell_bin * mis_mean * 1e12, color='red', label='Miscalibration angle')

		axs[0].set_ylabel(r'EB power spectrum, $\ell C_{\ell}^{EB}$ [$\mu K^2$]')
		axs[0].axhline(y = 0, color = 'grey', linestyle = '-')
		axs[0].set_xlim([30, 1511])
		axs[0].set_ylim([-4.2 * 1e-3, 4.2 * 1e-3])
		axs[0].tick_params(axis='y', rotation = 90)
		axs[0].set_title(r'Stacked observed $EB$ power spectrum')
		axs[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
		axs[0].text(1300, 3.53*1e-3, r'$f_{\textrm{sky}}=0.92$')

		best_fit_beta = best_fit_params[0]
		best_fit_alpha_correct_format = self.get_correct_alpha_format(best_fit_params[1:23])
		best_fit_A = best_fit_params[-4:]
		print('best fit - alpha', best_fit_alpha_correct_format*180/np.pi, 'beta', best_fit_beta*180/np.pi, 'A', best_fit_A)
		rms_mean, rms_std, chi_squared_all, chi_squared_stacked = tools.get_inverse_variance_mean_std(
			self.c_l_o_bin,
			self.c_l_th_bin,
			self.obs_cov_bin,
			bin_range=self.number_of_bins,
			nob=nob,
			alpha=best_fit_alpha_correct_format,
			beta=np.array([best_fit_beta]),
			psi_l=self.psi_l,
			A=best_fit_A,
			turn_off_for_lowest_indices = self.params['turn_off_for_lowest_indices'],
			with_band = False
		)
		axs[0].plot(ell_bin, ell_bin*(eb_mean-rms_mean)*1e12, color='limegreen', linewidth=4, label='Best-fit total')
		axs[0].legend(loc='lower left')

		print('chi squared:', chi_squared_all, chi_squared_stacked)
		axs[1].set_xlabel(r'$\ell$')
		axs[1].set_ylabel(r'Residual, $\ell C_{\ell}^{EB}$ [$\mu K^2$]')
		axs[1].errorbar(ell_bin, ell_bin * rms_mean*1e12, yerr = ell_bin * rms_std*1e12, fmt='.', color='black')
		axs[1].set_xlim([30, 1511])
		axs[1].set_ylim([-4.2 * 1e-3, 4.2 * 1e-3])
		axs[1].set_yticks([-4 * 1e-3, -2 * 1e-3, 0, 2 * 1e-3, 4 * 1e-3])
		axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
		axs[1].tick_params(axis='y', rotation = 90)
		axs[1].set_title(r'Residual with respect to the model')
		axs[1].axhline(y = 0, color = 'grey', linestyle = '-')
		plt.savefig('results/wmap_lfi_hfi_stacked_eb/stacked_eb_ell_5.pdf', bbox_inches='tight')

	def plot_cls(self):
		# Plot the power spectra
		ocs = self.observed_cross_spectra
		l_bin = self.params['l_bin']
		nob = self.params['nob']
		l_min = self.params['l_min']
		alpha_labels = self.params['alpha_labels']
		for i in range(nob):
			for j in range(nob):
				points = self.number_of_bins
				EB_o, EB_o_err, EE_BB_o, EE_BB_th = np.zeros(points), np.zeros(points), np.zeros(points), np.zeros(points)
				for l_b in range(points):
					EB_o_cur, EE_BB_o_cur, EE_BB_th_cur = np.zeros(l_bin), np.zeros(l_bin), np.zeros(l_bin)
					for l_ in range(l_bin):
						l = l_min + l_bin*l_b + l_
						EB_o_cur[l_] = l*ocs[i, j, 0, 1, 2, l]
						EE_BB_o_cur[l_] = l*(ocs[i, j, 0, 1, 0, l] - ocs[i, j, 0, 1, 1, l])
						EE_BB_th_cur[l_] = l*(self.c_l_th[l, i, j, 0, 1, 0] - self.c_l_th[l, i, j, 0, 1, 1])
					EB_o_err[l_b] = np.std(EB_o_cur)/np.sqrt(l_bin)
					EB_o[l_b] = np.sum(EB_o_cur)/l_bin
					EE_BB_o[l_b] = np.sum(EE_BB_o_cur)/l_bin
					EE_BB_th[l_b] = np.sum(EE_BB_th_cur)/l_bin
							
				l = np.arange(l_min, self.params['l_max'], l_bin)
				plt.figure()
				plt.title('{}Ax{}B'.format(alpha_labels[i], alpha_labels[j]))
				#plt.plot(l, EB_o, '.', label='EB_0 observed {}GHz'.format(self.alpha_labels[i]))		
				plt.errorbar(l, EB_o, fmt='.', color='black',yerr=EB_o_err, label='EB_0 observed')
				rotated_ee_bb_o = tan(4*0.3 * pi/180)/2*EE_BB_o
				plt.plot(l, rotated_ee_bb_o, color='red', label='Rotated EE-BB observed')
				plt.plot(l, sin(4*0.3 * pi / 180)/(2*cos(4*0.3*pi/180))*EE_BB_th, color='blue', label='Rotated EE-BB theory')
				plt.legend()
				plt.ylim([-max(1.3 * rotated_ee_bb_o), max(1.3 * rotated_ee_bb_o)])
				plt.savefig('tmp/freq_rotated_{}A_{}B.pdf'.format(alpha_labels[i], alpha_labels[j]))		
				plt.figure()
				plt.plot(l, EE_BB_o, color='red', label='EE-BB observed')
				plt.plot(l, EE_BB_th, color='blue', label='EE-BB theory')
				plt.plot(l, EB_o, color='green', label='EB observed')
				plt.title('{}Ax{}B'.format(alpha_labels[i], alpha_labels[j]))
				plt.legend()
				plt.savefig('tmp/freq_{}A_{}B.pdf'.format(alpha_labels[i], alpha_labels[j]))