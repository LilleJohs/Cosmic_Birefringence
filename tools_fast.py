from numba import njit, prange
import numpy as np
from numpy import cos, sin, tan

@njit
def get_A_B_unprimed(nob, alpha, beta):
	# This gets A and B matrices assuming no EB modeling of the foreground

	size_ = size(nob)
	A_cur = np.zeros((size_, size_*3))
	B_cur = np.zeros((size_, size_*2))
	v_i = 0
	for i in range(nob):
		for j in range(nob):
			for im in range(2):
				for jm in range(2):
					if do_not_count_this_eq(i, j, im, jm):
						continue

					cur_alpha_i = nob*im + i
					cur_alpha_j = nob*jm + j
					if len(beta) == 1:
						A_block_cur, B_block_cur = A_and_B(2*alpha[cur_alpha_i], 2*alpha[cur_alpha_j], 2*np.array([beta[0], beta[0]]))
					else:
						A_block_cur, B_block_cur = A_and_B(2*alpha[cur_alpha_i], 2*alpha[cur_alpha_j], 2*np.array([beta[i], beta[j]]))
					A_hori_start = v_i * 3
					B_hori_start = v_i * 2
					A_cur[v_i, A_hori_start:A_hori_start+3] = A_block_cur
					B_cur[v_i, B_hori_start:B_hori_start+2] = B_block_cur

					v_i += 1
	return A_cur, B_cur

@njit
def get_A_B_primed_k(nob, alpha, beta, psi_l, A, turn_off_for_lowest_indices, k):
	# Get the A and B matrices for one ell bin. Since psi_l depends on the multipole l, 
	# we need to calculate A and B matrices for each ell bin (k)
	size_ = size(nob)

	A_ = get_A_ell(k, A)

	A_p_k = np.zeros((size_, size_*3))
	B_p_k = np.zeros((size_, size_*2))
	
	v_i=0
	for i in range(nob):
		for j in range(nob):
			for im in range(2):
				for jm in range(2):
					if do_not_count_this_eq(i, j, im, jm):
						continue

					cur_alpha_i = nob*im + i
					cur_alpha_j = nob*jm + j

					A_hori_start = v_i * 3
					B_hori_start = v_i * 2

					# Are we sampling one beta or a frequency dependent beta?
					if len(beta) == 1:
						beta_sample = np.array([beta[0], beta[0]])
					else:
						beta_sample = np.array([beta[i], beta[j]])

					if i < turn_off_for_lowest_indices or j < turn_off_for_lowest_indices:
						# These frequencies do not have any EB from dust
						A_block_cur, B_block_cur = A_and_B(2*alpha[cur_alpha_i], 2*alpha[cur_alpha_j], 2*beta_sample)
						
						A_p_k[v_i, A_hori_start:A_hori_start+3] = A_block_cur
						B_p_k[v_i, B_hori_start:B_hori_start+2] = B_block_cur
					else:
						# These frequencies have EB from dust
						A_block_cur, B_block_cur = A_and_B_primed(2*alpha[cur_alpha_i], 2*alpha[cur_alpha_j], 2*beta_sample, psi_l[k], [A_])
						
						A_p_k[v_i, A_hori_start:A_hori_start+3] = A_block_cur
						B_p_k[v_i, B_hori_start:B_hori_start+2] = B_block_cur
					v_i += 1
	return A_p_k, B_p_k

@njit
def R_mat(ti2, tj2, inv = False):
	diag = cos(ti2)*cos(tj2)
	a_diag = sin(ti2)*sin(tj2)
	if inv:
		amp = 2/(cos(2*ti2)+cos(2*tj2))
		return amp * np.array([[diag, -a_diag],[-a_diag, diag]])
	else:
		return np.array([[diag, a_diag],[a_diag, diag]])

@njit
def D_vec(ai2, aj2):
	return np.array([cos(ai2)*cos(aj2), -sin(ai2)*sin(aj2)])

@njit
def D_mat(ai2, aj2):
	diag = cos(ai2)*sin(aj2)
	a_diag = cos(aj2)*sin(ai2)

	return np.array([[-diag, -a_diag], [a_diag, diag]])

@njit
def R_vec(ti2, tj2):	
	return np.array([cos(ti2)*sin(tj2), -sin(ti2)*cos(tj2)])

@njit
def A_and_B(ai2, aj2, b2):
	R_vec_alpha = R_vec(ai2, aj2)
	R_mat_inv = R_mat(ai2, aj2, inv=True)

	R_vec_alpha_beta = R_vec(ai2 + b2[0], aj2 + b2[1])
	R_mul = np.dot(R_mat_inv, R_mat(ai2 + b2[0], aj2 + b2[1]))
	
	block_A = -np.dot(R_vec_alpha, R_mat_inv)
	A_ = np.append(block_A, [1.])
	B_ = R_vec_alpha_beta - np.dot(R_vec_alpha, R_mul)

	return A_, B_

@njit
def get_minus_Lambda_T_inv_LambdaMatrix(ai2, aj2, A_sin_4_psi):
	# Get -Lambda^T * Lambda^{-1}
	# See Eq. 12 in Planck + WMAP joint analysis paper

	F = A_sin_4_psi * np.array([[1, 0],
															[1, 0]])
	lambda_vec = R_vec(ai2, aj2) + np.dot(D_vec(ai2, aj2), F)
	lambda_mat = R_mat(ai2, aj2) + np.dot(D_mat(ai2, aj2), F)
	return - np.dot(lambda_vec, np.linalg.inv(lambda_mat))

@njit
def A_and_B_primed(ai2, aj2, b2, psi_l, A):	
	# The A and B matrices when accounting for dust EB for a given subset of maps
	# ai2 = 2*alpha_i
	# aj2 = 2*alpha_j
	# b2 = 2*beta - 2 index array, one for i and one for j. Assuming frequency-independent beta
									# then b2[0] = b2[1]

	A_sin_4_psi = A[0] * sin(4*psi_l)
	minus_Lambda_T_inv_LambdaMatrix = get_minus_Lambda_T_inv_LambdaMatrix(ai2, aj2, A_sin_4_psi) 

	# See Eq. 12 of Planck+WMAP joint analysis paper
	A_ = np.append(minus_Lambda_T_inv_LambdaMatrix, [1.])

	R_vec_beta = R_vec(ai2 + b2[0], aj2 + b2[1])
	R_mat_beta = R_mat(ai2 + b2[0], aj2 + b2[1])
	B_ = R_vec_beta + np.dot(minus_Lambda_T_inv_LambdaMatrix, R_mat_beta)

	return A_, B_

@njit
def do_not_count_this_eq(i, j, im, jm):
	# An equation should not be counted if this function returns true.
	# We do not use auto spectra, hence why we use this function.
	# For example, dont use 100Ax100A
	return im == jm and i == j

@njit
def size(nob):
	# Number of maps
	return 4*nob**2 - nob*2

@njit(parallel=True)
def likelihood_prob(c_l_o_bin_a, c_l_th_bin_a, obs_cov_bin_a, ln_det_cov, bin_range, nob, alpha, beta):
	# Calculates the log likelihood when NOT accounting for dust EB
	A, B = get_A_B_unprimed(nob, alpha, beta)
	prob = 0
	for l_b in prange(bin_range):
		c_l_o, c_l_th, obs_cov = c_l_o_bin_a[l_b, :], c_l_th_bin_a[l_b, :], obs_cov_bin_a[l_b, :, :]
		v = np.dot(A, c_l_o) - np.dot(B, c_l_th)
		C = np.dot(np.dot(A, obs_cov), np.transpose(A))

		cur = np.dot(v, np.linalg.solve(C, v))

		prob += cur
	
		if ln_det_cov:
			_, logdet = np.linalg.slogdet(C)
			prob += logdet

	return -prob/2

@njit(parallel=True)
def likelihood_prob_model_eb(c_l_o_bin_a, c_l_th_bin_a, obs_cov_bin_a, ln_det_cov, bin_range, nob, alpha, beta, psi_l, A, turn_off_for_lowest_indices):
	# Calculates the log likelihood when accounting for dust EB
	# c_l_o_bin_a = observed binned power spectra
	# c_l_th_bin_a = theoretical LCDM binned power spectra
	# obs_cov_bin_a = binned covariance matrix
	# ln_det_cov = True if the log det M should be in the log likelihood
	# bin_range = how many bins do we have?
	# nob = number of bands (this counts K and Ka as one band with different detector splits. This is taken into account when doing frequency analysis)
	# alpha = miscalibration angles
	# beta = cosmic birefringence angles
	# psi_l = misaligment angle as a function of multipole
	# A = sampled dust EB amplitudes
	# turn_off_for_lowest_indices = The low frequency bands do not have substantial dust EB. 
																# so we order the frequency channels so that the 
																# synch channels come first, then dust channels later
																# Hence, we only apply the dust EB model to maps after
																# turn_off_for_lowest_indices
	prob = 0
	for l_b in prange(bin_range):
		A_mat, B_mat = get_A_B_primed_k(nob, alpha, beta, psi_l, A, turn_off_for_lowest_indices, l_b)
		
		c_l_o, c_l_th, obs_cov = c_l_o_bin_a[l_b, :], c_l_th_bin_a[l_b, :], obs_cov_bin_a[l_b, :, :]
		v = np.dot(A_mat, c_l_o) - np.dot(B_mat, c_l_th)

		A_Cov = np.dot(A_mat, obs_cov)
		C = np.dot(A_Cov, np.transpose(A_mat))
		
		cur = np.dot(v, np.linalg.solve(C, v))

		prob += cur
	
		if ln_det_cov:
			_, logdet = np.linalg.slogdet(C)
			prob += logdet

	return -prob/2

#@njit(parallel=True)
def get_inverse_variance_mean_std(
		c_l_o_bin_a,
		c_l_th_bin_a,
		obs_cov_bin_a,
		bin_range,
		nob,
		alpha,
		beta,
		psi_l,
		A,
		turn_off_for_lowest_indices,
		with_band=False
	):
	eb_mean, eb_std = np.zeros(bin_range), np.zeros(bin_range)
	eb_alpha, eb_beta = np.zeros(bin_range), np.zeros(bin_range)
	xi_squared = 0

	for l_b in prange(bin_range):
		A_mat, B_mat = get_A_B_primed_k(nob, alpha, beta, psi_l, A, turn_off_for_lowest_indices, l_b)

		c_l_o, c_l_th, obs_cov = c_l_o_bin_a[l_b, :], c_l_th_bin_a[l_b, :], obs_cov_bin_a[l_b, :, :]
		
		if with_band:
			A_no_A_no_alpha_no_beta, _ = get_A_B_primed_k(nob, np.zeros(len(alpha)), np.array([0.0]), psi_l, np.zeros(4), turn_off_for_lowest_indices, l_b)
			_, B_no_A_no_alpha = get_A_B_primed_k(nob, np.zeros(len(alpha)), beta, psi_l, np.zeros(4), turn_off_for_lowest_indices, l_b)
			A_no_A, _ = get_A_B_primed_k(nob, alpha, beta, psi_l, np.zeros(4), turn_off_for_lowest_indices, l_b)
			
			# Contribution from beta
			v_beta = np.dot(B_no_A_no_alpha, c_l_th)
			C_beta = np.dot(np.dot(A_no_A_no_alpha_no_beta, obs_cov), np.transpose(A_no_A_no_alpha_no_beta))
			C_beta_inv = np.linalg.inv(C_beta)

			# Contribution from alpha
			A_alpha = A_no_A_no_alpha_no_beta-A_no_A
			v_alpha = np.dot(A_alpha, c_l_o)
			C_alpha = np.dot(np.dot(A_no_A, obs_cov), np.transpose(A_no_A))
			C_alpha_inv = np.linalg.inv(C_alpha)
			vector_1 = np.ones(v_beta.size)

			v1_alpha_cov = np.dot(vector_1, C_alpha_inv)
			denominator_alpha = np.dot(v1_alpha_cov, vector_1)
			eb_alpha[l_b] = np.dot(v1_alpha_cov, v_alpha) / denominator_alpha

			v1_beta_cov = np.dot(vector_1, C_beta_inv)
			denominator_beta = np.dot(v1_beta_cov, vector_1)
			eb_beta[l_b] = np.dot(v1_beta_cov, v_beta) / denominator_beta
		else:
			v = np.dot(A_mat, c_l_o) - np.dot(B_mat, c_l_th)
			C = np.dot(np.dot(A_mat, obs_cov), np.transpose(A_mat))
			C_inv = np.linalg.inv(C)
			vector_1 = np.ones(v.size)

			v1_cov = np.dot(vector_1, C_inv)
			denominator = np.dot(v1_cov, vector_1)
			eb_mean[l_b] = np.dot(v1_cov, v) / denominator
			eb_std[l_b] = np.sqrt(1 / denominator)

			xi_squared += eb_mean[l_b] ** 2 / eb_std[l_b]**2

	if with_band:
		eb_mean = eb_alpha
		eb_std = eb_beta

	return eb_mean, eb_std, xi_squared

@njit
def covariance(l, i, j, p, q, im, jm, pm, qm, ocs, f):
	# Caluclates the covariance matrix for a given set of maps
	# i,j,p,q are frequency bands
	# im, jm, pm, qm are a given time/detector split
	# ocs = observed power spectra
	# f = sky fraction

	cov = np.zeros(3)
	cov[0] = ocs[i, p, im, pm, 0, l] * ocs[j, q, jm, qm, 0, l] + ocs[i, q, im, qm, 0, l] * ocs[p, j, pm, jm, 0, l] #EE
	cov[1] = ocs[i, p, im, pm, 1, l] * ocs[j, q, jm, qm, 1, l] + ocs[i, q, im, qm, 1, l] * ocs[p, j, pm, jm, 1, l] #BB
	cov[2] = ocs[i, p, im, pm, 0, l] * ocs[j, q, jm, qm, 1, l]# + ocs[i, q, im, qm, 2, l] * ocs[p, j, pm, jm, 2, l] #EB

	cov /= (2*l+1)
	cov /= (f[i, im]*f[j, jm]*f[p, pm]*f[q, qm])**(1/4)

	return cov


@njit(parallel=True)
def get_covariance(size, l_max, nob, ocs, f):
	# Gets the full covariance matrix for each multipole
	# The binning happens later

	obs_cov = np.zeros((size*3, size*3, l_max+1))

	for l in prange(2, l_max):
		v_i = 0
		for i in range(nob):
			for j in range(nob):
				for im in range(2):
					for jm in range(2):
						if do_not_count_this_eq(i, j, im, jm):
							continue
						h_i = 0
						for p in range(nob):
							for q in range(nob):
								for pm in range(2):
									for qm in range(2):
										if do_not_count_this_eq(p, q, pm, qm):
											# Make sure we dont use auto correlations 100Ax100A for example
											continue

										obs_cur = covariance(l, i, j, p, q, im, jm, pm, qm, ocs, f)
										obs_cov[3*v_i : 3*v_i+3, 3*h_i : 3*h_i+3, l] = np.diag(obs_cur)

										h_i += 1
						v_i += 1
	return obs_cov

@njit
def get_A_ell(k, A):
	# For a given bin, what A_\ell should we use?

	num_A = len(A)
	if num_A == 1:
		A_= A[0]

	elif num_A == 2:
		# 51 < l < 511, 511 < l
		if k < 23:
			A_ = A[0]
		else:
			A_ = A[1]
			
	elif num_A == 3:
		# 51 < l < 211, 211 < l < 511, 511 < l
		if k < 8:
			A_ = A[0]
		elif k < 23:
			A_ = A[1]
		else:
			A_ = A[2]

	elif num_A == 4:
		# 51 < l < 131, 131 < l < 211, 211 < l < 511, 511 < l
		if k < 4:
			A_ = A[0]
		elif k < 8:
			A_ = A[1]
		elif k < 23:
			A_ = A[2]
		else:
			A_ = A[3]

	elif num_A == 5:
		# 51 < l < 171, 171 < l < 311, 311 < l < 471, 471 < l < 811, 811 < l
		if k < 6:
			A_ = A[0]
		elif k < 13:
			A_ = A[1]
		elif k < 21:
			A_ = A[2]
		elif k < 38:
			A_ = A[3]
		else:
			A_ = A[4]

	elif num_A == 6:
		# 51 < l < 131, 131 < l < 211, 211 < l < 311, 311 < l < 431, 431 < l < 811, 811 < l
		if k < 4:
			A_ = A[0]
		elif k < 8:
			A_ = A[1]
		elif k < 13:
			A_ = A[2]
		elif k < 19:
			A_ = A[3]
		elif k < 38:
			A_ = A[4]
		else:
			A_ = A[5]

	return A_