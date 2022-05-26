import subprocess
import sys
sys.path.append('../parameter_files/')
from hfi_lfi_wmap_eb import maps_param


f = open('param_skeleton.txt', 'r+')
#f = open('param_experimental_skeleton.txt', 'r+')

mask_root = 'masks/combined_ps_co_{}_{}.fits'

cl_out_folder = 'cl_wmap_hfi_lfi/'


mask_list = ['0', '30']

sc = f.read()

AB_list = ['A', 'B']

def modify_text(ot, string_pos, insert_string):
	return ot[:ot.find(string_pos)+len(string_pos)] + insert_string + ot[ot.find(string_pos)+len(string_pos):]

for h, mask in enumerate(mask_list):
	print('Doing mask percent', h)
	i = 0
	
	for freq1, freq1_param in maps_param.items():
		print(freq1_param)
		for k in range(2):
			weightfile = mask_root.format(mask, str(freq1_param['n_side']))
			
			if freq1_param['data_set'] == 'npipe':
				mapfile =  freq1_param['map_root'].format(freq1_param['data_split'][k], freq1_param['data_split'][k])
			elif freq1_param['data_set'] == 'wmap' or freq1_param['data_set'] == 'bp':
				mapfile =  freq1_param['map_root'].format(freq1_param['data_split'][k])
			j = 0
			for freq2, freq2_param in maps_param.items():
				for l in range(2):
					if (j*2 + l < i*2 + k):
						#Only get upper triangle
						continue
					weightfile2 = mask_root.format(mask, str(freq2_param['n_side']))
					
					if freq2_param['data_set'] == 'npipe':
						mapfile2 =  freq2_param['map_root'].format(freq2_param['data_split'][l], freq2_param['data_split'][l])
					elif freq2_param['data_set'] == 'wmap' or freq2_param['data_set'] == 'bp':
						mapfile2 =  freq2_param['map_root'].format(freq2_param['data_split'][l])
					
					print(weightfile, weightfile2, mapfile, mapfile2)

					clout = '../' + cl_out_folder + 'mask_percent_{}/{}{}_{}{}.dat'.format(mask, freq1.replace("_", ""), AB_list[k], freq2.replace("_", ""), AB_list[l])

					new_file = modify_text(sc, 'clfile = ', clout)
					new_file = modify_text(new_file, 'mapfile = ', mapfile)
					new_file = modify_text(new_file, 'mapfile2 = ', mapfile2)
					new_file = modify_text(new_file, 'weightfile = ', weightfile)
					new_file = modify_text(new_file, 'weightfilep = ', weightfile)
					new_file = modify_text(new_file, 'weightfile2 = ', weightfile2)
					new_file = modify_text(new_file, 'weightfilep2 = ', weightfile2)
						
					cf = open("param_gen.txt", "w")
					cf.write(new_file)
					cf.close()

					bashCommand = '/path/to/polspice/PolSpice_v03-07-02/build/spice -optinfile param_gen.txt'
					print('Running', i, j, k, l)
					process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
					for line in process.stdout:
	    					print(line)
					process.wait()
					print(process.returncode)
				j += 1
		i += 1
		#print(skeleton_content)