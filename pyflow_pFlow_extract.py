import os
import time
import sys
from glob import glob
import pandas as pd
from tqdm import tqdm
import openbabel
from rdkit import Chem
from rdkit.Chem import Draw


# This is a script to extract the data from pyflow workflows with s0_solv and tddft extraction with cis/trans or individual molecules 
# It needs to refrence a csv file that points to the s0_solv/completed and the pm7/opt_pdbs directory for each wave of the workflow 

                     ######## IMPORTANT ########
#      Make sure the files in the directory are named as s0_solv        #
#      Still working out some erros when drawing the molecules          #
#      The data will still write but the molecule may not draw          #

def getVEE(log_file):
	with open(log_file, "r") as f:
		all_abs = {}
		for line in f: 
			if "Singlet" in line:
				vee = float(line.split()[4])
				osc = float(line.split()[8][2:])
				all_abs[vee]=osc
		print(all_abs)
		max_vee=0
		for i in range(len(all_abs)): 
			max_vee=max(all_abs, key=all_abs.get) 
			if max_vee < 4.1328 and max_vee > 1.54980: 
				max_vee=max_vee
			else: 
				all_abs[max_vee]=0
		print(all_abs)	
		# returns the key with the max value
		if max_vee == 0: 
			return 0,0
		else: 
			return max_vee, all_abs[max_vee]
	return 0, 0 

def getEnergy(log_file):
	with open(log_file, "r") as f:
		for line in f:
			line = line.strip()
			if "Sum of electronic and thermal Free Energies=" in line:
				energy = float(line.split()[-1]) * 27.2114
				return energy
	return None

def get_smiles(pdb_path, inchi_key):
	pdb_file = pdb_path + "/" + inchi_key + "_0.pdb"
	#pdb_file = flow_directory + "/" + inchi_key + "_0.pdb"

	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("pdb", "can")
	
	mol = openbabel.OBMol()
	obConversion.ReadFile(mol, pdb_file)
	smiles = obConversion.WriteString(mol).strip("\n").split("\t")[0]
	return smiles
	
######### Name csv file refrenced here #############
csv = "workflow_paths.csv" # change this! 
df = pd.read_csv(csv, index_col=None)
print("reading csv file: " + csv)

######### Name output excel file here #############
excel_writer =  pd.ExcelWriter('workflow-output.xlsx', engine='xlsxwriter')
print("new excel file: " + 'red-shifted-iter1-results.xlsx')


all_dicts = {}
n = 0
for i, row in df.iterrows():
	derivative = row['derivative']	
	s0_solv_path = row['s0_path']
	tddft_path = row['tddft_path']
	pdb_path = row['pdb_path']

	sp_tddft_source = os.path.join(tddft_path,"*_td-dft.log")
	
	log_files = glob(sp_tddft_source)
	all_dicts[derivative] = {}
	for f in log_files:
		temp_key = f.split('/')
		inchi_key = temp_key[len(temp_key) - 1].replace("_td-dft.log", "")
		short_inchi = inchi_key.split("-")[0]
		smiles = get_smiles(pdb_path, inchi_key)
		vee,osc = getVEE(f)
		if vee==0:
			print(inchi_key)
			print("file name with vee=999 : " + f)
		try:
			vee_nm = round(1239.84193/vee,0)
		except: 
			vee_nm = 0
		energy_log_file = s0_solv_path + "/" + inchi_key + "_s0_solv.log"
		try: 
			energy = getEnergy(energy_log_file)
		except FileNotFoundError: 
			energy = 0
		all_dicts[derivative][inchi_key]= [inchi_key, round(vee,2), vee_nm, round(vee,2), vee_nm, osc, energy, smiles]
		n += 1 
csv_header = ["2D Structure", "short/long_inchi", "smiles","cis_inchi", "cis_VEE(eV)","cis_VEE(nm)","cis_osc","trans_inchi", "trans_VEE(eV)","trans_VEE(nm)","trans_osc","cis-trans","reaction-energy","1 Mol","VEE(eV)", "VEE(nm)","osc","Energy (kcal/mol)"]
sheet_df = pd.DataFrame(columns=csv_header)
all_mols = {}
with tqdm(total=n, file=sys.stdout) as pbar: 
	for key_dict, value in all_dicts.items(): 
		for key, value_inchi in value.items():
			short_inchi = key.split("-")[0]
			if '/N=N/' in value_inchi[7]:
				trans=True
			else: 
				trans=False 
			if short_inchi in all_mols.keys():
				if not trans: 
					all_mols[short_inchi]['cis_inchi']= value_inchi[0]
					all_mols[short_inchi]['cis_VEE(eV)']=value_inchi[1]
					all_mols[short_inchi]['cis_VEE(nm)']=value_inchi[2]
					all_mols[short_inchi]['cis_osc']=value_inchi[5]
					all_mols[short_inchi]['trans_VEE(eV)']=all_mols[short_inchi]['VEE(eV)']
					all_mols[short_inchi]['trans_VEE(nm)']=all_mols[short_inchi]['VEE(nm)']
					all_mols[short_inchi]['trans_inchi'] = all_mols[short_inchi]['short/long_inchi']
					all_mols[short_inchi]['trans_osc']= all_mols[short_inchi]['osc']
					all_mols[short_inchi]['short/long_inchi'] = all_mols[short_inchi]['short/long_inchi'].split('_')[0]
				if trans:
					all_mols[short_inchi]['trans_inchi']= value_inchi[0]
					all_mols[short_inchi]['trans_VEE(eV)']=value_inchi[1]
					all_mols[short_inchi]['trans_VEE(nm)']=value_inchi[2]
					all_mols[short_inchi]['trans_osc']=value_inchi[5]
					all_mols[short_inchi]['cis_VEE(eV)']=all_mols[short_inchi]['VEE(eV)']
					all_mols[short_inchi]['cis_VEE(nm)']=all_mols[short_inchi]['VEE(nm)']
					all_mols[short_inchi]['cis_inchi'] = all_mols[short_inchi]['short/long_inchi']
					all_mols[short_inchi]['cis_osc']= all_mols[short_inchi]['osc']
					all_mols[short_inchi]['short/long_inchi'] = all_mols[short_inchi]['short/long_inchi'].split('_')[0]
					all_mols[short_inchi]['smiles']=value_inchi[7]
			
				if value_inchi[2] == None or all_mols[short_inchi]['VEE(nm)'] == None: 
					all_mols[short_inchi]['cis-trans'] = None
				else: 
					all_mols[short_inchi]['cis-trans'] = abs(float(value_inchi[2]) - float(all_mols[short_inchi]['VEE(nm)']))
				if value_inchi[6] == None or all_mols[short_inchi]['Energy'] == None: 
					all_mols[short_inchi]['reaction-energy'] = None
				else:
					all_mols[short_inchi]['reaction-energy']=round(23.060541945329 * (value_inchi[6] - all_mols[short_inchi]['Energy']), 3)
				all_mols[short_inchi]['1 Mol'] = 'False'
				all_mols[short_inchi]['Energy']=None 
				all_mols[short_inchi]['VEE(eV)'] = None
				all_mols[short_inchi]['VEE(nm)'] = None
				all_mols[short_inchi]['osc'] = None 
				pbar.update(1)
			else: 
				all_mols[short_inchi] = {'short/long_inchi': value_inchi[0], 'smiles' : value_inchi[7], 'cis_inchi': None, 'cis_VEE(eV)' : None,'cis_VEE(nm)' : None,'cis_osc': None, 'trans_osc': None, 'trans_inchi' : None, 'trans_VEE(eV)' : None,'trans_VEE(nm)' : None, 'cis-trans' : None,'reaction-energy': None, 'Energy' : value_inchi[6], 'VEE(eV)' : value_inchi[1],'VEE(nm)' : value_inchi[2],'1 Mol': 'True','osc': value_inchi[5]}
				
				pbar.update(1)

for key, value in all_mols.items():
	mol_row = {}
	for h in csv_header:
		mol_row[h] = None
	mol_row["short/long_inchi"] = value['short/long_inchi']
	mol_row["smiles"] = value['smiles']
	mol_row["cis_inchi"] = value['cis_inchi']
	mol_row["cis_VEE(eV)"] = value['cis_VEE(eV)']
	mol_row["cis_VEE(nm)"] = value['cis_VEE(nm)']
	mol_row["cis_osc"]=value['cis_osc']
	mol_row["trans_inchi"] = value['trans_inchi']
	mol_row["trans_VEE(eV)"] = value['trans_VEE(eV)']
	mol_row["trans_VEE(nm)"] = value['trans_VEE(nm)']
	mol_row["trans_osc"]=value['trans_osc']
	mol_row["cis-trans"] = value['cis-trans']
	mol_row["reaction-energy"] = value['reaction-energy']
	mol_row["1 Mol"] = value["1 Mol"]
	mol_row["VEE(eV)"] = value['VEE(eV)']
	mol_row["VEE(nm)"] = value['VEE(nm)']
	mol_row["osc"] = value['osc']
	mol_row["Energy (kcal/mol)"] = value['Energy']
	
	sheet_df = sheet_df.append(mol_row, ignore_index=True)

sheet_df.to_excel(excel_writer,sheet_name='all_mols', index=False)
worksheet = excel_writer.sheets['all_mols']
worksheet.set_default_row(110)
worksheet.set_column('A:A', 20.7)
worksheet.set_column('B:K', 10.0)

row = 2 
for j, mol in sheet_df.iterrows(): 
	smiles = mol["smiles"]
	if mol["1 Mol"] == 'True':
		inchi_key = mol["short/long_inchi"]
	else: 
		inchi_key = mol["trans_inchi"]
	mol_obj = Chem.MolFromSmiles(smiles)
	image_file = os.path.join("images", inchi_key + ".png")
	try: 
		Draw.MolToFile(mol_obj, image_file, size=(300,300))
	except Exception as e: 
		print(e)
	worksheet.insert_image("A{}".format(row), image_file, {'x_scale': 0.5, 'y_scale': 0.5})
	row += 1
excel_writer.save()  

