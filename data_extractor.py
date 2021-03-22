import os
import time
import sys
from glob import glob
import pandas as pd
from tqdm import tqdm
import openbabel
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)

def getVEE(log_file, dir):
    with open(log_file, "r") as f:
        all_abs = {}
        for line in f:
            if "Singlet" in line:
                vee = float(line.split()[4])
                osc = float(line.split()[8][2:])
                all_abs[vee] = osc
        max_vee = 0
        for i in range(len(all_abs)):
            max_vee = max(all_abs, key=all_abs.get)
            if max_vee < 4.1328 and max_vee > 1.54980:
                max_vee = max_vee
            else:
                all_abs[max_vee] = 0
        # returns the key with the max value
        if max_vee == 0:
            return 0, 0
        else:
            log_file_name = log_file.split("/")
            log_file_name = log_file_name[len(log_file_name) - 1].replace("_" + dir + ".log", "")
            return log_file_name, max_vee
    return 0, 0


def getEnergy(log_file, dir):
    with open(log_file, "r") as f:
        for line in f:
            line = line.strip()
            if "Sum of electronic and thermal Free Energies=" in line:
                energy = float(line.split()[-1]) * 27.2114
                log_file_name = log_file.split("/")
                log_file_name = log_file_name[len(log_file_name) - 1].replace("_" + dir + ".log", "")
                return log_file_name, energy
    return None


def get_smiles(log_files, dir, smiles):
    for file in glob(log_files):
        f = file.split("/")
        name = f[1].replace("_" + dir + ".log", "")
        if smiles.get(name) == None:
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("g98", "can")

            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, file)
            smile = obConversion.WriteString(mol).strip("\n").split("\t")[0]
            smiles[name] = smile
    return smiles

def write_to_csv(file_name, data, smiles, dirs):
    excel_writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
    csv_header = ["visual", "file_name", "smiles"]
    sheet_df = pd.DataFrame(columns=csv_header)
    for key, value in dirs.items():
        csv_header.append(key)
    for key, value in data.items():
        mol_row = {}
        for h in csv_header:
            mol_row[h] = None
        mol_row["file_name"] = key
        mol_row["smiles"] = smiles[key]
        # assigns each energy to the inchi
        for key2, value2 in value.items():
            mol_row[key2] = value2
        sheet_df = sheet_df.append(mol_row, ignore_index=True)
    sheet_df.to_excel(excel_writer, sheet_name='all_mols', index=False)
    worksheet = excel_writer.sheets['all_mols']
    worksheet.set_default_row(110)
    worksheet.set_column('A:A', 20.7)
    worksheet.set_column('B:K', 10.0)
    row = 2
    for j, mol in sheet_df.iterrows():
        smiles = mol["smiles"]
        mol_obj = Chem.MolFromSmiles(smiles, sanitize=False)
        Chem.SanitizeMol(mol_obj, sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP)
        image_file = os.path.join("images", mol["file_name"] + ".png")
        try:
            Draw.DrawingOptions.atomLabelFontSize = 55
            Draw.DrawingOptions.dotsPerAngstrom = 200
            Draw.DrawingOptions.bondLineWidth = 10.0
            Draw.DrawingOptions.colorBonds = False
            Draw.MolToFile(mol_obj, image_file, size=(300, 300))
        except Exception as e:
            print(e)

        worksheet.insert_image("A{}".format(row), image_file,
                               {'x_scale': 0.5, 'y_scale': 0.5})
        row += 1
    excel_writer.save()

def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) < 1:
        sys.exit("Usage: data_extractor.py excel-name.xlsx <dir1> <dir2> <dir3> etc")

    dirs = {}
    for i in range(2, len(argv)):
        dirs[argv[i]] = {}
# assigns files with their energies to the larger directory
    all_smiles = {}
    for key, value in dirs.items():
        file_dump = os.path.join(key, "*.log")
        all_smiles = get_smiles(file_dump, key, all_smiles)
        if 'tddft' in key:
            for file in glob(file_dump):
                file_name, vee = getVEE(file, key)
                dirs[key][file_name] = vee
        else:
            for file in glob(file_dump):
                file_name, energy = getEnergy(file, key)
                dirs[key][file_name] = energy

# finds the same files and combines the values {'inchi' : {'dir1': value}, ... }
    flow_data = {}
    for key, value in dirs.items():
        print(key)
        for key2, value2 in value.items():
            print("key2: " + key2)
            flow_data[key2] = {}
            flow_data[key2][key] = value2
    print(flow_data)
    write_to_csv(argv[1], flow_data, all_smiles, dirs)

if __name__ == '__main__':
   main()