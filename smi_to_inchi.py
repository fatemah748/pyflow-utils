import openbabel as ob
import sys

def getInChIAndInChIKey(smiles):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "inchi")
    ob.obErrorLog.SetOutputLevel(0)

    mol = ob.OBMol()
    conv.ReadString(mol, smiles)
    inchi = conv.WriteString(mol).rstrip()
    conv.SetOptions("K", conv.OUTOPTIONS)
    inchikey = conv.WriteString(mol).rstrip()
    return(inchi, inchikey)


def main(argv = None):
    if argv is None:
       argv = sys.argv

    if len(argv) != 4:
       sys.exit("Usage: smi_to_inchi.py recommended-smiles-file dft-labeled-inchi-file output-matched-file")

    rec_smiles = argv[1]
    inchi_labels = argv[2]

    with open(rec_smiles) as f: 
      iter_smiles = f.readlines()
    smi_inchi = {}
    for smi in iter_smiles: 
      inchi, inchikey = getInChIAndInChIKey(smi)
      smi_inchi[inchikey]=(smi.replace("\n",""), inchikey)
   
    final_labels = {} 
    inchi_dft = {}
    with open(inchi_labels) as f: 
      for line in f: 
       val= [x.strip() for x in line.split(",")] 
       inchi=val[0]
       label=val[1]
       inchi_dft[inchi] = label
    for key, value in inchi_dft.items():
      if smi_inchi[key][0] == None: 
        print(key)
      else: 
        print("dft inchi: " + key) 
        print("rec smi: " + smi_inchi[key][0].replace("\n", ""))
        print("rec smi inchi: " + smi_inchi[key][1] + "\n")   
      final_labels[smi_inchi[key][0]] = value
    with open(argv[3], "w") as w: 
      for key, value in final_labels.items(): 
      #for key, value in smi_inchi.items(): 
       # w.write("rec smi: " + value + "\n")
       # w.write("rec smi inchi: " + key + "\n") 
       
       w.write(str(key) + ", " + value + "\n")
      w.close() 
    
    return 0 
           



if __name__ == '__main__':
   main()
