import sys
import os
import json
import gemmi
import glob

def get_ligand(ligand_tlc: str, pdb_file_name: str):
   st = gemmi.read_pdb(pdb_file_name)
   for model in st:
      for chain in model:
         for residue in chain:
            if residue.name == ligand_tlc:
               return residue
   return False

def get_reso(pdb_file_name):
    st = gemmi.read_pdb(pdb_file_name)
    return st.resolution

def compare_ligands(pdb_file_name: str, refmac_pdb_file_name: str, accession_code:str, ligand_tlc: str, acedrg_cif_file_name: str):

   if not os.path.exists(pdb_file_name):
      print("Bad PDB file name", pdb_file_name)
   if not os.path.exists(refmac_pdb_file_name):
      print("Bad refmac file name", refmac_pdb_file_name)
   if not os.path.exists(acedrg_cif_file_name):
      print("bad acedrg cif file name", acedrg_cif_file_name)
   if len(ligand_tlc) > 3:
      print("bad ligand tlc", ligand_tlc)
   if len(ligand_tlc) == 0:
      print("bad tlc")

   lig_1 = get_ligand(ligand_tlc,        pdb_file_name)
   lig_2 = get_ligand(ligand_tlc, refmac_pdb_file_name)
   reso = get_reso(pdb_file_name)

   for atom_1, atom_2 in zip(lig_1, lig_2):
      pos_1 = atom_1.pos
      pos_2 = atom_2.pos
      d = pos_1.dist(pos_2)
      print("delta:", ligand_tlc, accession_code, str(reso), atom_1.name, d)
      if d > 1.0:
         print(refmac_pdb_file_name)

   return False

def multi_compare_ligands(refmac_dir):
   for d in glob.glob(os.path.join(refmac_dir, "*")):
      for info_file in glob.glob(os.path.join(d, "*.info.json")):
         print(info_file)
         with open(info_file, "r") as info_file_data:
            info = json.load(info_file_data)
            print(info)
            refmac_out_fn = info['refmac_pdb_file_name']
            # refmac_out_fn may not be there because refmac didn't properly complete
            if os.path.exists(refmac_out_fn):
               compare_ligands(info['pdb_file_name'], info['refmac_pdb_file_name'], info['accession_code'], info['ligand-code'], info['acedrg_cif_file_name'])
            else:
               print("missing refmac-out", refmac_out_fn, "for ligand-code", info['ligand-code'])

def compare_ligands_in_table_file(this_weeks_ligands_table_file_name):
   with open(this_weeks_ligands_table_file_name) as f:
      lines = f.readlines()
      for line in lines:
         parts = line.split()
         if parts:
               shell_script_file_name = parts[0]
               refmac_log_file_name   = parts[1]
               fn       = os.path.basename(refmac_log_file_name)
               dir_name = os.path.dirname(refmac_log_file_name)
               parts = fn.split("-")
               p = parts[1]
               parts_2 = p.split(".")
               code = parts_2[0]
               json_file_name = code.upper() + ".info.json"
               info_file = os.path.join(dir_name, json_file_name)
               print(info_file)
               with open(info_file, "r") as info_file_data:
                  info = json.load(info_file_data)
                  print(info)
                  refmac_out_fn = info['refmac_pdb_file_name']
                  # refmac_out_fn may not be there because refmac didn't properly complete
                  if os.path.exists(refmac_out_fn):
                     compare_ligands(info['pdb_file_name'], info['refmac_pdb_file_name'], info['accession_code'], info['ligand-code'], info['acedrg_cif_file_name'])
                  else:
                     print("missing refmac-out", refmac_out_fn, "for ligand-code", info['ligand-code'])

               # compare_ligands(info['pdb_file_name'], info['refmac_pdb_file_name'], info['accession_code'], info['ligand-code'], info['acedrg_cif_file_name'])



if __name__ == "__main__":
   refmac_dir = "refmac"
   dated_refmac_logs_table_list = []
   for arg in sys.argv[1:]:
      print(arg)
      table_file = arg
      compare_ligands_in_table_file(table_file)

   # multi_compare_ligands(refmac_dir)
