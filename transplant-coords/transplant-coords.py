#!/usr/bin/env python

############################
## transplant-coords.py
## 2022
## Author: Rob Nicholls
############################

import argparse
import datetime
import sys
import os
import gemmi

def transplantCoords(cif_dict,new_coords,opt={}):
  if not 'output_dir' in opt:
    opt['output_dir'] = "new_cifs"

  try:
      # Read PDB file containing coordinates, and create dict containing coordinates.
      if not os.path.isfile(new_coords):
        raise Exception("Cannot find file"+new_coords)
      st = gemmi.read_structure(new_coords)
      if len(st) != 1:
        raise Exception("Number of models: "+str(len(st)))
      coord_dict = {}
      for chain in st[0]:
        for res in chain:
          if res.name in coord_dict:
            raise Exception("Multiple residues found with same component ID in new coordinate file.")
          coord_dict[res.name] = {}
          for atom in res:
            coord_dict[res.name][atom.name] = atom.pos

      # Read CIF file and transplant coordinates.
      if not os.path.isfile(cif_dict):
        raise Exception("Cannot find file"+cif_dict)
      doc = gemmi.cif.read_file(cif_dict)
      comp_list = doc.find_block("comp_list")
      category = comp_list.find_mmcif_category("_chem_comp")
      comp_id_list = []
      for row in category:
        comp_id_list.append(row["id"])
        
      for comp_id in comp_id_list:
        print("Processing comp_id: "+comp_id)
        if not comp_id in coord_dict:
          raise Exception("Component ID "+comp_id+" not found in CIF dictionary.")
        block = doc.find_block("comp_"+comp_id)
        category = block.find_mmcif_category("_chem_comp_atom")
        atom_list = []
        for row in category:
          atom_list.append(row["atom_id"])
        if set(atom_list) != set(coord_dict[comp_id].keys()):
          print("Unique atoms in CIF dictionary: "+set(atom_list).difference(set(coord_dict[comp_id].keys())))
          print("Unique atoms in coordinate file: "+(set(set(coord_dict[comp_id].keys())).difference(set(atom_list))))
          raise Exception("Input CIF and coordinate file do not contain set of atoms.")
        for row in category:
          if not row["atom_id"] in coord_dict[comp_id]:
            raise Exception("Atom "+row["atom_id"]+" not found in new coordinates.")
          pos = coord_dict[comp_id][row["atom_id"]]
          #print(row["comp_id"],row["atom_id"],pos)
          row["x"] = str(pos.x)
          row["y"] = str(pos.y)
          row["z"] = str(pos.z)

      dir_new = opt["output_dir"]
      if not os.path.isdir(dir_new):
        os.mkdir(dir_new)
      cif_dict_new = os.path.join(dir_new,os.path.basename(cif_dict))
      doc.write_file(cif_dict_new)
      print("File written to: "+cif_dict_new)
      
  except Exception as e:
    print("Error:",e)
    pass

if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('cif_dict', help="Input CIF dictionary")
   parser.add_argument('coord_file', help="PDB/mmCIF containing coordinates to be transplanted")
   parser.add_argument('-o','--output_dir', dest='output_dir', required=False, help="Output directory for the new CIFs", default='new_cifs')
   args = parser.parse_args()
   
   opt = {}
   opt['output_dir'] = args.output_dir
      
   time_start = datetime.datetime.now()
   result = transplantCoords(args.cif_dict,args.coord_file,opt)
   time_end = datetime.datetime.now()

   print("Job completed ("+str(round((time_end-time_start).total_seconds(),3))+" seconds)")
