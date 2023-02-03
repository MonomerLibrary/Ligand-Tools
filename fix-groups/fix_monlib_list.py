"""
Author: "Keitaro Yamashita, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, generators
import argparse
import os
import gemmi

parser = argparse.ArgumentParser(prog="fix_monlib_list",
                                 description="Make group names consistent with monomer definitions")
parser.add_argument("--monlib",
                    help="Monomer library path. Default: $CLIBD_MON")
parser.add_argument("--do_not_keep_org", action="store_true",
                    help="By default it overwrites target files and keep .org files.")

def relative_monomer_path(mon): # ref: gemmi/monlib.hpp
    if mon in ("AUX", "COM", "CON", "LPT", "PRN"):
        fname = mon + "_" + mon + ".cif"
    else:
        fname = mon + ".cif"
    
    return os.path.join(mon[0].lower(), fname)

def main(args):
    monlib_list_path = os.path.join(args.monlib, "list", "mon_lib_list.cif")
    doc = gemmi.cif.read(monlib_list_path)
    groups = {}
    comp_list = doc.find_block("comp_list")
    changed = False
    for row in comp_list.find("_chem_comp.", ["id", "group"]):
        f = os.path.join(args.monlib, relative_monomer_path(row.str(0)))
        doc2 = gemmi.cif.read(f)
        b = doc2.find_block("comp_list")
        gr = b.find_values("_chem_comp.group")[0]
        groups[row.str(0)] = gr
        if gr != row.str(1):
            changed = True
            print("change: {} {} => {}".format(row.str(0), row.str(1), gr))
            row[1] = groups[row.str(0)]
            
    link_list = doc.find_block("link_list")
    for row in comp_list.find("_chem_link.", ["id", "comp_id_1", "group_comp_1", "comp_id_2", "group_comp_2"]):
        for i in (1, 3):
            gr = groups[row.str(i)]
            if row.str(i+1) != gr:
                changed = True
                row[i+1] = gr
                print("mon_lib_list link_list changed", row[0], row[i], gr)

    mod_list = doc.find_block("mod_list")
    for row in comp_list.find("_chem_mod.", ["id", "comp_id", "group_id"]):
        gr = groups[row.str(1)]
        if row.str(2) != gr:
            changed = True
            row[2] = gr
            print("mon_lib_list mod_list changed", row[0], row[1], newgr)

    if changed:
        if not args.do_not_keep_org: os.rename(monlib_list_path, monlib_list_path + ".org")
        doc.write_file(monlib_list_path)#, style=gemmi.cif.Style.Aligned)
    else:
        print("no changes")
    
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
