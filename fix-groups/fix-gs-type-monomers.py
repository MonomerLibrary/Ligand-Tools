"""
Author: "Keitaro Yamashita, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
# CCD:
# the current GS: OP1 (double), S2P-HSP2, OP3-HOP3
#             AS: OP1-HOP1, S2P (double), OP3-HOP3
# Monlib:
# the current GS: OP1 (double), S2P, OP3(-1)
#             AS: OP1(-1), S2P (double), OP3(-1)

from __future__ import absolute_import, division, print_function, generators
import argparse
import os
import glob
import json
import gemmi
import groups
from networkx.algorithms import isomorphism
import networkx

parser = argparse.ArgumentParser(prog="",
                                 description="Find GS-type wrong monomers (P-SH) and make them acedrg-ready cif files")
parser.add_argument("targets", nargs='+',
                    help="files or directories")
parser.add_argument("--monlib",
                    help="Monomer library path. Default: $CLIBD_MON")

def find_files(targets):
    ret = []
    for t in targets:
        if os.path.isdir(t):
            ret.extend(glob.glob(os.path.join(t, "*.cif"))) # want .gz also?
        elif os.path.isfile(t):
            ret.append(t)
        else:
            print("Error:", t, "is not file or directory")
    return ret

def graph_from_chemcomp(cc, keep_atoms=None):
    G = networkx.Graph()
    for atom in cc.atoms:
        if keep_atoms is not None and atom.id not in keep_atoms: continue
        # chem_type is not good; CA's CH1 or CH2 is not a problemm.
        G.add_node(atom.id, Z=str(atom.el.atomic_number))
    for bond in cc.rt.bonds:
        if keep_atoms is not None and bond.id1.atom not in keep_atoms: continue
        if keep_atoms is not None and bond.id2.atom not in keep_atoms: continue
        G.add_edge(bond.id1.atom, bond.id2.atom)
    return G

def make_ref_org(monlib_path):
    mon = "GS"
    exatoms = ["C3'", "C4'", "HO3'", "O3'", "O5'", 'OP1', 'S2P', 'HSP2', 'OP3', 'P', "C5'", "C1'", "C2'", "O4'"]
    f = os.path.join(monlib_path, mon[0].lower(), mon+".cif")
    cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(f)[-1])
    g = graph_from_chemcomp(cc, exatoms)
    #print(g, list(g), list(g.edges))
    return g
def make_ref(): #monlib_path):
    #mon = "GS"
    exatoms = ["C3'", "C4'", "HO3'", "O3'", "O5'", 'OP1', 'S2P', 'HSP2', 'OP3', 'P', "C5'", "C1'", "C2'", "O4'"]
    #f = os.path.join(monlib_path, mon[0].lower(), mon+".cif")
    #cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(f)[-1])
    s = """\
data_GS
_chem_comp.id                                    GS 
loop_
_chem_comp_atom.comp_id 
_chem_comp_atom.atom_id 
_chem_comp_atom.alt_atom_id 
_chem_comp_atom.type_symbol 
_chem_comp_atom.charge 
_chem_comp_atom.pdbx_align 
_chem_comp_atom.pdbx_aromatic_flag 
_chem_comp_atom.pdbx_leaving_atom_flag 
_chem_comp_atom.pdbx_stereo_config 
_chem_comp_atom.model_Cartn_x 
_chem_comp_atom.model_Cartn_y 
_chem_comp_atom.model_Cartn_z 
_chem_comp_atom.pdbx_model_Cartn_x_ideal 
_chem_comp_atom.pdbx_model_Cartn_y_ideal 
_chem_comp_atom.pdbx_model_Cartn_z_ideal 
_chem_comp_atom.pdbx_component_atom_id 
_chem_comp_atom.pdbx_component_comp_id 
_chem_comp_atom.pdbx_ordinal 
GS P      P    P 0 1 N N S -8.407 3.626  -5.514 -0.429 -0.556 4.767  P      GS 1  
GS OP1    O1P  O 0 1 N N N -9.521 3.327  -6.441 0.589  0.042  5.659  OP1    GS 2  
GS S2P    S2P  S 0 1 N N N -8.410 2.828  -3.539 -1.980 0.846  4.423  S2P    GS 3  
GS OP3    O3P  O 0 1 N Y N -8.276 5.236  -5.418 -1.041 -1.874 5.461  OP3    GS 4  
GS "O5'"  O5*  O 0 1 N N N -7.032 3.169  -6.227 0.253  -0.956 3.365  "O5'"  GS 5  
GS "C5'"  C5*  C 0 1 N N N -6.863 3.264  -7.630 0.775  0.248  2.801  "C5'"  GS 6  
GS "C4'"  C4*  C 0 1 N N R -5.401 3.011  -8.006 1.441  -0.063 1.460  "C4'"  GS 7  
GS "O4'"  O4*  O 0 1 N N N -4.539 3.822  -7.231 0.464  -0.510 0.493  "O4'"  GS 8  
GS "C3'"  C3*  C 0 1 N N S -4.939 1.577  -7.760 2.018  1.227  0.824  "C3'"  GS 9  
GS "O3'"  O3*  O 0 1 N N N -5.327 0.719  -8.822 3.336  1.491  1.309  "O3'"  GS 10 
GS "C2'"  C2*  C 0 1 N N N -3.433 1.797  -7.731 2.045  0.863  -0.680 "C2'"  GS 11 
GS "C1'"  C1*  C 0 1 N N R -3.254 3.230  -7.232 1.033  -0.291 -0.808 "C1'"  GS 12 
GS N9     N9   N 0 1 Y N N -2.661 3.263  -5.874 -0.020 0.072  -1.758 N9     GS 13 
GS C8     C8   C 0 1 Y N N -3.287 3.249  -4.652 -1.191 0.706  -1.460 C8     GS 14 
GS N7     N7   N 0 1 Y N N -2.482 3.393  -3.638 -1.892 0.871  -2.544 N7     GS 15 
GS C5     C5   C 0 1 Y N N -1.222 3.455  -4.223 -1.214 0.354  -3.598 C5     GS 16 
GS C6     C6   C 0 1 N N N 0.064  3.556  -3.613 -1.479 0.251  -4.984 C6     GS 17 
GS O6     O6   O 0 1 N N N 0.324  3.708  -2.425 -2.515 0.686  -5.456 O6     GS 18 
GS N1     N1   N 0 1 N N N 1.098  3.444  -4.540 -0.553 -0.338 -5.772 N1     GS 19 
GS C2     C2   C 0 1 N N N 0.917  3.288  -5.901 0.599  -0.826 -5.232 C2     GS 20 
GS N2     N2   N 0 1 N N N 2.015  3.174  -6.655 1.519  -1.424 -6.056 N2     GS 21 
GS N3     N3   N 0 1 N N N -0.294 3.233  -6.479 0.857  -0.736 -3.947 N3     GS 22 
GS C4     C4   C 0 1 Y N N -1.318 3.315  -5.585 -0.010 -0.161 -3.106 C4     GS 23 
GS HSP2   2HSP H 0 0 N N N -7.694 3.020  -2.943 -2.775 0.122  3.615  HSP2   GS 24 
GS HOP3   3HOP H 0 0 N N N -7.560 5.428  -4.822 -1.439 -1.591 6.295  HOP3   GS 25 
GS "H5'"  1H5* H 0 1 N N N -7.142 4.262  -7.970 1.511  0.678  3.481  "H5'"  GS 26 
GS "H5''" 2H5* H 0 0 N N N -7.496 2.527  -8.126 -0.036 0.959  2.647  "H5''" GS 27 
GS "H4'"  H4*  H 0 1 N N N -5.251 3.252  -9.060 2.224  -0.811 1.586  "H4'"  GS 28 
GS "H3'"  H3*  H 0 1 N N N -5.295 1.233  -6.788 1.361  2.078  1.006  "H3'"  GS 29 
GS "HO3'" *HO3 H 0 0 N Y N -5.039 -0.173 -8.668 3.671  2.248  0.809  "HO3'" GS 30 
GS "H2'"  1H2* H 0 1 N N N -2.952 1.083  -7.076 3.041  0.532  -0.974 "H2'"  GS 31 
GS "H2''" 2H2* H 0 0 N N N -3.022 1.724  -8.738 1.731  1.714  -1.284 "H2''" GS 32 
GS "H1'"  H1*  H 0 1 N N N -2.615 3.775  -7.927 1.542  -1.194 -1.146 "H1'"  GS 33 
GS H8     H8   H 0 1 N N N -4.353 3.140  -4.531 -1.492 1.024  -0.472 H8     GS 34 
GS HN1    HN1  H 0 1 N N N 2.040  3.473  -4.168 -0.712 -0.420 -6.725 HN1    GS 35 
GS HN21   1HN2 H 0 0 N N N 1.927  3.057  -7.654 2.345  -1.773 -5.688 HN21   GS 36 
GS HN22   2HN2 H 0 0 N N N 2.929  3.185  -6.222 1.341  -1.502 -7.007 HN22   GS 37 
# 
loop_
_chem_comp_bond.comp_id 
_chem_comp_bond.atom_id_1 
_chem_comp_bond.atom_id_2 
_chem_comp_bond.value_order 
_chem_comp_bond.pdbx_aromatic_flag 
_chem_comp_bond.pdbx_stereo_config 
_chem_comp_bond.pdbx_ordinal 
GS P     OP1    DOUB N N 1  
GS P     S2P    SING N N 2  
GS P     OP3    SING N N 3  
GS P     "O5'"  SING N N 4  
GS S2P   HSP2   SING N N 5  
GS OP3   HOP3   SING N N 6  
GS "O5'" "C5'"  SING N N 7  
GS "C5'" "C4'"  SING N N 8  
GS "C5'" "H5'"  SING N N 9  
GS "C5'" "H5''" SING N N 10 
GS "C4'" "O4'"  SING N N 11 
GS "C4'" "C3'"  SING N N 12 
GS "C4'" "H4'"  SING N N 13 
GS "O4'" "C1'"  SING N N 14 
GS "C3'" "O3'"  SING N N 15 
GS "C3'" "C2'"  SING N N 16 
GS "C3'" "H3'"  SING N N 17 
GS "O3'" "HO3'" SING N N 18 
GS "C2'" "C1'"  SING N N 19 
GS "C2'" "H2'"  SING N N 20 
GS "C2'" "H2''" SING N N 21 
GS "C1'" N9     SING N N 22 
GS "C1'" "H1'"  SING N N 23 
GS N9    C8     SING Y N 24 
GS N9    C4     SING Y N 25 
GS C8    N7     DOUB Y N 26 
GS C8    H8     SING N N 27 
GS N7    C5     SING Y N 28 
GS C5    C6     SING N N 29 
GS C5    C4     DOUB Y N 30 
GS C6    O6     DOUB N N 31 
GS C6    N1     SING N N 32 
GS N1    C2     SING N N 33 
GS N1    HN1    SING N N 34 
GS C2    N2     SING N N 35 
GS C2    N3     DOUB N N 36 
GS N2    HN21   SING N N 37 
GS N2    HN22   SING N N 38 
GS N3    C4     SING N N 39 
"""
    cc = gemmi.make_chemcomp_from_block(gemmi.cif.read_string(s)[-1])
    g = graph_from_chemcomp(cc, exatoms)
    #print(g, list(g), list(g.edges))
    return g

def main(args):
    if args.monlib is None: args.monlib = os.environ["CLIBD_MON"]
    #ref = make_ref_org(args.monlib)    
    ref = make_ref()
    #print([(x, networkx.get_node_attributes(ref, x)) for x in ref])
    #print(networkx.get_node_attributes(ref, "Z"))
    files = find_files(args.targets)
    ret = []
    for f in files:
        #print("processing", f)
        doc = gemmi.cif.read(f)
        if not doc[-1].name.startswith("comp_"):
            print("WARNING: {} is not a monomer dictionary".format(f))
            continue

        cc = groups.read_chemcomp(doc)
        g = graph_from_chemcomp(cc)
        node_match = isomorphism.categorical_node_match('Z', 0)
        GM = isomorphism.GraphMatcher(g, ref, node_match=node_match)
        if GM.subgraph_is_isomorphic():
            for mapping in GM.subgraph_isomorphisms_iter():
                rmapping = {mapping[k]: k for k in mapping}
                if all(x == rmapping["P"] or cc.find_atom(x).is_hydrogen() for a in ("OP1", "S2P", "OP3") for x in g.neighbors(rmapping[a])):
                    # Need fix
                    print(f)
                    tab = doc[-1].find("_chem_comp_atom.", ["atom_id", "charge"])
                    todel = None
                    for i, row in enumerate(tab):
                        if row.str(0) == rmapping["OP1"]:
                            row[1] = "-1"
                        elif row.str(0) == rmapping["S2P"]:
                            row[1] = "0"
                        elif row.str(0) == rmapping["OP3"]:
                            row[1] = "-1"
                        elif row.str(0) == rmapping["HSP2"]:
                            todel = i
                    del tab[todel]
                    tab = doc[-1].find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "type"])
                    todel = None
                    for i, row in enumerate(tab):
                        match_bond = lambda x, y: (row.str(0), row.str(1)) in ((x, y), (y, x))
                        if match_bond(rmapping["OP1"], rmapping["P"]):
                            row[2] = "SINGLE"
                        elif match_bond(rmapping["S2P"], rmapping["P"]):
                            row[2] = "DOUBLE"
                        elif match_bond(rmapping["OP3"], rmapping["P"]):
                            row[2] = "SINGLE"
                        elif match_bond(rmapping["S2P"], rmapping["HSP2"]):
                            todel = i
                    del tab[todel]
                    doc.write_file(f"{cc.name}.cif")

                #for a in ("OP1", "S2P", "OP3"):
                #    print(a, rmapping[a], all(x == rmapping["P"] or cc.find_atom(x).is_hydrogen() for x in g.neighbors(rmapping[a])))
                break

            #print(f)
# main()

# the current GS: OP1 (double), S2P, OP3(-1)
#             AS: OP1(-1), S2P (double), OP3(-1)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
