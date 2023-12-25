"""
Author: "Keitaro Yamashita, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, generators
import gemmi
import re
import os
import networkx
from networkx.algorithms import isomorphism

# gamma/delta linking peptides. they follow standard peptide nomenclature, but should remain non-polymer
keep_nonpolymer = ['ACB', 'FGA', 'GGL', 'IAS']
re_OPn_or_OnP = re.compile("^O(?:P([123])|([123])P)$")
non_supported_groups = (gemmi.ChemComp.Group.Pyranose,
                        gemmi.ChemComp.Group.Ketopyranose,
                        gemmi.ChemComp.Group.Furanose)
group_str = gemmi.ChemComp.group_str

def graph_from_chemcomp(cc, keep_atoms=None):
    G = networkx.Graph()
    for atom in cc.atoms:
        if keep_atoms is not None and atom.id not in keep_atoms: continue
        # chem_type is not good; CA's CH1 or CH2 is not a problemm.
        G.add_node(atom.id, Z="{}{}".format(atom.el.atomic_number, atom.charge))
    for bond in cc.rt.bonds:
        if keep_atoms is not None and bond.id1.atom not in keep_atoms: continue
        if keep_atoms is not None and bond.id2.atom not in keep_atoms: continue
        G.add_edge(bond.id1.atom, bond.id2.atom)
    return G

def read_chemcomp(doc):
    assert doc[-1].name.startswith("comp_")
    cc = gemmi.make_chemcomp_from_block(doc[-1])
    comp_list = doc.find_block("comp_list")
    if comp_list:
        tab = comp_list.find("_chem_comp.", ["id", "group"])
        if tab:
            cc.set_group(tab.find_row(cc.name).str(1))
    return cc

def prepare_references(monlib_path):
    defs = [("peptide", "GLY", ['N', 'CA', 'C', 'O', 'OXT', 'H', 'H2', 'H3']),
            ("P-peptide", "PRO", ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OXT', 'H']), # CB and CG are not needed in link
            ("M-peptide", "5JP", ['O', 'C', 'CA', 'N', 'CN', 'OXT', 'H']),
            ("DNA/RNA", "N", ["C3'", "C4'", "HO3'", "O3'", "O5'", 'OP1', 'OP2', 'OP3', 'P',
                              "C5'", "C1'", "C2'", "O4'"]), # not needed for link, but needed to be "DNA/RNA"
    ]
    groups = []
    for gr, mon, exatoms in defs:
        f = os.path.join(monlib_path, mon[0].lower(), mon+".cif")
        cc = gemmi.make_chemcomp_from_block(gemmi.cif.read(f)[-1])
        g = graph_from_chemcomp(cc, exatoms)
        groups.append((gr, g))
    return groups
# prepare_references()

def check_group(G1, G2, gr, cc):
    # ref: https://gemmi.readthedocs.io/en/latest/analysis.html#graph-isomorphism
    node_match = isomorphism.categorical_node_match('Z', 0)
    GM = isomorphism.GraphMatcher(G1, G2, node_match=node_match)
    ret = []
    genuine = False
    if GM.subgraph_is_isomorphic():
        short_diff = None
        for n, mapping in enumerate(GM.subgraph_isomorphisms_iter()):
            diff = {k: v for k, v in mapping.items() if k != v}
            if short_diff is None or len(diff) < len(short_diff):
                short_diff = diff
            if n == 10000:  # don't spend too much here
                print(' (it may not be the simplest isomorphism)')
                break
        for id1, id2 in short_diff.items():
            if gr == "DNA/RNA" and id2 in ("C5'", "C1'", "C2'", "O4'"): continue
            ret.append([id1, id2])
        if not ret:
            genuine = True

    if ret or genuine:
        if "peptide" in gr:
            N, chemtype = check_chemtype("N", cc, ret)
            if not chemtype.startswith("NT"):
                print(ret)
                print("INFO: {} is {} like but N ({}) is not NT* ({})".format(cc.name, gr, N, chemtype))
                return False, []
        if "P-peptide" == gr:
            # check all carbon is sp3
            for a in 'CA', 'CB', 'CG', 'CD':
                a2, chemtype = check_chemtype(a, cc, ret)
                if not chemtype.startswith("CH") and chemtype != "CT":
                    print(ret)
                    print("INFO: {} is P-peptide like but {} ({}) is not sp3 ({})".format(cc.name, a, a2, chemtype))
                    return False, []

            ret = [[id1, id2] for id1, id2 in ret if id2 not in ("CB", "CG")]

        if "M-peptide" == gr:
            # check if not forming a ring (P-peptide like)
            cn, cn_type = check_chemtype("CN", cc, ret)
            n = get_atom_name("N", ret)
            # PRO type ring will be 5. We can allow >6
            paths = [len(path) for path in networkx.all_simple_paths(G1, cn, n) if 2 < len(path) < 7]
            if paths:
                print(ret)
                print("INFO: {} is M-peptide like but forming {}-membered ring".format(cc.name, paths[0]))
                return False, []
            if not cn_type.startswith("CH") and cn_type != "CT":
                print(ret)
                print("INFO: {} is M-peptide like but CN ({}) is not sp3 ({})".format(cc.name, cn, cn_type))
                return False, []

        if "DNA/RNA" == gr:
            O3p, chemtype = check_chemtype("O3'", cc, ret)
            if chemtype != "OH1":
                print(ret)
                print("INFO: {} is {} like but O3' ({}) is not OH1 ({})".format(cc.name, gr, O3p, chemtype))
                return False, []

    if genuine:
        print("{} is genuinely {}".format(cc.name, gr))
            
    return genuine, ret
# check_group

def get_atom_name(name, alias):
    a = [a1 for a1, a2 in alias if a2 == name]
    return a[0] if a else name

def check_chemtype(name, cc, alias):
    a = get_atom_name(name, alias)
    chemtype = [x.chem_type for x in cc.atoms if x.id == a][0]
    return a, chemtype

def fix_OP3(cc, alias, doc):
    OP3 = get_atom_name("OP3", alias)
    P = get_atom_name("P", alias)
    OP3_charge = [x.charge for x in cc.atoms if x.id == OP3][0]
    P_OP3 = sorted((P, OP3))
    btype = [b.type for b in cc.rt.bonds if sorted((b.id1.atom, b.id2.atom)) == P_OP3][0]

    if OP3_charge != -1 or btype != gemmi.BondType.Single:
        print("WARNING {} OP3 has wrong charge and/or wrong valence ({}, {})".format(cc.name, OP3_charge, btype))
        OP1 = get_atom_name("OP1", alias)
        OP1_charge = [x.charge for x in cc.atoms if x.id == OP1][0]
        P_OP1 = sorted((P, OP1))
        btype1 = [b.type for b in cc.rt.bonds if sorted((b.id1.atom, b.id2.atom)) == P_OP1][0]
        if btype1 == gemmi.BondType.Single:
            print(" Swap OP1 and OP3 in cif")
            s1, s2 = OP3, OP1
        else:
            print(" Swap OP2 and OP3 in cif")
            OP1 = get_atom_name("OP2", alias)
            s1, s2 = OP3, OP2

        # FIXME: problematic if OP1/OP2/OP3 name is used other than atom names
        block = doc.find_block("comp_"+cc.name)
        rep = lambda v,x,y: y if v==x else v
        w = "*****"
        for item in block:
            if not item.loop: continue
            assert all(w != v for v in item.loop.values)
            vv = [[gemmi.cif.quote(rep(rep(rep(gemmi.cif.as_string(item.loop[j,i]),s1, w), s2, s1), w, s2))
                   for j in range(item.loop.length())]
                  for i in range(item.loop.width())]
            item.loop.set_all_values(vv)
        return True
# fix_OP3()

def fix_OXT(cc, alias, doc):
    OXT = get_atom_name("OXT", alias)
    O = get_atom_name("O", alias)
    C = get_atom_name("C", alias)
    OXT_charge = [x.charge for x in cc.atoms if x.id == OXT][0]
    C_OXT = sorted((C, OXT))
    btype = [b.type for b in cc.rt.bonds if sorted((b.id1.atom, b.id2.atom)) == C_OXT][0]
    
    if OXT_charge != -1 or btype != gemmi.BondType.Single:
        print("WARNING {} OXT has wrong charge and/or wrong valence ({}, {})".format(cc.name, OXT_charge, btype))
        s1, s2 = OXT, O
        block = doc.find_block("comp_"+cc.name)
        rep = lambda v,x,y: y if v==x else v
        w = "*****"
        for row in block.find("_chem_comp_atom.", ["atom_id"]):
            row[0] = gemmi.cif.quote(rep(rep(rep(row.str(0), s1, w), s2, s1), w, s2))
        for item in block:
            if not item.loop: continue
            # _chem_comp_atom.type_symbol etc has "O" that should not be replaced
            if "_chem_comp_atom.atom_id" in item.loop.tags: continue
            assert all(w != v for v in item.loop.values)
            vv = [[gemmi.cif.quote(rep(rep(rep(gemmi.cif.as_string(item.loop[j,i]),s1, w), s2, s1), w, s2))
                   for j in range(item.loop.length())]
                  for i in range(item.loop.width())]
            item.loop.set_all_values(vv)
        return True
# fix_OXT()

def fix_H3(alias): # they are totally equivalent, so just sort them alphabetically
    idxes = [i for i, a in enumerate(alias) if a[1] in ("H","H2","H3")]
    if not idxes: return
    a0 = sorted(alias[i][0] for i in idxes)
    a1 = sorted(alias[i][1] for i in idxes)
    for i,x,y in zip(idxes, a0, a1):
        alias[i] = [x, y]

def is_dna(cc, alias):
    _, chem = check_chemtype("C2'", cc, alias)
    return chem == "CH2"

def change_group(newgr, doc):
    b = doc.find_block("comp_list")
    found = b.find_values("_chem_comp.group")
    if found[0] != newgr:
        found[0] = newgr
        return True
    return False

#groups = prepare_references(monlib_path)

def fix_group_and_add_aliases(doc, references):
    # Limitations
    # 1. if monomer can belong to multiple groups (with no changes of atom names),
    #    the last matching group will be assigned.
    # 2. only check peptide/P-peptide/M-peptide/DNA/RNA. no support for sugars
    # 3. equivalent hydrogen atoms are not sorted (H/H2/H3 of peptides)
    
    cc = read_chemcomp(doc)
    g = graph_from_chemcomp(cc)

    peptide_found = False
    doc_changed = False
    newgr = "NON-POLYMER"
    aliases = []
    
    for gr, gref in references:
        if "peptide" in gr and peptide_found: continue
        genuine, alias = check_group(g, gref, gr, cc)
        if alias:
            if "peptide" in gr:
                peptide_found = True
                # O/OXT are essentially equivalent
                ooxt = {a[1]:i for i, a in enumerate(alias) if a[1] in ("O", "OXT")}
                if len(ooxt) == 2:
                    i, j = ooxt["O"], ooxt["OXT"]
                    org_O = alias[i][0]
                    org_OXT = alias[j][0]
                    if org_OXT == "O" or org_O == "OXT":
                        print("swap alias", alias[i], alias[j])
                        alias[i][1], alias[j][1] = alias[j][1], alias[i][1]
                fix_H3(alias)
                alias = [x for x in alias if x[0]!=x[1]]
                if not alias:
                    print("INFO: {} can be {} by just swapping atoms".format(cc.name, gr))
                    genuine = True      
                if fix_OXT(cc, alias, doc):
                    doc_changed = True

            if "DNA/RNA" in gr:
                # OP1/OP2/OP3 are essentially equivalent
                ops = {a[1]:i for i, a in enumerate(alias) if a[1].startswith("OP")}
                print(alias)
                for i in ops.values():
                    n2 = alias[i][1][2]
                    r = re_OPn_or_OnP.search(alias[i][0])
                    if r:
                        n1 = [x for x in r.groups() if x][0]
                        print(n1,n2)
                        if n1 != n2:
                            j = ops["OP"+n1]
                            print("swap alias", alias[i], alias[j])
                            alias[i][1], alias[j][1] = alias[j][1], alias[i][1]

                alias = [x for x in alias if x[0]!=x[1]]
                if not alias:
                    print("INFO: {} can be {} by just swapping atoms".format(cc.name, gr))
                    genuine = True
                if fix_OP3(cc, alias, doc):
                    doc_changed = True
            if not genuine:
                aliases.append((gr, alias))
        if genuine and cc.name not in keep_nonpolymer:
            # change group
            if gr == "DNA/RNA":
                if cc.group in (gemmi.ChemComp.Group.Dna, gemmi.ChemComp.Group.Rna, gemmi.ChemComp.Group.DnaRna):
                    newgr = group_str(cc.group)
                else:
                    newgr = "DNA" if is_dna(cc, alias) else "RNA"
            else:
                newgr = gr
                
    if cc.group in non_supported_groups:
        newgr = group_str(cc.group)
    elif change_group(newgr, doc):
        doc_changed = True
            
    if aliases:
        print(cc.name, cc.group, aliases)
        loop = doc[-1].init_loop("_chem_comp_alias.", ["comp_id", "group", "atom_id", "atom_id_standard"])
        for gr, alias in aliases:
            for id1, id2 in alias:
                loop.add_row([cc.name, gr, id1, id2])
        doc_changed = True
    else:
        for item in doc[-1]:
            if "_chem_comp_alias.comp_id" in item.loop.tags:
                print(cc.name, "alias no longer needed")
                item.erase()
                break
        
    #if doc_changed:
    #    doc.write_file(f, style=gemmi.cif.Style.Aligned)

    return doc_changed, cc.name, newgr, group_str(cc.group), aliases

# fix_group_and_add_aliases()
