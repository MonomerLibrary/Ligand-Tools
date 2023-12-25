#!/usr/bin/env python3

# This scripts changes atom names and adds _chem_comp_atom.alt_atom_id.
# It was written for renaming atoms after the remediation of atom names
# in peptide backbones, which was done in the PDB in Nov 2023.
# Here is what it does:
# - reads monomer codes from stdin
# - reads CCD file to obtain alt_atom_id -> atom_id mapping
# - finds monomers in MONOMER_LIBRARY directory
# - adds _chem_comp_atom.alt_atom_id as a copy of _chem_comp_atom.atom_id
# - uses the CCD mapping to rename all values with '.atom_id' in tag,
#   except _chem_comp_alias.atom_id_standard
# - fix links and modifications in mon_lib_list.cif

# Call it with two arguments and provide codes from stdin:
# python rename_atoms.py $CLIBD_MON components.cif.gz < monomer_codes.txt
# Warning: such call overwrites monomer files in $CLIBD_MON.

import sys
import os
import re
from typing import Dict
from gemmi import cif, MonLib

# two arguments are expected
MONOMER_LIBRARY, CCD = sys.argv[1:]

ID_TAG = '_chem_comp_atom.atom_id'
ALT_TAG = '_chem_comp_atom.alt_atom_id'

def main() -> None:
    print('reading CCD from', CCD, '...')
    ccd = cif.read(CCD)
    print('updating', MONOMER_LIBRARY, '...')
    monlib = MonLib()
    monlib.monomer_dir = MONOMER_LIBRARY
    path_list = os.path.join(MONOMER_LIBRARY, "list", "mon_lib_list.cif")
    doc_list = cif.read(path_list)
    links, mods = check_list(doc_list)
    for line in sys.stdin:
        code = line.strip()
        path = monlib.path(code)
        print(code, end='  ')
        renaming = read_old_new_mapping(ccd[code])
        if not renaming:
            print('--- nothing to rename')
            continue
        print(' '.join('%s->%s' % item for item in renaming.items()))
        if remediate(path, renaming, ccd[code]):
            remediate_list(doc_list, code, links, mods, renaming)
    doc_list.write_file(path_list)

def read_old_new_mapping(ccd_block : cif.Block) -> Dict[str,str]:
    renaming = {}
    for row in ccd_block.find([ALT_TAG, ID_TAG]):
        # don't unquote the new value, it will be copied
        old, new = row.str(0), row[1]
        if old and old != new:
            renaming[old] = new
    return renaming


def remediate(path : str, renaming : Dict[str, str], ccd_block : cif.Block) -> None:
    doc = cif.read(path)
    block = doc[-1]
    id_col : cif.Column = block.find_loop(ID_TAG)
    atom_loop = id_col.get_loop()
    assert atom_loop
    if ALT_TAG in atom_loop.tags:
        print(' --- already has', ALT_TAG, '- skipping')
        return
    if is_peptide(doc):
        renaming.update(check_peptide_n_h(ccd_block, block))

    pos = atom_loop.tags.index(ID_TAG)
    atom_loop.add_columns([ALT_TAG], '.', pos + 1)
    n = atom_loop.length()
    for i in range(n):
        atom_loop[i, pos+1] = atom_loop[i, pos]
    for item in block:
        table = block.item_as_table(item)
        for i in range(table.width()):
            remediate_column(table.column(i), renaming)

    # sanity checks
    new_names = [atom_loop[i, pos] for i in range(n)]
    if all(new_names[i] == atom_loop[i, pos+1] for i in range(n)):
        print(' --- atom names not changed - skipping')
        return
    if len(set(new_names)) != n:
        print(' --- skip, it would have duplicated atom names:',
              ' '.join(x for x in set(new_names) if new_names.count(x) > 1))
        return

    doc.write_file(path, style=cif.Style.Aligned)
    return True

def remediate_column(column : cif.Column, renaming : Dict[str, str]) -> None:
    tag = column.tag
    if '.atom_id' in tag and tag != '_chem_comp_alias.atom_id_standard':
        for n, atom_id in enumerate(column):
            unquoted_id = cif.as_string(atom_id)
            if unquoted_id in renaming:
                column[n] = renaming[unquoted_id]

def is_peptide(doc : cif.Document):
    if doc.find_block("comp_list").find_value("_chem_comp.group") == "peptide":
        return True
    if any(x == "peptide" for x in doc[-1].find_values("_chem_comp_alias.group")):
        return True
    return False

def check_peptide_n_h(ccd_block, monlib_block):
    names_monlib = {"N": "N", "H": "H", "H2": "H2", "H3": "H3"}
    for r in monlib_block.find("_chem_comp_alias.", ["group", "atom_id", "atom_id_standard"]):
        if cif.as_string(r[0]) == "peptide" and cif.as_string(r[2]) in names_monlib:
            names_monlib[cif.as_string(r[2])] = cif.as_string(r[1])

    # remake renaming, to include everything
    renaming = {row.str(0): row.str(1) for row in ccd_block.find([ALT_TAG, ID_TAG])}

    if names_monlib["N"] not in renaming:
        print("error {} not found in {} (ccd)".format(names_monlib["N"], ccd_block.name))
        return {}
    N_ccd = renaming[names_monlib["N"]]
    NH_ccd = []
    Hs_ccd = {cif.as_string(r[0]) for r in ccd_block.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
              if cif.as_string(r[1]) == "H"}
    for r in ccd_block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2"]):
        a1, a2 = cif.as_string(r[0]), cif.as_string(r[1])
        if a1 == N_ccd and a2 in Hs_ccd:
            NH_ccd.append(a2)
        elif a2 == N_ccd and a1 in Hs_ccd:
            NH_ccd.append(a1)
    #print("Hs_ccd", Hs_ccd, "NH_ccd", NH_ccd)
    ret = {}
    for h in "H", "H2", "H3":
        if names_monlib[h] not in renaming or renaming.get(names_monlib[h]) not in NH_ccd:
            print("unmatched H", ccd_block.name, names_monlib[h])
            if h not in Hs_ccd:
                ret[names_monlib[h]] = h
                print("debug", ccd_block.name, ret)
            break
    return ret

def check_list(doc: cif.Document):
    links, mods = {}, {}
    for r in doc.find_block("link_list").find("_chem_link.", ["id", "comp_id_1", "comp_id_2"]):
        for i in (1, 2):
           comp = cif.as_string(r[i])
           if comp:
               links.setdefault(comp, []).append((cif.as_string(r[0]), i))
    for r in doc.find_block("mod_list").find("_chem_mod.", ["id", "comp_id"]):
        comp = cif.as_string(r[1])
        if comp:
            mods.setdefault(comp, []).append(cif.as_string(r[0]))
    return links, mods

def remediate_list(doc : cif.Document, code : str, links, mods, renaming : Dict[str, str]):
    for link, idx in links.get(code, ()):
        block = doc.find_block("link_" + link)
        for item in block:
            table = block.item_as_table(item)
            cols = []
            for i, tag in enumerate(table.tags):
                r = re.search("atom(_[0-9]|)_comp_id", tag)
                if r:
                    cola = table.find_column(tag[:tag.rindex(".")] + ".atom_id" + r.group(1))
                    cols.append((table.column(i), cola))
            for colc, cola in cols:
                for i, c in enumerate(colc):
                    if cif.as_int(c) == idx:
                        unquoted_id = cif.as_string(cola[i])
                        if unquoted_id in renaming:
                            cola[i] = renaming[unquoted_id]
    for mod in mods.get(code, ()):
        block = doc.find_block("mod_" + mod)
        for item in block:
            table = block.item_as_table(item)
            for i in range(table.width()):
                remediate_column(table.column(i), renaming)

main()
