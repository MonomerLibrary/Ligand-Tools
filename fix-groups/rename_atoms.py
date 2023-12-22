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

# Call it with two arguments and provide codes from stdin:
# python rename_atoms.py $CLIBD_MON components.cif.gz < monomer_codes.txt
# Warning: such call overwrites monomer files in $CLIBD_MON.

import sys
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
    for line in sys.stdin:
        code = line.strip()
        path = monlib.path(code)
        print(code, end='  ')
        renaming = read_old_new_mapping(ccd[code])
        if not renaming:
            print('--- nothing to rename')
            continue
        print(' '.join('%s->%s' % item for item in renaming.items()))
        remediate(path, renaming)

def read_old_new_mapping(ccd_block : cif.Block) -> Dict[str,str]:
    renaming = {}
    for row in ccd_block.find([ALT_TAG, ID_TAG]):
        # don't unquote the new value, it will be copied
        old, new = row.str(0), row[1]
        if old and old != new:
            renaming[old] = new
    return renaming


def remediate(path : str, renaming : Dict[str, str]) -> None:
    doc = cif.read(path)
    block = doc[-1]
    id_col : cif.Column = block.find_loop(ID_TAG)
    atom_loop = id_col.get_loop()
    assert atom_loop
    if ALT_TAG in atom_loop.tags:
        print(' --- already has', ALT_TAG, '- skipping')
        return
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

    doc.write_file(path)

def remediate_column(column : cif.Column, renaming : Dict[str, str]) -> None:
    tag = column.tag
    if '.atom_id' in tag and tag != '_chem_comp_alias.atom_id_standard':
        for n, atom_id in enumerate(column):
            unquoted_id = cif.as_string(atom_id)
            if unquoted_id in renaming:
                column[n] = renaming[unquoted_id]

main()
