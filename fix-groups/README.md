# Post-acedrg tools

Group name of monomer (`_chem_comp.group`) specifies what links may be applied in the refinement.
For example, monomers of "peptide" group can have TRANS or CIS link in the polymer chain.
However, it only works if atom names follow a standard nomenclature.
Otherwise `_chem_comp.group` needs to be "NON-POLYMER" and an alias may be added, which is a translation table of atom names to the standard.

You can fix group names and add aliases (if applicable):
```sh
python post_acedrg.py AcedrgOut.cif
```

At the moment, peptide/P-peptide/M-peptide/DNA/RNA are only supported. In other words, there is no support for carbohydrates.

## Example: CIR
CIR is a peptide monomer, but the atom names does not follow the standard.
```sh
wget https://files.rcsb.org/ligands/download/CIR.cif
acedrg -c CIR.cif 
python post_acedrg.py AcedrgOut.cif 
```
This results in NON-POLYMER group with the following alias:
```
loop_
_chem_comp_alias.comp_id
_chem_comp_alias.group
_chem_comp_alias.atom_id
_chem_comp_alias.atom_id_standard
CIR peptide N2   N
CIR peptide C2   CA
CIR peptide C1   C
CIR peptide O1   O
CIR peptide O2   OXT
CIR peptide HN22 H
CIR peptide H23  H2
CIR peptide HN21 H3
```

# monomer library maintenance

You can update all monomers:
```sh
python post_acedrg.py --do_not_keep_org /path-to-monlib/?/
```

Then update mon_lib_list.cif 
```
python fix_monlib_list.py --monlib /path-to-monlib/
```
This script check a group name from each cif file and update mon_lib_list.cif to make group names consistent. 
