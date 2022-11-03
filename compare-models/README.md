# compare_models.py

Tool to compare two models, which are assumed to have the same atomic composition and nomenclature. Various per-atom statistics are extracted - B-factor, density peak height, RSCC, etc.

At present, only a single MTZ file is used for calculation of density-based statistics corresponding to the positions of atoms in the two models.

## For full list of options:

python3 compare_models.py --help

## Example of usage:

python3 compare_models.py --f1 3b50.cif --f2 3b50_besttls.pdb --mtz 3b50_besttls.mtz

Where:
* 3b50.cif is taken from the PDB
* 3b50_besttls.pdb is taken from PDB_REDO
* 3b50_besttls.mtz is taken from PDB_REDO
