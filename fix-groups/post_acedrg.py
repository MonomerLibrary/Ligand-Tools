"""
Author: "Keitaro Yamashita, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, generators
import argparse
import os
import glob
import json
import gemmi
import groups

parser = argparse.ArgumentParser(prog="post_acedrg",
                                 description="Fix group name and add alias table if needed")
parser.add_argument("targets", nargs='+',
                    help="files or directories")
parser.add_argument("--monlib",
                    help="Monomer library path. Default: $CLIBD_MON")
parser.add_argument("--json_out", default="post_acedrg.json",
                    help="output json file name")
parser.add_argument("--do_not_keep_org", action="store_true",
                    help="By default it overwrites target files and keep .org files.")

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

def main(args):
    if args.monlib is None: args.monlib = os.environ["CLIBD_MON"]
    references = groups.prepare_references(args.monlib)    
    files = find_files(args.targets)
    ret = []
    for f in files:
        print("processing", f)
        doc = gemmi.cif.read(f)
        if not doc[-1].name.startswith("comp_"):
            print("WARNING: {} is not a monomer dictionary".format(f))
            continue

        doc_changed, cc_name, newgr, oldgr, aliases = groups.fix_group_and_add_aliases(doc, references)
        ret.append((f, cc_name, {"updated": doc_changed, "new_group": newgr,
                                 "old_group": oldgr.replace("non-polymer", "NON-POLYMER"),
                                 "aliases": aliases}))
        if doc_changed:
            if not args.do_not_keep_org: os.rename(f, f + ".org")
            doc.write_file(f, style=gemmi.cif.Style.Aligned)
            
    json.dump(ret, open(args.json_out, "w"), indent=True)
# main()

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
