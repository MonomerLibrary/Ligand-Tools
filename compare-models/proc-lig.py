
import os
import json
import urllib.request
import time
import subprocess

def get_latest_ligands_json(json_file_name: str) -> bool:
    # Thank you Sameer
    s = "select?q=(q_document_type%3Alatest_chemistry)&fl=new_revised_ligand&rows=10000&omitHeader=true&wt=json"
    s = "pdb/select?group=true&group.field=pdb_id&group.ngroups=true&facet.field=new_revised_ligand&facet=true&rows=0&facet.limit=100000&facet.mincount=1&json.nl=map&group.facet=true&q=(q_document_type%3Alatest_chemistry)&omitHeader=true&wt=json"
    url = os.path.join("https://www.ebi.ac.uk/pdbe/search/", s)
    urllib.request.urlretrieve(url, json_file_name)
    # add a test that that overwrote the old one (if there was an old one)
    success = True
    return success

def mkdir_maybe(dir: str, sub_dir: str) -> None:
    dir_name = os.path.join(dir, sub_dir)
    if os.path.exists(dir_name):
        pass
    else:
        if os.path.exists(dir):
            pass
        else:
            os.mkdir(dir)
        os.mkdir(dir_name)

def tlc_to_ccd_cif_file_name(cif_dir: str, tlc: str) -> str:
    letter = tlc[0].upper()
    tlc_local_file_name = tlc + ".cif"
    fn = os.path.join(cif_dir, letter, tlc_local_file_name)
    return fn

def accession_code_to_pdb_file_name(code: str, pdb_dir: str) -> str:
    sub_dir_local = code[1:3].lower()
    pdb_fn_local = code.lower() + ".pdb"
    sub_dir = os.path.join(pdb_dir, sub_dir_local)
    fn = os.path.join(pdb_dir, sub_dir_local, pdb_fn_local)
    if os.path.exists(pdb_dir):
        pass
    else:
        os.mkdir(pdb_dir)
    if os.path.exists(sub_dir):
        pass
    else:
        os.mkdir(sub_dir)
    return fn

def fetch_ccd_cif_for(tlc: str, cif_dir: str) -> str:
    fn = tlc_to_ccd_cif_file_name(cif_dir, tlc)
    if os.path.exists(fn):
        return fn
    else:
        url = os.path.join("https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/", tlc + ".cif")
        letter = tlc[0].upper()
        mkdir_maybe(cif_dir, letter)
        urllib.request.urlretrieve(url, fn)

def get_model_db_code(tlc_cif: str) -> str:
    if os.path.exists(tlc_cif):
        with open(tlc_cif, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if "_chem_comp.pdbx_model_coordinates_db_code" in line:
                    parts = line.split()
                    if len(parts) == 2:
                        model_db_code = parts[1]
                        if model_db_code == "?":
                            return False
                        return model_db_code
            return False
    else:
        print("No file", tlc_cif)
        return False

def fetch_pdb_for_model_db_code(accession_code: str, pdb_dir: str) -> str:
    # print("fetch_pdb_for_model_db_code", accession_code)
    fn = accession_code_to_pdb_file_name(accession_code, pdb_dir)
    if os.path.exists(fn):
        return fn
    else:
        s = "https://www.ebi.ac.uk/pdbe/entry-files/download"
        fn_server = "pdb" + accession_code.lower() + ".ent"
        url = os.path.join(s, fn_server)
        print("Fetching pdb", url, fn)
        try:
            # the PDB file might not (yet) exist on the server
            fn_tmp = fn + ".partial_download"
            urllib.request.urlretrieve(url, fn_tmp)
            os.rename(fn_tmp, fn)
            return fn
        except urllib.error.HTTPError as e:
            print(e, url)
            with open("missing-models.table", "w") as f:
                f.write(url, fn)
                f.write("\n")
            return False
    return False

def ligand_has_only_organic_set_atoms(ccd_cif_file_name):

    # the Right Way to do this is with GEMMI, but this simple
    # hack will do now - so that I don't have to mess about
    # with python and installing the gemmi module

    def is_acceptable_element(ele):
        goodies = ["H", "C", "N", "I", "CL", "BR", "F", "O", "P", "S", "B"]
        return ele in goodies

    with open(ccd_cif_file_name, 'r') as f:
        lines = f.readlines()
        atoms_block = False
        for line in lines:
            if "_chem_comp_atom.pdbx_ordinal" in line:
                atoms_block = True
            if "loop_" in line:
                atoms_block = False
            if atoms_block:
                parts = line.split()
                if len(parts) > 10:
                    element = parts[3]
                    if not is_acceptable_element(element):
                        print("non-organic-set element: %s in %s" % (element, ccd_cif_file_name))
                        return False
    return True

def run_acedrg(ccd_cif_file_name: str, ligand_tlc: str, acedrg_cif_root: str) -> bool:
    args = ["acedrg", "-c", ccd_cif_file_name, "--coords", "-z", "-o", acedrg_cif_root]
    print("   ", args)
    comp_proc = subprocess.run(args, capture_output=True)
    log_dir = "acedrg-logs"
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    log_file_name = os.path.join(log_dir, "acedrg-" + ligand_tlc + ".log")
    with open(log_file_name, "bw") as f:
        f.write(comp_proc.stdout)
    return False

def convert_to_mtz(accession_code: str, sub_dir: str, data_cif_file_name: str,
                   data_mtz_file_name: str) -> bool:
    args = ["cif2mtz", "HKLIN", data_cif_file_name, "HKLOUT", data_mtz_file_name]
    s=b"END\n"
    comp_proc = subprocess.run(args, capture_output=True, input=s)
    log_file_name = os.path.join(sub_dir, "cif2mtz-" + accession_code + ".log")
    with open(log_file_name, "bw") as f:
        # print("writing cif2mtz log file: ", log_file_name)
        f.write(comp_proc.stdout)
    return comp_proc.returncode

def fetch_and_make_mtz_data(accession_code: str,  pdb_data_dir: str) -> str:
    sub_dir_local = accession_code[1:3].lower()
    data_cif_fn_local = accession_code.lower() + ".cif"
    data_mtz_fn_local = accession_code.lower() + ".mtz"
    data_cif_file_name = os.path.join(pdb_data_dir, sub_dir_local, data_cif_fn_local)
    data_mtz_file_name = os.path.join(pdb_data_dir, sub_dir_local, data_mtz_fn_local)
    data_cif_file_name_tmp = data_cif_file_name + ".partial_download"

    # note to self" check here for EXPDATA in the PDB file here. Don't try to download sfs
    # for "ELECTRON MICROSCOPY" structures

    if not os.path.exists(pdb_data_dir):
        os.mkdir(pdb_data_dir)
    if not os.path.exists(os.path.join(pdb_data_dir, sub_dir_local)):
        os.mkdir(os.path.join(pdb_data_dir, sub_dir_local))
    if os.path.exists(data_mtz_file_name):
        return data_mtz_file_name
    else:
        sub_dir = os.path.join(pdb_data_dir, sub_dir_local)
        server_dir = "https://www.ebi.ac.uk/pdbe/entry-files/download/"
        sf_file_name_local = "r" + accession_code.lower() + "sf.ent"
        url = os.path.join(server_dir, sf_file_name_local)
        try:
            result,headers = urllib.request.urlretrieve(url, data_cif_file_name_tmp)
            if os.path.exists(result):
                os.rename(data_cif_file_name_tmp, data_cif_file_name)
                print("Download", data_cif_file_name, "created")
                r = convert_to_mtz(accession_code, sub_dir, data_cif_file_name, data_mtz_file_name)
                if r != 0:
                    print("debug:: in fetch_and_make_mtz_data() with r:", r)
                if r == 0: # shell program completed successfully
                    return data_mtz_file_name
        except urllib.error.HTTPError as e:
            print("Error: in fetch_and_make_mtz_data() ", e, "for", url, data_mtz_file_name)
        return False

def write_refmac_shell_script(pdb_file_name: str, data_mtz_file_name: str, acedrg_cif_file_name: str,
               refmac_sub_dir: str, accession_code: str, ligand_tlc: str) -> bool:

    print("making refmac script for:", pdb_file_name, data_mtz_file_name, acedrg_cif_file_name)
    if data_mtz_file_name == False:
        print("No data_mtz_file - nothing doing")
        return
    if len(pdb_file_name) == 0:
        return False
    if len(data_mtz_file_name) == 0:
        return False
    if len(acedrg_cif_file_name) == 0:
        return False
    pdb_out = os.path.join(refmac_sub_dir, accession_code + "-refmac.pdb")
    hkl_out = os.path.join(refmac_sub_dir, accession_code + "-refmac.mtz")
    args = ["refmac5", "XYZIN", pdb_file_name, "HKLIN", data_mtz_file_name,
            "XYZOUT", pdb_out, "HKLOUT", hkl_out, "LIBIN", acedrg_cif_file_name]
    script_file_name_local = "refmac-" + accession_code.lower() + ".sh"
    shell_script_file_name = os.path.join(refmac_sub_dir, script_file_name_local)
    with open(shell_script_file_name, 'w') as f:
        f.write("\n")
        for arg in args:
            f.write(arg)
            f.write(" ")
        f.write("<< !\n")
        f.write("NCYCLES 5\n")
        f.write("END\n")
        f.write("!\n")

    # now the info file: # it's used for validation of the ligand
    info_file = os.path.join(refmac_sub_dir, accession_code + ".infox")
    with open(info_file, 'w') as f:
        f.write("accession_code ")
        f.write(accession_code)
        f.write("\n")
        f.write("ligand-code ")
        f.write(ligand_tlc)
        f.write("\n")
        f.write("pdb_file_name ")
        f.write(pdb_file_name)
        f.write("\n")
        f.write("data_mtz_file_name ")
        f.write(data_mtz_file_name)
        f.write("\n")
        f.write("ligand-cif ")
        f.write(acedrg_cif_file_name)
        f.write("\n")

    info_file = os.path.join(refmac_sub_dir, accession_code + ".info.json")
    dict = {}
    dict['accession_code'] = accession_code
    dict['ligand-code'] = ligand_tlc
    dict['pdb_file_name'] = pdb_file_name
    dict['refmac_pdb_file_name'] = pdb_out
    dict['data_mtz_file_name'] = data_mtz_file_name
    dict['refmac_mtz_file_name'] = hkl_out
    dict['acedrg_cif_file_name'] = acedrg_cif_file_name
    j = json.dumps(dict, indent = 4)
    with open(info_file, 'w') as f:
        f.write(j)
        f.write("\n")

def is_em(pdb_file_name: str) -> bool:
    if os.path.exists(pdb_file_name):
        f = open(pdb_file_name)
        lines = f.readlines()
        for line in lines:
            if "EXPDTA" in line:
                return "ELECTRON MICROSCOPY" in line
    return False

def proc_model_and_ligand(accession_code: str, pdb_file_name: str, ligand_tlc: str,
                          cif_dir: str, acedrg_dir: str, refmac_dir: str) -> bool:

    # Don't continue if the pdb_file_name is False or empty
    try:
        if len(pdb_file_name) == 0:
            return False
    except Exception as e:
        print("in proc_model_and_ligand()", e)
        return False

    if os.path.exists(pdb_file_name):
        letter = ligand_tlc[0].lower()
        acedrg_cif_local  = ligand_tlc + "-acedrg.cif"
        acedrg_root_local = ligand_tlc + "-acedrg"
        acedrg_cif_root      = os.path.join(acedrg_dir, letter, acedrg_root_local)
        acedrg_cif_file_name = os.path.join(acedrg_dir, letter, acedrg_cif_local)
        if os.path.exists(acedrg_cif_file_name):
            pass
        else:
            ccd_cif_file_name = tlc_to_ccd_cif_file_name(cif_dir, ligand_tlc)
            if ligand_has_only_organic_set_atoms(ccd_cif_file_name):
                run_acedrg(ccd_cif_file_name, ligand_tlc, acedrg_cif_root)
        is_em_flag = is_em(pdb_file_name)
        if is_em_flag:
            print("EM data for", accession_code)
        else:
            data_mtz_file_name = fetch_and_make_mtz_data(accession_code, pdb_data_dir)
            # print("data_mtz_file_name", data_mtz_file_name)
            sub_dir = os.path.join(refmac_dir, accession_code[1:3].lower())
            if not os.path.exists(refmac_dir):
                os.mkdir(refmac_dir)
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)

            try:
                # print("data_mtz_file_name: ", data_mtz_file_name, len(data_mtz_file_name))
                if os.path.exists(acedrg_cif_file_name):
                    if os.path.exists(pdb_file_name):
                        if os.path.exists(data_mtz_file_name):
                            write_refmac_shell_script(pdb_file_name, data_mtz_file_name,
                                                      acedrg_cif_file_name,
                                                      sub_dir, accession_code, ligand_tlc)
                        else:
                            print("missing mtz", data_mtz_file_name)
            except TypeError as e:
                print(e, "for in write_refmac_shell_script")
        return False
    else:
        print("pdb file does not exist: " + pdb_file_name)
        return False

def parse_latest_ligands_json(j: json, cif_dir: str, pdb_dir: str, acedrg_dir: str, refmac_dir: str):

    grouped = j["grouped"]
    fc = j['facet_counts']
    ff = fc["facet_fields"]
    nrl = ff["new_revised_ligand"]
    for lig_descr in nrl:
        s = lig_descr
        tlc = s[:3]
        print("Processing: -----------", tlc, "--------------")
        tlc_cif = tlc_to_ccd_cif_file_name(cif_dir, tlc)
        fetch_ccd_cif_for(tlc, cif_dir)
        try:
            model_db_code = get_model_db_code(tlc_cif)
            pdb_fn = fetch_pdb_for_model_db_code(model_db_code, pdb_dir)
            if len(pdb_fn) > 0:
                if os.path.exists(pdb_fn):
                    proc_model_and_ligand(model_db_code, pdb_fn, tlc, cif_dir, acedrg_dir, refmac_dir)
        except TypeError as e:
            print(e, "for", tlc)

def new_ligands_from_pdbe_json(json_file_name, cif_dir, pdb_dir: str, acedrg_dir: str, refmac_dir: str):
    try:
        with open(json_file_name, 'r') as f:
            j = json.load(f)
            parse_latest_ligands_json(j, cif_dir, pdb_dir, acedrg_dir, refmac_dir)
    except TypeError as e:
        print(e)

def need_new_json_for_ligands(json_file_name: str, cif_dir: str) -> bool:
    if os.path.exists(json_file_name):
        # if it's more than a day old, get a new one
        st=os.stat(json_file_name)
        mtime=st.st_mtime
        t = time.time()
        delta_time = t - mtime
        return delta_time > 24 * 60 * 60
    else:
        status = get_latest_ligands_json(json_file_name)
        return False

if __name__ == '__main__':

    cif_dir = "CCD-cifs"
    pdb_dir = "pdb"
    pdb_data_dir = "pdb-data"
    acedrg_dir = "acedrg"
    refmac_dir = "refmac"
    json_file_name = "new_revised_ligands.json"

    if need_new_json_for_ligands(json_file_name, cif_dir):
        get_latest_ligands_json(json_file_name) # download the json file

    new_ligands_from_pdbe_json(json_file_name, cif_dir, pdb_dir, acedrg_dir, refmac_dir)
