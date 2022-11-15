
import argparse
import datetime
import gemmi
import os
import sys
import math
import numpy
import json
from scipy.stats import pearsonr

def print_to_file(text,fout,mode='a'):
  stdout = sys.stdout
  with open(fout,mode) as sys.stdout:
    print(text)
    sys.stdout.flush()
  sys.stdout = stdout
  return

def get_rscc(map_o,map_c,pos,mask_radius):
  print("---- get_rsccs -start ---")
  mask = gemmi.FloatGrid(numpy.array(map_o, copy=False),map_o.unit_cell,map_o.spacegroup)
  mask.fill(0)
  mask.set_points_around(pos, radius=mask_radius, value=1)
  #mask.symmetrize_max()
  maskarray = numpy.array(mask, copy=False)
  points = numpy.argwhere(maskarray == 1)
  print("points.size", points.size)
  if points.size < 4:
    return None
  list_obs = []
  list_calc = []
  list_diff = []
  for point in points:
    obs =  map_o.get_point(point[0],point[1],point[2]).value
    calc = map_c.get_point(point[0],point[1],point[2]).value
    list_obs.append(obs)
    list_calc.append(calc)
    list_diff.append(obs-calc)
  n_points = len(list_obs)
  sum_obs =  round(sum(list_obs), 2)
  sum_calc = round(sum(list_calc),2)
  sum_diff = round(sum(list_diff),2)
  sum_absdiff = round(sum(abs(el) for el in list_diff),2)
  rmsd_obs_calc = round(math.sqrt(sum(el**2 for el in list_diff)/n_points),2)
  cor_obs_calc,cor_obs_calc_p = pearsonr(list_obs,list_calc)
  return cor_obs_calc

def modelComparison(model1, model2, mtzin, ligand_comp_id, opt):

  if model1 == None:
    print("modelComparison: model1 is None")
    return

  if model2 == None:
    print("modelComparison: model2 is None")
    return

  print("debug:: in modelComparison() model1", model1, "model2", model2)

  try:
    if not 'sample_rate' in opt:
      opt['sample_rate'] = 1.0
    if not 'mask_radius' in opt:
      opt['mask_radius'] = 1.5

    print("Reading model file:",model1)
    st1 = gemmi.read_structure(model1)
    print("Reading model file:",model2)
    st2 = gemmi.read_structure(model2)
    print("Reading MTZ file:",mtzin)
    mtz = gemmi.read_mtz_file(mtzin)

    print("Converting FWT/PHWT to map")
    map_obs = mtz.transform_f_phi_to_map('FWT','PHWT', sample_rate=opt['sample_rate'])
    print("Converting FC_ALL/PHIC_ALL to map")
    map_calc = mtz.transform_f_phi_to_map('FC_ALL','PHIC_ALL', sample_rate=opt['sample_rate'])

    print("Comparing equivalent atoms in all polymeric entities")
    result = []
    for model in st1:
      for chain in model:
        for residue in chain:
          if residue.name != ligand_comp_id:
            if ligand_comp_id:
              continue
          print(" Found", residue, residue.entity_type)
          if residue.entity_type == gemmi.EntityType.NonPolymer:
            for atom in residue:
              if atom.element.name == "H":
                  continue
              addr = gemmi.make_address(chain, residue, atom)
              cra1 = model.find_cra(addr)
              cra2 = st2[0].find_cra(addr)  # This always uses the first model in st2 - probably worth checking that both st1 and st2 only have one model, otherwise fail.
              print("---\ncra1: ", cra1, "cra2: ", cra2)

              if cra1 and cra2:
                print("cra1 pos:", cra1.atom.pos)
                print("cra2 pos:", cra2.atom.pos)
                peak1 = map_obs.interpolate_value(cra1.atom.pos)
                peak2 = map_obs.interpolate_value(cra2.atom.pos)
                rscc1 = get_rscc(map_obs,map_calc,cra1.atom.pos,opt['mask_radius'])
                rscc2 = get_rscc(map_obs,map_calc,cra2.atom.pos,opt['mask_radius'])

                if rscc1 == None: continue
                if rscc2 == None: continue

                result.append({'chain':chain.name,
                               'residue':str(residue.seqid.num) + (residue.seqid.icode if residue.seqid.icode != ' ' else ''),
                               'tlc':residue.name,
                               'atom':atom.name,
                               'dist':cra1.atom.pos.dist(cra2.atom.pos),
                               'b1':cra1.atom.b_iso,
                               'b2':cra2.atom.b_iso,
                               'peak1':peak1,
                               'peak2':peak2,
                               'rscc1':rscc1,
                               'rscc2':rscc2})

    print("Processed",len(result),"atoms")
    return result

  # except Exception as e:
  except FloatingPointError as e: #  I want to see the backtrace
    print("keyError: %s" % e)
    print("Cannot continue - program terminated.")
    sys.exit(1)

def printModelComparison(analysis,output_file):
  try:
    output = '\t'.join(["chain","residue","TLC","atom","dist","b1","b2","peak1","peak2","rscc1","rscc2"])
    for atom in analysis:
      output = output + '\n' + '\t'.join([atom['chain'],
                          atom['residue'],
                          atom['tlc'],
                          atom['atom'],
                          str(round(atom['dist'],2)),
                          str(round(atom['b1'],2)),
                          str(round(atom['b2'],2)),
                          str(round(atom['peak1'],2)),
                          str(round(atom['peak2'],2)),
                          str(round(atom['rscc1'],2)),
                          str(round(atom['rscc2'],2))])

    print("Writing output to file:",output_file)
    print_to_file(output,output_file,'w')

  except Exception as e:
    print("Error: %s" % e)
    print("Cannot continue - program terminated.")
    sys.exit(1)

def proc_using_code(accession_code, ligand_analysis_base):

    two_code    = accession_code[1:3]
    upcase_code = accession_code.upper()
    refmac_dir  = os.path.join(ligand_analysis_base, 'refmac')
    sub_dir     = os.path.join(refmac_dir, two_code)
    info_file   = os.path.join(sub_dir, upcase_code + ".info.json")
    print("debug:: in proc_using_code() accession_code", accession_code, "info_file", info_file)
    if os.path.exists(info_file):
        print("INFO file", info_file)
        with open(info_file, "r") as info_file_data:
            info = json.load(info_file_data)
            refmac_out_fn = info['refmac_pdb_file_name']
            ligand_comp_id = info['ligand-code']
            if os.path.exists(refmac_out_fn):
                pdb_file_name = info['pdb_file_name']
                mtz  = info['refmac_mtz_file_name']
                dict = info['acedrg_cif_file_name']
                opt = {}
                result = modelComparison(pdb_file_name, refmac_out_fn, mtz, ligand_comp_id, opt)
                results_output_file_name = "results.out"
                printModelComparison(result, results_output_file_name)
            else:
                print('Missing refmac pdb output file "{}"'.format(refmac_out_fn))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--f1', dest='model1', help="First model (PDB or mmCIF format)")
    parser.add_argument('--f2', dest='model2', help="Second model (PDB or mmCIF format)")
    parser.add_argument('--mtz', dest='mtzin', help="MTZ file containing phases (columns required: FWT, PHWT,  FC_ALL, PHIC_ALL). At the minute just one MTZ is used for model-map analysis... in future this should be extended to a second.")
    parser.add_argument('-o','--output', dest='output', required=False, help="Output file containing results of comparative analysis", default='analysis.txt')
    parser.add_argument('-s','--samplerate', dest='samplerate', required=False, help="Map sampling rate", default=1.0, type=float)
    parser.add_argument('-r','--maskradius', dest='maskradius', required=False, help="Mask radius around atomic position", default=1.5, type=float)
    parser.add_argument("--code", dest='accession_code')
    parser.add_argument('--ligand-analysis-base', dest='ligand_analysis_base')
    args = parser.parse_args()

    opt = {}
    if args.samplerate:
        opt['sample_rate'] = args.samplerate
    if args.maskradius:
        opt['mask_radius'] = args.maskradius
    if args.accession_code:
        print("INFO:: Compare models using this accession code", args.accession_code, "with base", args.ligand_analysis_base)
        proc_using_code(args.accession_code, args.ligand_analysis_base)
    else:
        time_start = datetime.datetime.now()
        ligand_code = "" # special value, means "all residue types"
        result = modelComparison(args.model1,args.model2,args.mtzin,ligand_code,opt)
        printModelComparison(result,args.output)
        time_end= datetime.datetime.now()
        print("Job completed ("+str(round((time_end-time_start).total_seconds(),3))+" seconds)")
