import argparse
import datetime
import gemmi
import sys
import math
import numpy
from scipy.stats import pearsonr

def print_to_file(text,fout,mode='a'):
  stdout = sys.stdout
  with open(fout,mode) as sys.stdout:
    print(text)
    sys.stdout.flush()
  sys.stdout = stdout
  return
   
def get_rscc(map_o,map_c,pos,mask_radius):
  mask = gemmi.FloatGrid(numpy.array(map_o, copy=False),map_o.unit_cell,map_o.spacegroup)
  mask.fill(0)
  mask.set_points_around(pos, radius=mask_radius, value=1)
  #mask.symmetrize_max()
  maskarray = numpy.array(mask, copy=False)
  points = numpy.argwhere(maskarray == 1)
  list_obs = []
  list_calc = []
  list_diff = []
  for point in points:
    obs = map_o.get_point(point[0],point[1],point[2]).value
    calc = map_c.get_point(point[0],point[1],point[2]).value
    list_obs.append(obs)
    list_calc.append(calc)
    list_diff.append(obs-calc)
  n_points = len(list_obs)
  sum_obs = round(sum(list_obs),2)
  sum_calc = round(sum(list_calc),2)
  sum_diff = round(sum(list_diff),2)
  sum_absdiff = round(sum(abs(el) for el in list_diff),2)
  rmsd_obs_calc = round(math.sqrt(sum(el**2 for el in list_diff)/n_points),2)
  cor_obs_calc,cor_obs_calc_p = pearsonr(list_obs,list_calc)
  return cor_obs_calc

def modelComparison(model1,model2,mtzin,opt):
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
          if residue.entity_type == gemmi.EntityType.Polymer:
            for atom in residue:
              addr = gemmi.make_address(chain, residue, atom)
              cra1 = model.find_cra(addr)
              cra2 = st2[0].find_cra(addr)  # This always uses the first model in st2 - probably worth checking that both st1 and st2 only have one model, otherwise fail.
              if cra1 and cra2:
                peak1 = map_obs.interpolate_value(cra1.atom.pos)
                peak2 = map_obs.interpolate_value(cra2.atom.pos)
                rscc1 = get_rscc(map_obs,map_calc,cra1.atom.pos,opt['mask_radius'])
                rscc2 = get_rscc(map_obs,map_calc,cra2.atom.pos,opt['mask_radius'])
              
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

  except Exception as e:
    print("Error: %s" % e)
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

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--f1', dest='model1', required=True, help="First model (PDB or mmCIF format)")
  parser.add_argument('--f2', dest='model2', required=True, help="Second model (PDB or mmCIF format)")
  parser.add_argument('--mtz', dest='mtzin', required=True, help="MTZ file containing phases (columns required: FWT, PHWT,  FC_ALL, PHIC_ALL). At the minute just one MTZ is used for model-map analysis... in future this should be extended to a second.")
  parser.add_argument('-o','--output', dest='output', required=False, help="Output file containing results of comparative analysis", default='analysis.txt')
  parser.add_argument('-s','--samplerate', dest='samplerate', required=False, help="Map sampling rate", default=1.0, type=float)
  parser.add_argument('-r','--maskradius', dest='maskradius', required=False, help="Mask radius around atomic position", default=1.5, type=float)
  args = parser.parse_args()

  opt = {}
  if args.samplerate:
    opt['sample_rate'] = args.samplerate
  if args.maskradius:
    opt['mask_radius'] = args.maskradius

  time_start = datetime.datetime.now()
  result = modelComparison(args.model1,args.model2,args.mtzin,opt)
  printModelComparison(result,args.output)
  time_end= datetime.datetime.now()

  print("Job completed ("+str(round((time_end-time_start).total_seconds(),3))+" seconds)")
