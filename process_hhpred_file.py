#! /usr/bin/python

import sys
import subprocess
import os
import cctbx
import iotbx.bioinformatics
import re
import iotbx.pdb
from iotbx.pdb import pdbe

# Script to process an output from HHPred server
# Requires cctbx and phenix

if len(sys.argv)==1:
    sys.exit("Usage:\n process_hhpred_file.py file.hhr")
elif len(sys.argv)>=2: # I migth be interested in adding more parameters
    name_file=sys.argv[1]

min_ident=0.1  # minimum of sequence identity required to be taken into account

wd=os.getcwd()

hhr=open(name_file).read()
hsearch=iotbx.bioinformatics.hhsearch_parser(output=hhr)
count_hit=0
# Prepare a table to fill with the information from all the hits
table_file=open("table_hits.txt",'w')
table_file.write('%-4s %-10s %-15s %-15s %-10s %-25s %-10s %-10s %-10s %-10s %-10s %-10s\n' % ('PDB','Identity','E-value','Align lenght','PDB chain','Experimental method','A','B','C','Alpha','Beta','Gamma'))

del table_file


list_models_processed=[]

for hit in hsearch.hits():
  try:
	  print("*************************************************************************************************************************")
	  print("\n\n\nProcessing...\n")
	  print(hit.annotation)
	  count_hit=count_hit+1

	  # Get alignement information
	  id_pdb=hit.identifier
	  chain_pdb=hit.chain
	  aln=hit.alignment
          # First of all check that evalue is significative at all
          p = subprocess.Popen('grep '+id_pdb + ' ' + name_file, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=wd)
          complete_output, errors = p.communicate()
          line_values=complete_output.split("\n")[0]
          evalue_string=line_values[41:48]
          evalue=float(evalue_string)
          if evalue>1:
              print '\n e-value for this alignement is larger than 1, this PDB will not be further processed'
              continue 
	  #print '\nAlignement',aln
	  #print '\nTarget aligned sequence',aln.alignments[0]
	  target_aln=aln.alignments[0]
	  #print '\nHit aligned sequence',aln.alignments[1]
	  hit_aln=aln.alignments[1]
	  #print '\nNumber of identities',aln.identity_count()
	  nid=aln.identity_count()
	  #print '\nFraction identities',aln.identity_fraction()
	  identity=round(aln.identity_fraction(),3)
	  if identity<min_ident:
	  	print '\n Identity is less than ',min_ident,' this PDB will not be further processed'
	  	continue
	  #print '\nLength of the alignement',aln.length()
	  length_aln=aln.length()
	  #print '\nAlignment multiplicity',aln.multiplicity() 
	  aln_multiplicity=aln.multiplicity()
          # Check if this model and chain are already downloaded
          if id_pdb+chain_pdb in list_models_processed:
              print 'This model is already processed ',id_pdb+chain_pdb
              continue
          else:
              list_models_processed.append(id_pdb+chain_pdb)
	  # Fetch the pdb of the hit and do some modifications according to the alignement
	  hit_folder=os.path.join(os.getcwd(),'hhpred_hit_'+str(count_hit))
	  os.makedirs(hit_folder)   
	  p = subprocess.Popen(["phenix.fetch_pdb",id_pdb], stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=hit_folder)
	  out, err = p.communicate()
	  print '\n',out
          # TODO: Remove non aligning parts and chains
	  # Get crystallographic information if experimental method is crystallographic
	  obj=iotbx.pdb.pdbe.get_request( 'http://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/',id_pdb )
	  experimental_method=obj[0]['experimental_method']
	  if experimental_method=='X-ray diffraction':
              # Keys are: [u'resolution_low', u'starting_model', u'r_free_selection_details', u'resolution_high', u'r_free_percent_reflections', u'r_factor', u'r_free', u'refinement_software', u'completeness', u'cell', u'percent_reflections_observed', u'diffraction_experiment', u'expression_host_scientific_name', u'crystal_growth', u'experimental_method', u'num_reflections', u'phasing_method', u'experiment_data_available', u'experimental_method_class', u'r_work', u'spacegroup', u'resolution', u'structure_determination_method']
              
	      resolution=obj[0]['resolution_high']
	      cell=obj[0]['cell']
	  else:
	      resolution=None
	      cell=None
	  # Write the information in the table
	  table_file=open("table_hits.txt",'a')
          table_file.write('%-4s %-10s %-15s %-15s %-10s %-25s %-10s %-10s %-10s %-10s %-10s %-10s\n' % (id_pdb,identity,evalue,length_aln,chain_pdb,experimental_method,cell['a'],cell['b'],cell['c'],cell['alpha'],cell['beta'],cell['gamma']))
	  del table_file
  except:
  	print '\nThere was some error while trying to fetch this pdb and its information'
        print 'Continuing with the next pdb from the list\n'
  	continue



