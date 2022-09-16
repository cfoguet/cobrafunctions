
import getopt
import sys
import operator

import cobra
import os
import json
from cobrafunctions.gim3e import integrate_omics_gim3e_and_remove
from cobrafunctions.cobra_functions import sampling,sampling_matrix_get_mean_sd,get_equation
from cobrafunctions.read_spreadsheets import read_spreadsheets
from cobrafunctions.write_spreadsheet import write_spreadsheet


"""
gim3e_and_sampling.py
Average organ transcript abundances are mapped to the organ-specific subnetworks using the gene reaction annotations in each network. More in detail, transcripts abundances of isoenzymes and enzymes subunits catalysing each reaction or transport processes are added and, subsequently, LOG2 transformed. The resulting values are used as input to apply the GIM3E algorithm. The GIM3E algorithm applies a flux minimization weighted by transcript abundance allowing to identify solutions that are both enzymatically efficient and consistent with gene expression data. Subsequently, the flux ranges within 99% of the GIM3E optimal solution are identified and the resulting solution space is sampled using the Artificially Centered hit-and-run (ACHR) algorithm. The average of these flux samples can be used as the reference or average flux map of each organ.

Usage: gim3e_and_sampling.py [INPUTS...] 

INPUTS:
-m, --organ_specific_model_directory  : Directory containing organ-specific genome-scale metabolic models in SBML format. All SBML models will be analyzed in succession. 
-r, --reference_transcript_abundance :  Path to the CSV or XLSX file defining the average organ gene expression in TPM or FPKM.
-o, --output_directory : Output directory. Will be created if it does not exist. 

OUTPUTS:
Reference_fluxes_ORGAN_NAME.csv : A CSV file with the GIM3E optimal solution, GIM3E flux ranges and average flux distributions for all reactions in all analysed organs. 
Reference_fluxes_ORGAN_NAME_dict.json : A json file with the average flux distributions for all reactions in all analysed organs. It is used as input for run_qMTA.py
ORGAN_NAME_gim3e__constrained_model.sbml : Organ-specific model constrained to the Gim3e flux ranges. Used as input for run_qMTA.py. 

"""
opts, args = getopt.getopt(sys.argv[1:],"r:m:o:",["reference_transcript_abundance=","organ_specific_model_directory=","output_directory="])
for opt, arg in opts:
      #print opt,arg
      if opt in ("-r", "--reference_transcript_abundance"):
         gene_expression_file = arg
      elif opt in ("-m", "--organ_specific_model_directory"):
          sbml_folder= arg
          os.chdir(sbml_folder)
          model_dict={}
          SBML_files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f)) and (".sbml" in f or ".xml" in f)]
          for SBML_file in SBML_files:
              tissue_key=SBML_file.replace(".sbml","").replace(".SBML","").replace(".XML","").replace(".xml","").replace("organ_metabolic_network_","").replace("organ_specific_metabolic_network_","")
              model_dict[tissue_key]=cobra.io.read_sbml_model(SBML_file)         
      elif opt in ("-o", "--output_directory"):
          output_folder= arg

#Change to output directory
if not os.path.exists(output_folder):
             os.makedirs(output_folder)
os.chdir(output_folder)

gene_expression_file=file_name=gene_expression_file



conditions_of_interest=model_dict.keys() #Conditions that will be analyzed with GIMME
conditions_to_sample=conditions_to_fva=conditions_of_interest # #Conditions where with FVA and sampling will be performed

###########Gene expression inpputs
output_prefix="Reference_fluxes_"
gpr_mode="full" #By default leave it to full to use full GPR rules. Alternatively, you can set it to the "average" to use to average gene expression for the genes associated to each reaction 
or_mode="sum" #When there are multiple isoenzymes indicated with OR the gene expression of each isoform will be added. ONLY USE SUM IF GENE EXPRESSION DATA IS NOT IN LOG SCALE in gene_expression_file
convert_log2=True #After computing the gene expression associated to each reaction, it converts it to log2 Scale. It uses the expression math.log(reaction_expression_dict[x]+0.1,2) #WARNING SET TO FALSE IF DATA ALREADY IN LOG2 
absent_gene_expression=5 #Percentile of expression of metabolic genes assigned to genes without experimental measurments (Note that in RNA-SEQ all genes should have experimental measures) 
absent_reaction_expression=100 #Percentile assigned to reactions without genes
##Gene prefix and sufix
#Example: If in the SBML model genes are defined as G_10101.1 were 10101 is gene ID, prefix="G_" , sufix="."
gene_prefix="" 
gene_sufix=""
replace_and_with_or=True


##Gimme parameters
#Minimization_weight=max(0,low_expression_threshold-reaction_expression)+base_penalty
low_expression_threshold=95 #Gene expression pecentile that sets threshold for minimization
base_penalty=1 #Base penalty to ensure that even reactions without gene expression are also minimized
penalty_precision=2 #Precision used for weights in GIMME
correct_for_complexes=True #correct for reactions like '3.0 6pgc_c + 3.0 nadp_c --> 3.0 co2_c + 3.0 nadph_c + 3.0 ru5p__D_c' . When applying GIMME it will give this reaction x3 the minimization weight. This prevents GIMME from prioritizing this kind of reactions
fraction_of_optimum_biomass=0 #Minimumm fraction of biomass when minimizing reactions with GIMME
force_run=True #Force GIMME to run even if the output models exist
gim3e_fraction_optimum=0.99#How much the GIMME objective can deviate when running with FVA. It shoult be set close to 1 to prevent large flux through loops
#####Sampling parameters
n_samples=1000 #How many flux samples will be computed
thinning=1000   #Thining parameter in artifical hit and run algorythm. Larger values provide more reprsentative flux samples



gene_expression=read_spreadsheets(gene_expression_file)
sheet_name=gene_expression.keys()[0]
header=gene_expression[sheet_name][0]

#Find the columns of interest
column_dict={}
for n,col in enumerate(header):
    col_feeding=col+"_feeding"
    col_fasting=col+"_fasting"
    #print col_feeding,col_fasting
    if n==0:
       continue
    elif col in conditions_of_interest:
       column_dict[col]=n 
    if col_feeding in conditions_of_interest:
       column_dict[col_feeding]=n 
    if col_fasting in conditions_of_interest:
       column_dict[col_fasting]=n 

print column_dict
del(gene_expression)


out_list=[]
low_reaction_dict={}
sample_sampling_dict={}
fva_dict={}
reaction_expression_dict_dict={}
optimal_solution_dict_dict={}
penalty_dict_dict={}
th_dict={}


output_sheet=[["Reaction id","Reaction name","Reaction","Reaction","Optimal solution","Minimum","Maximum","Reaction expression","Sampling mean","Sampling SD"]]

for sample in sorted(conditions_of_interest):
    if sample in conditions_to_fva:
        run_fva=True
    else:
        run_fva=False
    n=column_dict[sample]
    model=model_dict[sample].copy()
    #Remove empty reactions
    for reaction in list(model.reactions):
        if len(reaction.metabolites)==0:
           reaction.remove_from_model()
    if replace_and_with_or:
     for reaction in model.reactions:
         gene_str=""
         for n_reaction_gene,gene in enumerate(reaction.genes):
             if n_reaction_gene==0:
                gene_str=gene.id
             else:
                gene_str+=" or "+gene.id  
         reaction.gene_reaction_rule=gene_str
    reaction_ids_to_omit=[x.id for x in model.reactions.query("RGROUP_")] 
    penalty_dict,gene_expression_model,objective,solution_dict,fva, low_reactions,core_reactions,low_th,reaction_expression_dict, gene_expression_dict, value_list=integrate_omics_gim3e_and_remove(model,file_name,fraction_of_optimum=fraction_of_optimum_biomass,low_expression_threshold=low_expression_threshold,absent_gene_expression=5,absent_reaction_expression=100,percentile=True,gene_method="average",gene_prefix=gene_prefix,gene_sufix=gene_sufix,metabolite_list_fname=None,label_model=None,epsilon=0.0001,gim3e_fraction_optimum=gim3e_fraction_optimum,run_fva=run_fva,add_as_constraints=True,boundaries_precision=0.00001,all_gene_exp_4_percentile=False,base_penalty=base_penalty,convert_solution_dict_to_reversible=True,gpr_mode=gpr_mode,or_mode=or_mode,omit_0=True,remove_reactions_below=0,gene_value_col=n,penalty_precision=penalty_precision,reactions_to_keep=[],log2=convert_log2,correct_for_complexes=correct_for_complexes,reaction_ids_to_omit=reaction_ids_to_omit)
    optimal_solution_dict_dict[sample]=solution_dict
    penalty_dict_dict[sample]=penalty_dict
    th_dict[sample]=low_th
    low_reaction_dict[sample]=low_reactions
    fva_dict[sample]=fva
    reaction_expression_dict_dict[sample]=reaction_expression_dict
    cobra.io.write_sbml_model(model,sample+"_gim3e__constrained_model.sbml")
    if sample in conditions_to_sample:
       aggregated_results, reaction_ids=sampling(model,n=n_samples,processes=1,objective=None,starts=1,return_matrix=True,method="achr",thinning=thinning)
       reaction_n_dict={}
       for n,rid in enumerate(reaction_ids):
           reaction_n_dict[rid]=n            
       reversible_matrix=[]
       reversible_reaction_ids=[]
       reversible_matrix=aggregated_results
       reversible_reaction_ids=reaction_ids
       stat_dict=sampling_matrix_get_mean_sd(reversible_matrix,reversible_reaction_ids,include_absolute_val_stats=True)
       sample_sampling_dict[sample]=stat_dict
    else:
        stat_dict={}
    for reaction in model.reactions: 
        rid=reaction.id
        try:
           fva_max=fva.get(reaction.id)["maximum"] 
           fva_min=fva.get(reaction.id)["minimum"]
        except:
           fva_max=""
           fva_min="" 
        try:
           mean=stat_dict[reaction.id]["mean"]
           std=stat_dict[reaction.id]["std"]
        except:
           mean=""
           std=""
        try:
           equation=get_equation(model,rid,include_compartment=True)
        except:
           equation=""
        row=[rid,reaction.name,reaction.reaction,equation,solution_dict.get(rid),fva_min,fva_max,reaction_expression_dict.get(reaction.id),mean,std] 
        output_sheet.append(row)
    keys=column_dict.keys()
    sample_names="__".join(keys)
    out_prefix=output_prefix+sample_names
    write_spreadsheet(out_prefix+".csv",{1:output_sheet})
    with open(out_prefix+"_dict.json","w") as f:
     json.dump(sample_sampling_dict,f)
