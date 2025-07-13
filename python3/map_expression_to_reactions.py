import getopt
import sys
 
import cobra
import os
import copy
import math
from cobrafunctions.read_spreadsheets import read_spreadsheets
from cobrafunctions.write_spreadsheet import write_spreadsheet
from cobrafunctions.gim3e import get_gene_exp
from cobrafunctions.qMTA import read_gene_data,find_nearZeroVar
import pandas


"""
map_expression_to_reactions.py
This script maps genetically imputed patient-specific expression patterns to organ-specific models using the gene reaction annotations in these models. Imputed values are expressed as Log2 fold changes relative to average gene expression in a given organ and then mapped to reactions in the organ-specific model considering the relative transcript abundance of isoenzymes and enzyme subunits (e.g. in a reaction catalysed by multiple isoenzymes genetic variation on the isoenzyme with the highest expression will have a stronger effect on putative reaction activity). The script must be run for each organ under study. 
Usage: map_expression_to_reactions.py [INPUTS...] 

INPUTS:
-i, --imputed_transcript_abundance : Path to the CSV or XLSX file with the Organ-specific transcript abundance patterns imputed from genotype data. Genes must be in rows and individuals/samples in columns. 
-m, --organ_specific_model : path to the organ-specific model in SBML format. Can also take the gim3e__constrained_model.sbml model as input. 
-r, --reference_transcript_abundance :  Path to the CSV or XLSX file defining the average organ gene expression in TPM or FPKM.
-o, --output_directory : Working and output directory. Will be created if it does not exist. 
-s, --sample_list : Optional, list of sample/individual IDs that should be analysed. If not provided all samples will be analysed. Each row should contain a sample/individuals ID.
-t, --organ_name :  Optional, Organ or tissue to be analysed. Has to match a column in the reference_transcript_abundance file. If not provided it will take organ name from the  organ_specific_model file name
-g, --gene_id_column_name* : Optional, defines the column name in imputed_transcript_abundance that defines the gene identifiers used in the model. If it is not provided, it will be set to "NCBI.gene..formerly.Entrezgene..ID". If it is not present it will be assumed to be the first column in the file. 

OUTPUTS:
reaction_expression: CSV file containing putative reaction activity fold changes for each individual. Used as input for run_qMTA.py.

"""
tissue_key_defined_flag=False
sample_output=""
sample_list_file=None
gene_id_str="gene_id"
remove_reactions_with_low_variance=True
gz=False
opts, args = getopt.getopt(sys.argv[1:],"t:i:m:r:s:o:g:kz",["organ_name=","imputed_transcript_abundance=","organ_specific_model=","reference_transcript_abundance=","sample_list=","output_directory=","gene_id_column_name=","keep_reactions_with_low_variance","gz"])
for opt, arg in opts:
      #print opt,arg
      if opt in ("-i", "--imputed_transcript_abundance"):
          imputed_file= arg
      elif opt in ("-s", "--sample_list"):
         sample_list_file = arg
         sample_output="_"+os.path.basename(sample_list_file).replace(".csv","").replace(".xlsx","")
      elif opt in ("-r", "--reference_transcript_abundance"):
         gene_expression_file = arg
      elif opt in ("-m", "--organ_specific_model"):
          sbml_file= arg
          if not tissue_key_defined_flag: #If tissue key has been manually defined, don't take it from here
             tissue_key=os.path.basename(sbml_file).replace(".sbml","").replace(".xml","").replace(".SBML","").replace(".XML","").replace("organ_metabolic_network_","").replace("organ_specific_metabolic_network_","").replace("_gim3e__restricted_model","").replace("_gim3e__constrained_model","")
      elif opt in ("-o", "--output_directory"):
          output_folder= arg
          if not os.path.exists(output_folder):
             os.makedirs(output_folder)
          os.chdir(output_folder)
      elif opt in ("-t", "--organ_name"):
          tissue_key= arg
          tissue_key_defined_flag=True 
      elif opt in ("-g", "--gene_id_column_name"):
         gene_id_str = str(arg)
      elif opt in ("-k", "--keep_reactions_with_low_variance"):
         remove_reactions_with_low_variance=False
      elif opt in ("-z", "--gz"):
         gz=True


#Build dict with all data
model_dict={tissue_key:{"model":cobra.io.read_sbml_model(sbml_file), "imputed_file":imputed_file}}
conditions_of_interest=list(model_dict.keys())#["Adipose_Subcutaneous"] #Conditions that will be analyzed with GIMME


###########Gene expression inpputs

gpr_mode="full" #By default leave it to full to use full GPR rules. Alternatively, you can set it to the "average" to use to average gene expression for the genes associated to each reaction 
or_mode="sum" #When there are multiple isoenzymes indicated with OR the gene expression of each isoform will be added. ONLY USE SUM IF GENE EXPRESSION DATA IS NOT IN LOG SCALE in gene_expression_file
absent_gene_expression=5 #Percentile of expression of metabolic genes assigned to genes without experimental measurments (Note that in RNA-SEQ most genes should have experimental measures) 
replace_and_with_or=True #In GPR rules and will be replaced by or
##Gene prefix and sufix
#Example: If in the SBML model genes are defined as G_10101.1 were 10101 is gene ID, prefix="G_" , sufix="."
gene_prefix="" 
gene_sufix=""


gene_expression=read_spreadsheets(gene_expression_file)
sheet_name=list(gene_expression.keys())[0]
header=gene_expression[sheet_name][0]

#Get tissue position
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

del(gene_expression)

#differential_gene_file="/rds/user/cf545/hpc-work/results/INTERVAL/imputed_tissue_expression/predicted_expression/Muscle_Skeletal/csv/aggregated/annotated_prediction_muscle.csv"





#1 get base and reaction gene expression
#tissue_key="Muscle_Skeletal"

for tissue_key in conditions_of_interest:
 model=model_dict[tissue_key]["model"]
 if replace_and_with_or:
     """     for reaction in model.reactions:
         reaction.gene_reaction_rule=reaction.gene_reaction_rule.replace("and","or")
     """
     for reaction in model.reactions:
         gene_str=""
         for n_reaction_gene,gene in enumerate(reaction.genes):
             if n_reaction_gene==0:
                gene_str=gene.id
             else:
                gene_str+=" or "+gene.id  
         reaction.gene_reaction_rule=gene_str
 n_col=column_dict[tissue_key]
 reaction_expression_dict_base,value_list_base,expression_dict_base=get_gene_exp(model,absent_gene_expression=absent_gene_expression,percentile=True,file_name=gene_expression_file,gene_method="average",gene_prefix=gene_prefix,gene_sufix=gene_sufix,omit_reflections=True,omit_0=False,gene_value_col=n_col,verbose=False,or_mode=or_mode,expression_dict={},round_reaction_expression=None)
 reaction_list=sorted(reaction_expression_dict_base)
 ##In reactions with 2 or more absent genes, there might be discrepancies between get_gene_exp and the method below. The method below should be more accccurate
 if replace_and_with_or:
        reaction_expression_dict_base={}
        for reaction_id in reaction_list:
            reaction_expression_value=0
            gene_reaction_object=model.reactions.get_by_id(reaction_id)
            for gene in gene_reaction_object.genes:
                if gene.id in expression_dict_base:
                   reaction_expression_value+=expression_dict_base[gene.id]
            reaction_expression_dict_base[reaction_id]=reaction_expression_value
 #2
 #differential_gene_sheet=read_spreadsheets(model_dict[tissue_key]["imputed_file"])
 
 if sample_list_file!=None:
   sample_file=read_spreadsheets(sample_list_file)
   sample_file=sample_file[list(sample_file.keys())[0]]
   #Assume samples are in rows
   samples=[row[0] for row in sample_file]
   differential_gene_data=pandas.read_csv(model_dict[tissue_key]["imputed_file"],usecols=[gene_id_str]+samples)
 else:
    #sheet_name=list(differential_gene_sheet.keys())[0]
    differential_gene_data=pandas.read_csv(model_dict[tissue_key]["imputed_file"])
    names_to_omit=[gene_id_str,"dummy","",None,"V1","Row.names","Gene.stable.ID","NCBI.gene..formerly.Entrezgene..ID","Gene.type","Gene.name","Gene.description","pathway","metabolic","dummy","tissue"]
    samples=[x for x in differential_gene_data.columns if x not in names_to_omit]  
 rows=[]#[["id"]+[x.replace("feeding_","").replace("fasting_","") for x in reaction_list]]
 differential_gene_data.set_index(gene_id_str,inplace=True)
 for n,key in enumerate(samples):
  try:
    print(str(n)+" "+key)
    reference_model=model
    log2_str=key
    #up_genes, down_genes, log2fold_change_dict,   p_value_dict ,  gene_weight_dict, gene_weight_normalized_dict,=read_gene_data(fname=None,model=reference_model,log2_str=log2_str,log2_factor=1,padj_str="ignored",p_th=1,log2fc_th=0,gene_str=gene_id_str,p_weight_formula="1",sheet_dict=differential_gene_sheet,ignore_p_value=True)
    #Modify genes according to the fold changedifferential_gene_sheet
    expression_dict_sample=copy.deepcopy(expression_dict_base)
    #print log2fold_change_dict
    log2fold_change_dict=differential_gene_data[key].to_dict()
    genes_in_both = expression_dict_base.keys() & log2fold_change_dict.keys()
    for gene_id in genes_in_both:
        expression_dict_sample[gene_id]=expression_dict_base[gene_id]*math.pow(2,log2fold_change_dict[gene_id])
    reaction_expression_dict={}
    if replace_and_with_or:
        for reaction_id in reaction_list:
            reaction_expression_value=0
            gene_reaction_object=model.reactions.get_by_id(reaction_id)
            for gene in gene_reaction_object.genes:
                if gene.id in expression_dict_sample:
                   reaction_expression_value+=expression_dict_sample[gene.id]
            reaction_expression_dict[reaction_id]=reaction_expression_value
    else: reaction_expression_dict,value_list,expression_dict=get_gene_exp(model,absent_gene_expression=absent_gene_expression,percentile=True,file_name=gene_expression_file,gene_method="average",gene_prefix=gene_prefix,gene_sufix=gene_sufix,omit_reflections=True,omit_0=False,gene_value_col=n_col,verbose=False,or_mode=or_mode,expression_dict=expression_dict_sample,round_reaction_expression=None)
    #Get fold change to reation expression
    log2fold_change_dict_reactions={}
    row=[]
    for reaction_id in reaction_list:
        if reaction_expression_dict[reaction_id]!=0:
           fc= float(reaction_expression_dict[reaction_id])/float(reaction_expression_dict_base[reaction_id])
           log2fold_change_dict_reactions[reaction_id]=math.log(fc,2)
           row.append(log2fold_change_dict_reactions[reaction_id])
    rows.append(row)
  except Exception as e:
        print(str(key+": "+str(e)))
 #Add dummy row
 #write_spreadsheet(tissue_key+"_reactions_expression_base"+sample_output+".csv",{1:rows})
 #new_row=["dummy"]+[0]*(len(rows[0])-1)
 #rows.append(new_row)
 rows = pandas.DataFrame(rows)
 rows.index = samples       # Setting row names
 rows.columns = reaction_list    # Setting column names
 #Transpose
 rows=rows.transpose()
 #Remove features with near zero variance
 if(remove_reactions_with_low_variance):
  rows_to_drop, rows_to_keep,mask_to_keep=find_nearZeroVar(rows,uniqueCut=10 ,freqCut=95.0/10.0)
  rows=rows[mask_to_keep]
 out_name=tissue_key+"_reactions_expression"+sample_output+".csv"
 if(gz): out_name+=".gz"
 rows.to_csv(out_name,index_label="id")       



        
