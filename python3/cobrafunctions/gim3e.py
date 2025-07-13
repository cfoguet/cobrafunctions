import re
import numpy as np
import copy
from warnings import warn
import math
import cobra
from cobra.flux_analysis.variability import flux_variability_analysis as cobra_flux_variability_analysis
from cobra import Reaction, Metabolite
import traceback
import logging
from .write_spreadsheet import write_spreadsheet
from .read_spreadsheets import read_spreadsheets
import math
from numpy import abs

try:
  cobra_config = cobra.Configuration()
  cobra_config.solver = "glpk"
  #print("glpk set as default solver")   
except:
  print("could not set glpk to default solver")   


def round_sig(x, sig=2):
  if x==0:
    value=0
  else:
     value=round(x, sig-int(math.floor(math.log10(abs(x))))-1)
  return value

def round_up(number,positions):
    exponent=pow(10,positions)
    new_number=math.ceil(number*exponent)/exponent
    """if new_number==number:
       new_number=number+1.0/exponent"""
    return new_number



def round_down(number,positions):
    if number==0.0:
       return 0
    exponent=pow(10,positions)
    return math.floor(number*exponent-0.0001)/exponent
    """if new_number==number:
       new_number=number-1.0/exponent"""
    return new_number


from cobra.flux_analysis.variability import flux_variability_analysis as cobra_flux_variability_analysis

def flux_variability_analysis(model,fraction_of_optimum=0,tolerance_feasibility=1e-6,reaction_list=None):
       fva={}
       if reaction_list!=None:
          if isinstance(reaction_list[0], str):
             reaction_list=[model.reactions.get_by_id(x) for x in reaction_list]
       try:
          pandas_fva=cobra_flux_variability_analysis(model,fraction_of_optimum=fraction_of_optimum,reaction_list=reaction_list)
       except:
          cobra.io.write_sbml_model(model,"failed_model.sbml")
          raise Exception('FVA failed, error model saved as failed_model.sbml')
       for reaction in pandas_fva.index:
           fva[reaction]={"maximum":pandas_fva.loc[reaction]["maximum"],"minimum":pandas_fva.loc[reaction]["minimum"]}
       
       return fva

#from cobra.core.Reaction import Reaction
#from cobra.core.Metabolite import Metabolite
import copy




def convert_to_irreversible_with_indicators(cobra_model,reaction_id_list,metabolite_list, mutually_exclusive_directionality_constraint = False,label_model=None,reactions_with_no_indicators=[]):
    #Function modified from the work by : """Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23."""
    """Will break all of the reversible reactions into two separate irreversible
     reactions with different directions.  This function call modified from
     a version in the core cobra to facilitate the MILP formulation and
     include gene_reaction_rules with the reverse reaction
   
     Arguments:
      cobra_model: A model object which will be modified in place.
      mutually_exclusive_directionality_constraint: Boolean.  If True, turnover 
       reactions are constructed to serve as MILP constraints to prevent loops.
      
     Returns:
      None, cobra_model is modified in place
    
    
    
    
     """
    reactions_to_add = []
    #from cobra.core.Reaction import Reaction
    #from cobra.core import Metabolite
    reactions_to_make_irreversible=[]
    for x in reaction_id_list:
        reactions_to_make_irreversible.append(cobra_model.reactions.get_by_id(x))
    """for x in lexs:
        reactions_to_make_irreversible.append(cobra_model.reactions.get_by_id(x))"""
    
    #If a label model object  is provided make sure all experimentally measured metabolites (or at least one of the metabolites in the pool) is produced
    full_metabolite_list=copy.copy(metabolite_list)
    print(full_metabolite_list)
    if label_model!=None:
       emus=[]
       for condition in label_model.experimental_dict:
           for emu in label_model.experimental_dict[condition]:
               if emu not in emus:
                  emus.append(emu)
       measured_metabolite_dict={}
       
       for emu in emus:
           iso_id=str(label_model.emu_dict[emu]["met_id"])
           #print label_model.id_isotopomer_object_dict
           #isotopomer_object=label_model.id_isotopomer_object_dict[iso_id]
           metabolites=label_model.isotopomer_id_metabolite_id_dict[iso_id]
           print([iso_id,label_model.isotopomer_id_metabolite_id_dict[iso_id]])
           if isinstance(metabolites,list):
              for metabolite in metabolites:
                  full_metabolite_list.append(metabolite)
           else:
              full_metabolite_list.append(metabolites)
    
    for metabolites in full_metabolite_list:
       print(metabolites)
       if not isinstance(metabolites,list):
          metabolites=[metabolites]
       for metabolite in metabolites:
          print(metabolite)
          the_metabolite=cobra_model.metabolites.get_by_id(metabolite)
          for x in the_metabolite.reactions:
             if x not in reactions_to_make_irreversible:
              reactions_to_make_irreversible.append(x)    
                  
    for reaction in reactions_to_make_irreversible:
        # Potential artifact because a reaction might run backwards naturally
        # and this would result in adding an empty reaction to the
        # model in addition to the reverse reaction.
        if reaction.lower_bound < 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            #reverse_reaction = copy.deepcopy(reaction)
            print("adding reverse reaction")
            reverse_reaction.gene_reaction_rule=reaction.gene_reaction_rule
            reverse_reaction.id = reaction.id + "_reverse"
            reverse_reaction.lower_bound = max(0,-1*reaction.upper_bound)
            reverse_reaction.upper_bound = reaction.lower_bound * -1.
            #reaction.lower_bound = 0
            if reaction.upper_bound<0:
               reaction.upper_bound=0
            reaction.lower_bound = 0   
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {}
            current_metabolites = [x for x in reaction.metabolites]
            for the_metabolite in current_metabolites:
                reaction_dict[the_metabolite] = -1 * reaction.get_coefficient(the_metabolite.id)
            reverse_reaction.add_metabolites(reaction_dict)
            reactions_to_add.append(reverse_reaction)
            # Also: GPRs should already copy
            # reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            # reverse_reaction._genes = reaction._genes#
            
            if mutually_exclusive_directionality_constraint and reaction.id not in reactions_with_no_indicators:
                # A continuous reaction bounded by 0., 1.
                # Serves as a source for the indicator metabolites
                tmp_source = Reaction('IRRMILP_direction_constraint_source_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_source.upper_bound = 1.
                tmp_source.lower_bound = 0.
                # The reverse indicator reaction is
                # an integer-valued reaction bounded by 0,1
                # that activates flux to the reverse reaction
                # and deactivates the forward reaction only when it is
                # turned on to 1
                tmp_indicator = Reaction('IRRMILP_reverse_indicator_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_indicator.upper_bound = 1
                tmp_indicator.lower_bound = 0
                tmp_indicator.variable_kind = 'integer'                    
                flux_constraint_forward = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reaction.id)
                flux_constraint_reverse = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reverse_reaction.id)
                flux_constraint_reverse._constraint_sense = 'G'
                flux_constraint_reverse._bound = 0.
                
                tmp_source.add_metabolites({flux_constraint_forward: 1})
                
                tmp_indicator.add_metabolites({flux_constraint_forward: -1,
                                      flux_constraint_reverse: 1})
                """NEW"""                      
                reverse_reaction.add_metabolites({flux_constraint_reverse: -1e-5})
                reaction.add_metabolites({flux_constraint_forward: -1e-5})
                """ENDNEW"""
                                      
                """OLD if reaction.upper_bound != 0:
                        reaction.add_metabolites({flux_constraint_forward: -1./reaction.upper_bound})
                else:
                    # could put 1.01 X the tolerance here,
                    # This is arbitrary.  Use 0.001
                    # since 1000 is a typical upper bound
                    reaction.add_metabolites({flux_constraint_forward: -0.00001})
                if reverse_reaction.upper_bound != 0:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -1./reverse_reaction.upper_bound})
                else:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -0.00001})"""
                reactions_to_add.append(tmp_indicator)
                reactions_to_add.append(tmp_source)
    cobra_model.add_reactions(reactions_to_add)



def get_innactive_reactions(model,optimal_solution_dict,percentile=25,file_name="",gene_method="average",gene_prefix="",gene_sufix="",absent_gene_expression=50,or_mode="max",omit_0=True,gene_value_col=1,log2=False,gpr_mode="full"):
    reaction_expression_dict,value_list,gene_expression_dict=get_gene_exp(model,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,absent_gene_expression=absent_gene_expression,or_mode=or_mode,omit_0=omit_0,omit_reflections=False,gene_value_col=gene_value_col)
    if log2=="reverse":
       for x in reaction_expression_dict:
           reaction_expression_dict[x]=pow(2,reaction_expression_dict[x]) 
       value_list=[pow(2,x) for x in value_list]
    if log2==True:
       print("Log2")
       for x in reaction_expression_dict:
           reaction_expression_dict[x]=math.log(reaction_expression_dict[x]+0.1,2)
       value_list=[math.log(x+0.1,2) for x in value_list]
    
    core_reactions=[x for x in optimal_solution_dict if abs(optimal_solution_dict[x])>1e-8]
    reflection_dict={}
    for x in model.reactions:
        if "reflection" in x.notes:
            reflection_dict[x.id]=x.notes["reflection"]
    low_th=np.percentile(value_list,percentile)
    low_reactions=[x for x in reaction_expression_dict if reaction_expression_dict[x]< low_th]
    low_reactions=[x for x in low_reactions if x not in core_reactions and reflection_dict.get(x) not in core_reactions]
    
    print(low_th, len(low_reactions))
    return low_reactions, core_reactions,low_th,reaction_expression_dict, gene_expression_dict,value_list




def integrate_omics_gim3e_and_remove(metabolic_model,gene_expression_file,fraction_of_optimum=1,low_expression_threshold=25,absent_gene_expression=50,absent_reaction_expression=100,percentile=True,gene_method="average",gene_prefix="",gene_sufix="",metabolite_list_fname=None,label_model=None,epsilon=0.0001,gim3e_fraction_optimum=0.75,run_fva=True,add_as_constraints=True,boundaries_precision=0.001,all_gene_exp_4_percentile=False,base_penalty=0,convert_solution_dict_to_reversible=True,gpr_mode="full",or_mode="max",omit_0=True,gene_value_col=1,log2=False,reaction_list=[],penalty_mode="normal",abs_max_objective=None,lp_tolerance_feasibility=1e-6,remove_reactions_below=33,remove_percentile=True,gene_expression_model=None,penalty_precision=4,reaction_ids_to_omit=[],reactions_to_keep=[],correct_for_complexes=False):
   """
   Functions that reads gene expression data and optionally qualtitaive metabolomics data and uses the GIM3E algorythm to find the flux distribution that is most consistent with such data. Optionally,you can use the results to constraint the cobra model.
   metabolic_model
             Model where the gene expression data and metabolomics data will be integrated
   fraction_of_optimum: float
             Fraction of the model objective reactions that must be fulfilled
   low_expression_threshold: float
             Threshold at which a gene is considered lowly expressed- If percentile is set to True this will be a percentile
   absent_gene_expression: Gene expression level given to genes that are not measured.  If percentile is set to True this will be a percentile
   percentile: bool
              If True the low_expression_threshold and absent_gene_expression will be percintiles
    gene_method: str
          Determines wich gene expression value should be given to a gene when there are multiples entries for the same gene. It can either be "average" (it will use the average of all values), "maximum" (it will use the maximum) or "minimum" (it will use the minimum)
    gene_prefix: str
          Prefix used by genes in the cobra model but not present in the gene expression file
    gene_sufix= str
          Sufix  used by genes in the cobra model but not present in the gene expression file. In Recon 1 alternative transtricpts are indicated by appending _AT1, _AT2 , AT3_ at the end of gene. If in the gene expression file alternative scripts are not indicated in that case "_AT" should be defined as Sufix
   metabolite_list_fname: string, optional
          The path of a XLSX or CSV file that indicates the metabolites that have been detected in the study conditions. Metabolites can be indicated either from metabolite name or metabolite ID. 
   label_model: label_model object, optional
          If a label_model object is provided the algorythm will ensure that all metabolites whose isotopologues are quantified can be produced 
   epsilon: float
          Maximum flux to consider a reaction inaactive. Used for setting the lower bound on turnover reactions
   gim3e_fraction_optimum: float
          Fraction of gim3e optimum that must be fulfilled. Must have a value bewteen 0 and 1
   add_as_constraints: bool
          If True the gim3e results will be used to constraint the input model
   boundaries_precision: float 
          If add_as_constraints is set to True this determines the precision of the constraints that will be added
   """
   if percentile in (True,"true","True",1,"1","yes"):
      percentile=True
   else:
      percentile=False
   low_reactions=[]
   core_reactions=[]   
   precision=int(-1*(math.log10(boundaries_precision)))
   if gene_expression_model==None:
      gene_expression_model=copy.deepcopy(metabolic_model)
      if metabolite_list_fname!=None and metabolite_list_fname!="": 
         metabolite_list=read_metabolomics_data(gene_expression_model,metabolite_list_fname)
      else:
         metabolite_list=[]
      """if len(gene_expression_model.objective)>1:
      raise Exception("Multiple objectives not supported")"""
      for reaction in gene_expression_model.reactions: #TODO Add support for multiple objectives
       if reaction.objective_coefficient!=0:
          fva1=flux_variability_analysis(gene_expression_model,fraction_of_optimum=fraction_of_optimum,reaction_list=[reaction])
          reaction.lower_bound=max(round_down(fva1[reaction.id]["minimum"]-boundaries_precision/2.0,precision),reaction.lower_bound)
          reaction.upper_bound=min(round_up(fva1[reaction.id]["maximum"]+boundaries_precision/2.0,precision),reaction.upper_bound)  
          reaction.objective_coefficient=0
      penalty_dict=create_gim3e_model(gene_expression_model,file_name=gene_expression_file,metabolite_list=metabolite_list,label_model=label_model,epsilon=epsilon,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,low_expression_threshold=low_expression_threshold,absent_gene_expression=absent_gene_expression,absent_reaction_expression=absent_reaction_expression,percentile=percentile,all_gene_exp_4_percentile=all_gene_exp_4_percentile,base_penalty=base_penalty,milp=False,gpr_mode=gpr_mode,or_mode=or_mode,omit_0=omit_0,gene_value_col=gene_value_col,log2=log2,penalty_mode=penalty_mode,penalty_precision=penalty_precision,reaction_ids_to_omit=reaction_ids_to_omit,correct_for_complexes=correct_for_complexes)
   else:
      gene_expression_model=gene_expression_model.copy()
      penalty_dict={}  
   solution=gene_expression_model.optimize()
   objective=-1*solution.objective_value
   solution_dict={}
   irreversible_optimal_solution_dict=solution.fluxes.to_dict()
   reverese_reactions=gene_expression_model.reactions.query("_reverse")
   if convert_solution_dict_to_reversible: 
     for reaction in gene_expression_model.reactions:
        reaction_id=reaction.id
        if reaction_id not in metabolic_model.reactions:
           continue 
        """reaction=gene_expression_model.reactions.get_by_id(reaction_id)
        if reaction in reverese_reactions or "IRRMILP_" in reaction.id:
           continue"""
        solution_dict[reaction_id]=irreversible_optimal_solution_dict[reaction_id]
        if "reflection" in reaction.notes:
            reverse_reaction= reaction.notes["reflection"]
            solution_dict[reaction_id]-=irreversible_optimal_solution_dict[reverse_reaction]
            #print reversible_fva
   else:         
     solution_dict=irreversible_optimal_solution_dict
   if remove_reactions_below>0:
      low_reactions, core_reactions,low_th, reaction_expression_dict, gene_expression_dict,value_list=get_innactive_reactions(gene_expression_model,irreversible_optimal_solution_dict,percentile=remove_reactions_below,file_name=gene_expression_file,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,absent_gene_expression=absent_gene_expression,or_mode=or_mode,omit_0=omit_0,gene_value_col=gene_value_col,log2=log2,gpr_mode=gpr_mode)
      for rid in reactions_to_keep:
          if rid in low_reactions:
             low_reactions.remove(rid)
      for rid in low_reactions:
          gene_expression_model.reactions.get_by_id(rid).bounds=(0,0)
   else:
       low_reactions, core_reactions,low_th, reaction_expression_dict, gene_expression_dict,value_list=get_innactive_reactions(gene_expression_model,irreversible_optimal_solution_dict,percentile=remove_reactions_below,file_name=gene_expression_file,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,absent_gene_expression=absent_gene_expression,or_mode=or_mode,omit_0=omit_0,gene_value_col=gene_value_col,log2=log2,gpr_mode=gpr_mode)
       #reaction_expression_dict={}
       #gene_expression_dict={}
       #value_list={} 
       low_th=remove_reactions_below
   print(solution)    
   if run_fva==True:
     if  abs_max_objective!=None:
         gene_expression_model.reactions.get_by_id("gim3e_objective").upper_bound=abs_max_objective
         gim3e_fraction_optimum=None
         print(gene_expression_model.optimize(),"max objective" , abs_max_objective)
         #print "max objective" , abs_max_objective
     fva,irrevfva=flux_minimization_fva(gene_expression_model,solver=None,reaction_list=reaction_list,lp_tolerance_feasibility=lp_tolerance_feasibility,objective="gim3e_objective",fraction_optimum=gim3e_fraction_optimum,reaction_ids_to_omit=reaction_ids_to_omit)
     if add_as_constraints==True:
       for reaction_id in fva: 
         if reaction_id in metabolic_model.reactions:
            if "RGROUP_" in reaction_id:
                continue
            reaction=metabolic_model.reactions.get_by_id(reaction_id)
            lower_bound=fva[reaction_id]["minimum"]#round_down(fva[reaction_id]["minimum"],precision)  
            upper_bound=fva[reaction_id]["maximum"]#round_up(fva[reaction_id]["maximum"],precision)  
            reaction.lower_bound=max(round_down(lower_bound,precision),reaction.lower_bound)
            reaction.upper_bound=min(round_up(upper_bound,precision),reaction.upper_bound)       
            reaction.objective_coefficient=0 #The new constraints enforce the optimalty of the solution, objective coefficient can be
   else:
     fva={}           
   return penalty_dict,gene_expression_model,objective,solution_dict,fva, low_reactions,core_reactions,low_th, reaction_expression_dict, gene_expression_dict, value_list




def create_gim3e_model(cobra_model,file_name="gene_expression_data.xlsx",metabolite_list=[],label_model=None,epsilon=0.0001,gene_method="average",gene_prefix="",gene_sufix="",low_expression_threshold=25,absent_gene_expression=100,absent_reaction_expression=100,percentile=True,all_gene_exp_4_percentile=False,base_penalty=0,milp=True,or_mode="max",gpr_mode="full",omit_0=True,gene_value_col=1,log2=False,penalty_mode="normal", penalty_precision=4,reaction_ids_to_omit=[],correct_for_complexes=False):
    """
    Creates a Gim3e model
    """
    if gpr_mode=="average":
       reaction_expression_dict,value_list,c=get_average_gene_expression(cobra_model,file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,omit_0=omit_0,gene_value_col=gene_value_col) 
    else:    
       reaction_expression_dict,value_list,expression_dict=get_gene_exp(cobra_model,absent_gene_expression=absent_gene_expression,percentile=percentile,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,or_mode=or_mode,omit_0=omit_0,gene_value_col=gene_value_col)
    if log2=="reverse":
       for x in reaction_expression_dict:
           reaction_expression_dict[x]=pow(2,reaction_expression_dict[x]) 
       value_list=[pow(2,x) for x in value_list]
    if log2==True:
       print("Log2")
       for x in reaction_expression_dict:
           reaction_expression_dict[x]=math.log(reaction_expression_dict[x]+0.1,2)
       value_list=[math.log(x+0.1,2) for x in value_list]
    if absent_reaction_expression!=100:
       absent_reaction_expression_value=np.percentile(value_list,absent_reaction_expression)
       for reaction in cobra_model.reactions:
         if reaction.id not in reaction_expression_dict:
            reaction_expression_dict[reaction.id]=absent_reaction_expression_value
            print(reaction_expression_dict)
    #return
    if percentile==True: #If percentile is set to True, low_expression_threshold is assumed to be a percintile
       """if not all_gene_exp_4_percentile:
         value_list=[] 
         for x in expression_dict:
           value_list.append(expression_dict[x])
       else:
         value_list
       """
       low_expression_threshold=lower_percentile=np.percentile(value_list,low_expression_threshold)
    print(low_expression_threshold)
    #f=open("gene_log.txt","w")
    #f.write("low_threshol:"+str(low_expression_threshold)+"\n")
    #f.write(str(log2)+"\n")
    #f.write(str(reaction_expression_dict))
    #f.close()
    penalty_dict={}
    if base_penalty!=0:
        for reaction in cobra_model.reactions:
            penalty_dict[reaction.id]=base_penalty
    for reaction_id in  reaction_expression_dict:
        gene_expression_value=reaction_expression_dict[reaction_id]
        if reaction_id not in penalty_dict:
              penalty_dict[reaction_id]=0
        if penalty_mode=="division":
           penalty_dict[reaction_id]=round(1000/gene_expression_value,3) 
        else:    
         if gene_expression_value< low_expression_threshold:
           if reaction_id not in penalty_dict:
              penalty_dict[reaction_id]=0
           gene_expression_penalty=round(-gene_expression_value+low_expression_threshold,4)
           penalty_dict[reaction_id]+=gene_expression_penalty
    print(penalty_dict)
    if correct_for_complexes: #correct for reactions like '3.0 6pgc_c + 3.0 nadp_c --> 3.0 co2_c + 3.0 nadph_c + 3.0 ru5p__D_c'
       for reaction in cobra_model.reactions:
           min_coef=min([abs(reaction.metabolites[x]) for x in reaction.metabolites])     
           if reaction.id in penalty_dict:
              penalty_dict[reaction.id]*=min_coef 
               
    for rid in reaction_ids_to_omit:
        del(penalty_dict[rid])
    #return     penalty_dict
    convert_to_irreversible_with_indicators( cobra_model,list(penalty_dict.keys()),metabolite_list=[], mutually_exclusive_directionality_constraint = milp,label_model=label_model)
    """if len(metabolite_list)>0:
       add_turnover_metabolites(cobra_model, metabolite_id_list=metabolite_list, epsilon=epsilon,label_model=label_model)"""
    
    objective_reaction = Reaction('gim3e_objective')
    gim3e_indicator = Metabolite('gim3e_indicator',formula='',name='',compartment='')
    objective_reaction.add_metabolites({gim3e_indicator:-1})
    cobra_model.add_reactions([objective_reaction])
    objective_reaction.objective_coefficient=-1
    total_bound=0         
    for reaction_id in penalty_dict:
           gene_expression_penalty=round(penalty_dict[reaction_id],penalty_precision)
           reaction=cobra_model.reactions.get_by_id(reaction_id)
           reaction.add_metabolites({gim3e_indicator:gene_expression_penalty})
           total_bound+=gene_expression_penalty*reaction.upper_bound
           if "reflection" in reaction.notes:
               reflection_id=reaction.notes["reflection"]
               reflection=cobra_model.reactions.get_by_id(reflection_id)
               reflection.add_metabolites({gim3e_indicator:gene_expression_penalty})
               total_bound+=gene_expression_penalty*reflection.upper_bound
    print(total_bound)
    objective_reaction.lower_bound=0.0
    objective_reaction.upper_bound=total_bound
    return penalty_dict


def flux_minimization_fva(model,solver=None,reaction_list=[],lp_tolerance_feasibility=1e-6,objective=None,fraction_optimum=1,reaction_ids_to_omit=[]):
  #with irreversible_model: 
  irreversible_model=model.copy() 
  if objective!=None and fraction_optimum!=None:
     objetctive_reaction=irreversible_model.reactions.get_by_id(objective)
     objetctive_reaction_ub=objetctive_reaction.upper_bound
     objetctive_reaction.objective_coefficient=-1 
     irreversible_model.optimize()
     xval=objetctive_reaction.x
     objetctive_reaction.objective_coefficient=0 
     objetctive_reaction.upper_bound=round_sig(xval/fraction_optimum,6)
  elif objective!=None:
      objetctive_reaction=irreversible_model.reactions.get_by_id(objective)
      objetctive_reaction.objective_coefficient=0
  reaction2test=[]
  irreversible_fva={}
  model_reaction_n_dict={}
  if reaction_list!=[]:
       for reaction_id in reaction_list:
           reaction=irreversible_model.reactions.get_by_id(reaction_id)
           if reaction.bounds==(0,0):
              irreversible_fva[reaction] 
              continue 
           reaction2test.append(reaction_id)  
           if "reflection" in reaction.notes:
               reflection= reaction.notes["reflection"]
               if reflection not in reaction2test:
                  reaction2test.append(reflection)
       
  else:
     for reaction in irreversible_model.reactions:
         if  "IRRMILP_" not in reaction.id:
             reaction2test.append(reaction.id) 
  for reaction in reaction2test:
    if reaction in reaction_ids_to_omit:
          reaction2test.remove(reaction)
  #print reaction2test
  #print len(reaction2test)
  #irreversible_fva=flux_variability_analysis(irreversible_model,reaction_list=reaction2test,fraction_of_optimum=0,solver=solver)
  for n,reaction in enumerate(irreversible_model.reactions):
      model_reaction_n_dict[reaction.id]=n
  normal_fva_list=[]
  counter=0
  print("starting analysis")
  for reaction_id in reaction2test:
      reaction=irreversible_model.reactions.get_by_id(reaction_id)
      if "reflection" in reaction.notes :
               counter+=1
               irreversible_model2=irreversible_model
               reaction=irreversible_model2.reactions.get_by_id(reaction_id)
               reflection_id= reaction.notes["reflection"]
               n=model_reaction_n_dict[reaction_id]
               n_reflection=model_reaction_n_dict[reflection_id]
               reflection=irreversible_model2.reactions.get_by_id(reflection_id)
               with irreversible_model:
                   #Test max
                   irreversible_model2.reactions.get_by_id(reflection_id).objective_coefficient=-100.0 #make the reflection inacive 
                   irreversible_model2.reactions.get_by_id(reaction_id).objective_coefficient=1.0
                   #print irreversible_model2.reactions.get_by_id(reflection_id).objective_coefficient,irreversible_model.reactions.get_by_id(reaction_id).objective_coefficient
                   max_sol=irreversible_model2.optimize()#solver=solver,tolerance_feasibility=lp_tolerance_feasibility)
                   max_value=max_sol.fluxes.to_dict()[reaction_id]
                   #print max_value, irreversible_model2.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict[reflection_id]
                   #solution_dict=irreversible_model2.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict
                   #print [[reaction_id,solution_dict[reaction_id]],[reflection_id,solution_dict[reflection_id]]]
                   #reflection.objective_coefficient=1
                   status1=irreversible_model2.optimize()    
               with irreversible_model:
                    irreversible_model3=irreversible_model
                    irreversible_model3.reactions.get_by_id(reflection_id).objective_coefficient=0 #make the reflection inacive 
                    irreversible_model3.reactions.get_by_id(reaction_id).objective_coefficient=-1.0
                    #print irreversible_model.optimize()
                    #print irreversible_model3.reactions.get_by_id(reflection_id).objective_coefficient,irreversible_model.reactions.get_by_id(reaction_id).objective_coefficient
                    min_sol=irreversible_model3.optimize()#irreversible_model3.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility)
                    min_value=min_sol.fluxes.to_dict()[reaction_id]
                    #print min_value, irreversible_model2.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict[reflection_id]
                    #solution_dict=irreversible_model3.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict
               #print [[reaction_id,solution_dict[reaction_id]],[reflection_id,solution_dict[reflection_id]]]
               irreversible_fva[reaction_id]={"minimum":min_value,"maximum":max_value}
               print(counter,reaction_id,irreversible_fva[reaction_id],max_sol.status,min_sol.status)
               #print "sdfghjk",reaction_id,max_value,min_value
               #print irreversible_fva
               #Restore
               #irreversible_model.reactions.get_by_id(reflection_id).objective_coefficient=0
               #irreversible_model.reactions.get_by_id(reaction_id).objective_coefficient=0
      else:
        normal_fva_list.append(reaction_id)
  if normal_fva_list!=[]:
     print(irreversible_model.optimize())
     normal_fva=flux_variability_analysis(irreversible_model,reaction_list=normal_fva_list,fraction_of_optimum=0.0)
     """try:
       print("running FVA")
       normal_fva=flux_variability_analysis(irreversible_model,reaction_list=normal_fva_list,fraction_of_optimum=0.0)
       #normal_fva=flux_variability_analysis(irreversible_model,reaction_list=normal_fva_list,fraction_of_optimum=0.0,tolerance_feasibility=lp_tolerance_feasibility)
     except:
       print("using cplex",  lp_tolerance_feasibility)
       normal_fva=flux_variability_analysis(irreversible_model,reaction_list=normal_fva_list,fraction_of_optimum=0.0,tolerance_feasibility=lp_tolerance_feasibility,solver="cplex")""" #changed
     for reaction_id in normal_fva:
         irreversible_fva[reaction_id]=normal_fva[reaction_id]
  reversible_fva={}
  reverese_reactions=irreversible_model.reactions.query("_reverse") 
  #print irreversible_fva
  for reaction_id in irreversible_fva:
        reaction=irreversible_model.reactions.get_by_id(reaction_id)
        if reaction in reverese_reactions or "IRRMILP_" in reaction.id:
           continue
        reversible_fva[reaction.id]=copy.deepcopy(irreversible_fva[reaction.id])
        if "reflection" in reaction.notes:
            reverse_reaction= reaction.notes["reflection"]
            if irreversible_fva[reverse_reaction]["maximum"]!=0:
               reversible_fva[reaction.id]["minimum"]=-irreversible_fva[reverse_reaction]["maximum"]
            if irreversible_fva[reverse_reaction]["minimum"]!=0: 
               reversible_fva[reaction.id]["maximum"]=-irreversible_fva[reverse_reaction]["minimum"] 
  #print reversible_fva
  #objetctive_reaction.upper_bound=objetctive_reaction_ub  
  return reversible_fva,irreversible_fva  


def get_expression(model,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix="",omit_0=False,gene_value_col=1,verbose=False):
    """
    Reads the gene expression file
    model: cobra model object
    file_name: str
          Name of the file with the gene expression data. It should be either a CSV or a XLSX file with gene names in the first column and gene expression value in the second column
    gene_method: str
          Determines wich gene expression value should be given to a gene when there are multiples entries for the same gene. It can either be "average" (it will use the average of all values), "maximum" (it will use the maximum) or "minimum" (it will use the minimum)
    gene_prefix: str
          Prefix used by genes in the cobra model but not present in the gene expression file
    gene_sufix= str
          Sufix  used by genes in the cobra model but not present in the gene expression file. In Recon 1 alternative transtricpts are indicated by appending _AT1, _AT2 , AT3_ at the end of gene. If in the gene expression file alternative scripts are not indicated in that case _AT should be defined as Sufix
    """
    genexpraw_dict={}
    spreadsheet_dict=read_spreadsheets(file_names=file_name,csv_delimiter=',')
    gene_expression_dict=gene_expression_dict={}
    #wb = load_workbook(file_name, read_only=True)
    #ws=wb.active
    value_list=[]
    for sheet in spreadsheet_dict:
      for n_row,row in enumerate(spreadsheet_dict[sheet]):
       try:   
        #geneid=str(row[0])
        if str(row[0])==None or str(row[0])=="":
             continue
        if(n_row==0):
           if(isinstance(row[gene_value_col], str)):
             print(row[gene_value_col])
             continue    
        #print row
        gene_expression_value=float(row[gene_value_col])
        if gene_expression_value==None:
           continue 
        geneids=str(row[0]).split("///")
        
        #value_list.append(gene_expression_value) 
        for geneid in geneids:
          geneid=geneid.rstrip().lstrip()
          if geneid=="NA":
             continue   
          gene_matches=model.genes.query(geneid)
          regular_expression=""
          if gene_prefix=="":
             regular_expression+="^"
          else:
             regular_expression+=gene_prefix 
          regular_expression+=geneid
          if gene_sufix=="":
             regular_expression+="$"
          elif gene_sufix==".":
             regular_expression+="\." 
          else:
             regular_expression+=gene_sufix
          if verbose:
             print(regular_expression)
          for gene in gene_matches:
            if re.search(regular_expression,gene.id)!=None: 
               if gene.id not in genexpraw_dict:
                  genexpraw_dict[gene.id]=[]
               genexpraw_dict[gene.id].append(gene_expression_value)
       except Exception as e:
          error=traceback.format_exc()
          logging.error(error)
        
    for geneid in genexpraw_dict:
        list_of_values=genexpraw_dict[geneid]
        if gene_method=="average" or gene_method=="mean":
           value=np.mean(list_of_values)
        elif gene_method in("max","maximum"):
           value=max(list_of_values)
        elif gene_method in ("min","minimum"):
           value=min(list_of_values)
        else:
            print("Error: Method not supoprted")
            return {}
        gene_expression_dict[geneid]=value#round(value,4)
        if not(omit_0 and value==0):
           value_list.append(value)      
    print(gene_expression_dict)      
    return value_list, gene_expression_dict


def get_average_gene_expression(model,file_name,gene_method="average",gene_prefix="",gene_sufix="",omit_0=True,gene_value_col=1,verbose=True):
    value_list,expression_dict=get_expression(model,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,omit_0=omit_0,gene_value_col=gene_value_col,verbose=gene_value_col)
    reaction_expression_dict={} #We will make average
    for reaction in model.reactions: 
     for gene in reaction.genes: 
        if gene.id in expression_dict:
           if reaction.id not in reaction_expression_dict:
              reaction_expression_dict[reaction.id]=[]
           reaction_expression_dict[reaction.id].append(expression_dict[gene.id])
    
    average_expression_dict={x : np.average(reaction_expression_dict[x]) for x in reaction_expression_dict}
    return average_expression_dict,value_list, reaction_expression_dict 
           




def get_gene_exp(model,absent_gene_expression=50,percentile=True,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix="",omit_reflections=True,omit_0=False,gene_value_col=1,verbose=True,or_mode="max",expression_dict={},reactions_to_analyze=None,round_reaction_expression=4):
    """
    Assigns each reaction a expression value based on the gene_expression file and the GPR rules
    """
    value_list=[]
    if expression_dict=={}:
      print("loading expression")  
      value_list,expression_dict=get_expression(model,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,omit_0=omit_0,gene_value_col=gene_value_col,verbose=verbose)
      if percentile==True: #If percentile is set to True, low_expression_threshold and high_expression_threshold
          
         absent_gene_expression=np.percentile(value_list,absent_gene_expression)
         #f=open("gene_log.txt","a")
         #f.write("absent_threshol:"+str(absent_gene_expression)+"\n")
         #f.close()
      
      
      for gene in model.genes:
        if gene.id not in expression_dict:
           expression_dict[gene.id]=absent_gene_expression
    reaction_expression_dict={}
    #f=open("ReactionExpression","w")
    if reactions_to_analyze in (None,"",False):
       reaction_list=model.reactions
    else:
        reaction_list=[model.reactions.get_by_id(x) for x in reactions_to_analyze]
    for the_reaction in reaction_list:
        genes_present=[x.id for x in the_reaction.genes if x.id in expression_dict ]
        if len(genes_present)==0:
            continue
        if "reflection" in the_reaction.notes:
            if the_reaction.notes["reflection"] in reaction_expression_dict and omit_reflections:
               continue
        if the_reaction.gene_reaction_rule != "" :
           if verbose:
              print(the_reaction.id,"-------------------------------------")
              print(the_reaction.gene_reaction_rule) 
           the_gene_reaction_relation = copy.deepcopy(the_reaction.gene_reaction_rule)
           for the_gene in the_reaction.genes: 
               the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
               the_gene_reaction_relation = the_gene_re.sub(str(expression_dict[the_gene.id]), the_gene_reaction_relation) 
           #print the_gene_reaction_relation
           expression_value=evaluate_gene_expression_string( the_gene_reaction_relation,or_mode=or_mode)
           if(not round_reaction_expression is None): 
                  expression_value=round(expression_value,round_reaction_expression)
           reaction_expression_dict[the_reaction.id]=expression_value
           #f.write(the_reaction.id+" "+str(expression_value)+"\n")
           #print(the_reaction.id+" "+str(expression_value))    
    #f.close()
    return (reaction_expression_dict,value_list,expression_dict)


def gene_exp_classify_reactions(model,low_expression_threshold=25,high_expression_threshold=75,percentile=True,absent_gene_expression=50,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix="",omit_0=False,or_mode="max",debug=False,gene_value_col=1,gpr_mode="full"):
    if gpr_mode=="average":
       reaction_expression_dict,value_list,c=get_average_gene_expression(model,file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,omit_0=omit_0,gene_value_col=gene_value_col) 
    else:
       reaction_expression_dict,value_list,expression_dict=get_gene_exp (model,absent_gene_expression=absent_gene_expression,percentile=percentile,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,omit_0=omit_0,or_mode=or_mode,gene_value_col=gene_value_col)
    if percentile==True: #If percentile is set to True, low_expression_threshold and high_expression_threshold
       high_expression_threshold=upper_percentile=np.percentile(value_list,high_expression_threshold)
       low_expression_threshold=lower_percentile=np.percentile(value_list,low_expression_threshold)
    print([high_expression_threshold,low_expression_threshold])
    hex_reactions=[]
    lex_reactions=[]
    iex_reactions=[] #Inconclusive reactions
    for reaction_id in reaction_expression_dict:
           expression_value=reaction_expression_dict[reaction_id]
           if expression_value > high_expression_threshold: 
              hex_reactions.append(reaction_id) 
           elif expression_value < low_expression_threshold:      
                lex_reactions.append(reaction_id)
           else :      
                iex_reactions.append(reaction_id)
    if debug:
       return hex_reactions, lex_reactions, iex_reactions,reaction_expression_dict,expression_dict             
    return (hex_reactions, lex_reactions, iex_reactions )  


number_finder = re.compile("[\d]+\.?[\d]*")

#TODO add and mode
class gene_expression_score:
    #Based on the work by : #Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23.
    def __init__(self, value,or_mode="max",added_values=None,verbose=False):
        self.str_value = str(value)
        self.value = float(value)
        self.or_mode=or_mode.lower()
        self.verbose=verbose
        if added_values==None:
           self.added_values=set()
        else:
           self.added_values=added_values 
    
    def __add__(self, other):
        # addition is like AND
        #self.added_values=None #Reset when taking an and
        #New get the list of added_values only for the minim values
        if self.value>other.value:
           added_values=other.added_values
        elif self.value<other.value: 
           added_values=self.added_values
        else: # if both values are equal
           added_values=self.added_values
           added_values.update(other.added_values)
        
        return gene_expression_score(min(self.value, other.value),self.or_mode,added_values)
    
    def __mul__(self, other):
        # multiplication is like OR
        if self.or_mode=="sum":
           if self.value==other.value: #If they are isoforms with the same gene expression value don't add them. With RNA SEQ having two different genes with the same expression value is highly unlikely
              value= self.value
              if self.verbose:
               print(value, "equal")
           elif other.value in self.added_values:
               #print other.value#, self.added_values
               value= self.value
               if self.verbose:
                 print(value, "already present", other.value)
           else:
               value=self.value+other.value
               self.added_values.update([self.value,other.value])
               if self.verbose:
                print(self.value,"+", other.value)#, self.added_values
           return gene_expression_score(value,self.or_mode,self.added_values) 
        else:    
           return gene_expression_score(max(self.value, other.value),self.or_mode,None)
    
    def __neg__(self): #CFC Added
        return gene_expression_score (-self.value,self.or_mode,self.added_values) #CFC Added







def evaluate_gene_expression_string(gene_expression_string,or_mode="max"):
    """ penalty string will have:
        * 'or' statements which need to be converted to min
        * 'and' statements which need to be converted to max
    >>> evaluate_penalty("(1 and 2 and 3)")
    max(1, 2, 3)"""
    # if there are no ands or ors, we are done
    gene_expression_string = gene_expression_string.lower()  # don't want to deal with cases
    
    if "and" not in gene_expression_string and "or" not in gene_expression_string:
        return eval(gene_expression_string)
    # we will replace AND and OR with addition and multiplication
    # equivalent to min/max
    #gene_expression_string = gene_expression_string.replace("or", "+").replace("and", "*") #Changed from GIMME
    gene_expression_string = gene_expression_string.replace("or", "*").replace("and", "+")
    # replace the numbers with the custom class which have overloaded AND/OR
    values = [gene_expression_score(i,or_mode) for i in number_finder.findall(gene_expression_string)]
    values_strings = tuple("values[%i]" % i for i in range(len(values)))
    gene_expression_string = number_finder.sub("%s", gene_expression_string)
    gene_expression_string = gene_expression_string % values_strings
    return eval(gene_expression_string).value

