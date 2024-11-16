# -*- coding: utf-8 -*-
import cobra
from cobra.flux_analysis import flux_variability_analysis
from .read_spreadsheets import read_spreadsheets
from .write_spreadsheet import write_spreadsheet
import numpy as np
import re
import copy
from cobra.flux_analysis import pfba
import math
try:
    from cobra.flux_analysis import sample
except:
    from cobra.sampling import sample

from cobra import Reaction, Metabolite
from cobra.core.gene import parse_gpr, eval_gpr

try:
  cobra_config = cobra.Configuration()
  cobra_config.solver = "glpk"
  #print("glpk set as default solver")   
except:
  print("could not set glpk to default solver")   


import cobra.util.solver as sutil
from pprint import pprint

def round_sig(x, sig=2):
  if x==0:
    value=0
  else:
     value=round(x, sig-int(math.floor(math.log10(abs(x))))-1)
  return value


from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

from cobra.manipulation.delete import  remove_genes

def list_to_str(input_list,separator="; "):
    output=[str(x) for x in input_list]
    output = separator.join(output)
    return(output)

def clean_str(string_to_process,elements_to_remove=["[u'","'","[","]"]):
    for x in elements_to_remove:
        string_to_process=string_to_process.replace(x,"") 
    return(string_to_process)


def reaction_from_string(model,rid,r_string,bounds=None,gpr=""):
     if rid in model.reactions:
        new_reaction=model.reactions.get_by_id(rid)
     else:    
       new_reaction=Reaction(str(rid))
       model.add_reactions([new_reaction])
       new_reaction.build_reaction_from_string(r_string)
     if bounds!=None:
        new_reaction.bounds=bounds
     if gpr!=None:
        new_reaction.gene_reaction_rule=gpr
     return new_reaction


def define_reaction_group(model,reaction_dict,group_reaction_id=None,lower_bound=None,upper_bound=None,objective_coefficient=0):
    new_reaction_id="RGROUP"
    if group_reaction_id!=None:
       if "RGROUP" in group_reaction_id:
           new_reaction_id=group_reaction_id
       else:
           new_reaction_id+="_"+group_reaction_id
    else:
       for reaction_id in reaction_dict:
           new_reaction_id+="_"+reaction_id
    if new_reaction_id in model.reactions:
       model.reactions.get_by_id(new_reaction_id).remove_from_model()
    new_reaction_name="Reaction Group:"
    for reaction_id in reaction_dict:
        if  reaction_dict[reaction_id]>0:
            new_reaction_name+="+"+reaction_id
        else:
            new_reaction_name+="-"+reaction_id
    metabolite = Metabolite("m"+new_reaction_id,formula='',name="mGROUP"+new_reaction_id,compartment='gr')
    group_reaction = Reaction(new_reaction_id)
    group_reaction.name = new_reaction_name
    group_reaction.subsystem = 'Reaction group'
    if upper_bound!=None:
       group_reaction.upper_bound=upper_bound
    group_reaction.add_metabolites({metabolite:-1})
    if objective_coefficient==None:
        group_reaction.objective_coefficient=0
    model.add_reactions([group_reaction])
    group_reaction.objective_coefficient=objective_coefficient
    theoretical_lower_bound=0
    theoretical_upper_bound=0
    for reaction_id in reaction_dict:
        coef=reaction_dict[reaction_id]
        reaction=model.reactions.get_by_id(reaction_id)
        reaction.add_metabolites({metabolite:coef})
        if coef>=0:
           theoretical_upper_bound+=reaction.upper_bound
           theoretical_lower_bound+=reaction.lower_bound
        else:
           theoretical_upper_bound-=reaction.lower_bound
           theoretical_lower_bound-=reaction.upper_bound
    if lower_bound==None:
        group_reaction.lower_bound=min(round_down(theoretical_lower_bound,2),0)
    else:
        group_reaction.lower_bound=lower_bound
    if upper_bound==None:
        group_reaction.upper_bound=max(round_up(theoretical_upper_bound,2),1000)
    else:
        group_reaction.upper_bound=upper_bound
    return group_reaction


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

def remove_innactive(model,remove=True,fva=None,reaction_id_remove=None):
  if reaction_id_remove==None:
    reaction_to_remove=[]
    reaction_to_test=[]
    for reaction in model.reactions:
        if reaction.upper_bound==0 and reaction.lower_bound==0:
           reaction_to_remove.append(reaction)
        else:
           reaction_to_test.append(reaction)
    print(len(reaction_to_remove), "are already blocked")
    if fva==None:
       fva=fva2remove=flux_variability_analysis(model,fraction_of_optimum=0.000,reaction_list=reaction_to_test)
    else:
       fva2remove=fva
    original_fva=fva
    for reaction_id in fva2remove:
        if abs(fva2remove[reaction_id]["maximum"])<1e-8 and abs(fva2remove[reaction_id]["minimum"])<1e-8:
            #print fva2remove[reaction_id]
            reaction_to_remove.append(model.reactions.get_by_id(reaction_id))
    print(len(reaction_to_remove), "to remove")
    if not remove:
       return  fva, [x.id for x in reaction_to_remove]
  
  else:
     reaction_to_remove=[model.reactions.get_by_id(x) for x in reaction_id_remove]
     original_fva={} 
  for reaction in reaction_to_remove:
        reaction.remove_from_model()
    
    #model.reactions.get_by_id("EX_hdcea(e)").lower_bound=0 #Set the 
    #Remove empty genes and reactions
  genes_to_remove=[]
  for gene in model.genes:
      if len(gene.reactions)==0:
        print(gene)
        genes_to_remove.append(gene)
  for gene in genes_to_remove:
        try:
          gene.remove_from_model()  
        except:
           print("Gene "+ gene.id+" could not be removed")  
  metabolites_to_remove=[]
  for metabolite in model.metabolites:
       if len(metabolite.reactions)==0:
          metabolites_to_remove.append(metabolite) 
  for metabolite in metabolites_to_remove:
        metabolite.remove_from_model()
  return original_fva, [x.id for x in reaction_to_remove]


def find_gene_knockout_reactions(cobra_model, gene_list,
                                 compiled_gene_reaction_rules=None):
    """
    Copied from older cobra.py releases
    
    identify reactions which will be disabled when the genes are knocked out

    cobra_model: :class:`~cobra.core.Model.Model`

    gene_list: iterable of :class:`~cobra.core.Gene.Gene`

    compiled_gene_reaction_rules: dict of {reaction_id: compiled_string}
        If provided, this gives pre-compiled gene_reaction_rule strings.
        The compiled rule strings can be evaluated much faster. If a rule
        is not provided, the regular expression evaluation will be used.
        Because not all gene_reaction_rule strings can be evaluated, this
        dict must exclude any rules which can not be used with eval.
    """
    potential_reactions = set()
    for gene in gene_list:
        if isinstance(gene, str):
            gene = cobra_model.genes.get_by_id(gene)
        potential_reactions.update(gene._reaction)
    gene_set = {str(i) for i in gene_list}
    if compiled_gene_reaction_rules is None:
        compiled_gene_reaction_rules = {r: parse_gpr(r.gene_reaction_rule)[0]
                                        for r in potential_reactions}
    return [r for r in potential_reactions
            if not eval_gpr(compiled_gene_reaction_rules[r], gene_set)]



from cobra import Reaction

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



def remove_isoforms_information(model,separator="\."):
    genes_to_delete=[]
    for reaction in model.reactions:
        replace_dict={}
        for gene in reaction.genes:
            gene_match=re.match("^(.+)"+separator, gene.id)
            if gene_match==None:
               continue
            replace_dict[gene.id]=gene_match.group(1)
            print(gene.id+"->"+gene_match.group(1))
        gene_reaction_rule=reaction.gene_reaction_rule
        for gene_id in replace_dict:
               #gene_reaction_rule=gene_reaction_rule.replace("("+gene_id,"("+replace_dict[gene_id]).replace(" "+gene_id," "+replace_dict[gene_id])
               gene_reaction_rule=gene_reaction_rule.replace("("+gene_id,"("+replace_dict[gene_id]).replace(" "+gene_id," "+replace_dict[gene_id])
               gene_reaction_rule=re.sub("^"+gene_id, replace_dict[gene_id], gene_reaction_rule, count=0, flags=0)
               if len(reaction.genes)==1:
                   gene_reaction_rule=gene_reaction_rule.replace(gene_id,replace_dict[gene_id])
               if gene_id not in genes_to_delete:
                  genes_to_delete.append(gene_id)
        print(reaction.gene_reaction_rule)
        print(gene_reaction_rule)
        reaction.gene_reaction_rule=gene_reaction_rule
    genes_to_remove=[]
    print(genes_to_remove)
    for gene in model.genes:
      if len(gene.reactions)==0:
        print(gene)
        genes_to_remove.append(gene)
    for gene in genes_to_remove:
        #print gene.id
        try:
          model.genes.get_by_id(gene.id).remove_from_model()
          #gene.remove_from_model()"""
        except:
          print(gene.id + "could not be removed")



def sampling(model,n=100,processes=6,objective=None,starts=1,return_matrix=False,method="optgp",thinning=100):
    print(method, thinning)
    reaction_ids=[x.id for x in model.reactions]
    if objective!=None:
        print(model.reactions.get_by_id(objective).lower_bound)
    flux_dict_list=[]
    for i in range(0,starts):
       result_matrix = sample(model, n,processes=processes,method=method,thinning=thinning).to_numpy() #Valid methods are optgp and achr. Process is only used in optgp. Thinning (“Thinning” means only recording samples every n iterations) is only used in achr
       result_matrix=np.asmatrix(result_matrix)
       if not return_matrix:
         for row in result_matrix:
          flux_dict={}
          for n_flux,flux in enumerate(row):
            flux_dict[reaction_ids[n_flux]]=flux
          flux_dict_list.append(flux_dict)
          if objective!=None:
             print(flux_dict[objective])
       elif return_matrix:
            if i==0:
               aggregated_results=result_matrix
            else:
               aggregated_results=np.vstack((aggregated_results,result_matrix)) 
    if not return_matrix:
       return flux_dict_list
    else:
       return np.transpose(aggregated_results), reaction_ids


    


def sampling_matrix_get_mean_sd(aggregated_results,reaction_ids,include_absolute_val_stats=False,percentiles=[25,50,75]):
    stat_dict={}
    for n,row in enumerate(aggregated_results): 
        mean=np.mean(row).item()
        std=np.std(row).item()
        percentile_dict={}
        for percentile in percentiles:
            percentile_dict[str(percentile)]=np.percentile(row,percentile,axis=1).item()
            
        stat_dict[reaction_ids[n]]={"mean":mean,"std":std,"max":np.max(row),"min":np.min(row),"percentile":percentile_dict}
        if include_absolute_val_stats:
           percentile_dict_abs={}
           abs_row=np.abs(row)
           for percentile in percentiles:
            percentile_dict_abs[str(percentile)]=np.percentile(abs_row,percentile,axis=1).item()
           mean=np.mean(abs_row).item()
           std=np.std(abs_row).item() 
           stat_dict[reaction_ids[n]]["abs_percentile"]=percentile_dict_abs
           stat_dict[reaction_ids[n]]["abs_mean"]=mean
           stat_dict[reaction_ids[n]]["abs_std"]=std
    return stat_dict


def met_explorer(met_id,model,exclude_transporters=True,solve=True,hide0=True,corereactions=[],fva=None,pfba_flag=True,solution_dict={},reaction_expression_dict={}):
    metabolite=model.metabolites.get_by_id(met_id)
    print(metabolite.name, metabolite.formula)
    fva_out=""
    sol={}
    if solution_dict not in ({},None):
        sol= solution_dict
        print("A")
    elif solve:
       if pfba_flag:
         sol=pfba(model).fluxes.to_dict()
       else:
         sol=model.optimize().fluxes.to_dict()
    
    formula=metabolite.formula
    for reaction in metabolite.reactions:
        gene_expression=reaction_expression_dict.get(reaction.id)
        if gene_expression==None:
           gene_expression="" 
        is_transporter=False
        if exclude_transporters:
           for reaction_metabolite in reaction.metabolites:
               if reaction_metabolite.id!=metabolite.id:
                  if reaction_metabolite.formula==formula: #If the metabolite ID is different but it has the same formula it is a transporter reaction
                     is_transporter=True
        if not is_transporter :
           if fva!=None:
                 fva_out=fva[reaction.id]
           if solve:
              if not hide0 or  abs(sol[reaction.id])>1e-6:
                 print(round(sol[reaction.id],5),reaction.id, reaction.reaction, reaction.id in corereactions, fva_out,gene_expression)
           else:
              print(reaction.id, reaction.reaction,fva_out,gene_expression, reaction.bounds,sol.get(reaction.id))



def remove_innactive(model,remove=True,fva=None,reaction_id_remove=None):
  if reaction_id_remove==None:
    reaction_to_remove=[]
    reaction_to_test=[]
    for reaction in model.reactions:
        if reaction.upper_bound==0 and reaction.lower_bound==0:
           reaction_to_remove.append(reaction)
        else:
           reaction_to_test.append(reaction)
    print(len(reaction_to_remove), "are already blocked")
    if fva==None:
       fva=fva2remove=flux_variability_analysis(model,fraction_of_optimum=0.000,reaction_list=reaction_to_test)
    else:
       fva2remove=fva
    original_fva=fva
    for reaction_id in fva2remove:
        if abs(fva2remove[reaction_id]["maximum"])<1e-8 and abs(fva2remove[reaction_id]["minimum"])<1e-8:
            #print fva2remove[reaction_id]
            reaction_to_remove.append(model.reactions.get_by_id(reaction_id))
    print(len(reaction_to_remove), "to remove")
    if not remove:
       return  fva, [x.id for x in reaction_to_remove]
  
  else:
     reaction_to_remove=[model.reactions.get_by_id(x) for x in reaction_id_remove]
     original_fva={} 
  for reaction in reaction_to_remove:
        reaction.remove_from_model()
    
    #model.reactions.get_by_id("EX_hdcea(e)").lower_bound=0 #Set the 
    #Remove empty genes and reactions
  genes_to_remove=[]
  for gene in model.genes:
      if len(gene.reactions)==0:
        print(gene)
        genes_to_remove.append(gene)
  for gene in genes_to_remove:
        try:
          gene.remove_from_model()  
        except:
           print("Gene "+ gene.id+" could not be removed")  
  metabolites_to_remove=[]
  for metabolite in model.metabolites:
       if len(metabolite.reactions)==0:
          metabolites_to_remove.append(metabolite) 
  for metabolite in metabolites_to_remove:
        metabolite.remove_from_model()
  return original_fva, [x.id for x in reaction_to_remove]


def get_equation(model,reaction_id,include_compartment=False,get_compartment_from_met_id=True):
    reaction=model.reactions.get_by_id(reaction_id)
    reaction_str=" "+reaction.reaction+" "
    #print reaction_str
    for x in reaction.metabolites:
            met_name=x.name
            if get_compartment_from_met_id:
               compartment=re.search(".+\[(.+)\]$",x.id).group(1)
               met_name=re.sub(" \[.+\]$","",met_name).rstrip()
            else:    
               compartment=x.compartment
            if compartment!=None and include_compartment:
              compartment_str="["+compartment+"]"
              #print x.id, compartment_str
            else:
              compartment_str=""
            if x.name==None:
               x.name=x.id
            #print x.id, x.name      
            reaction_str=reaction_str.replace(" "+x.id+" "," "+met_name+compartment_str+" ")
    return reaction_str.strip()
