# -*- coding: utf-8 -*-
import cobra
import pandas as pd
#from cobra.flux_analysis import flux_variability_analysis
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

#try:
#  cobra_config = cobra.Configuration()
#  cobra_config.solver = "glpk"
  #print("glpk set as default solver")   
#except:
#  print("could not set glpk to default solver")   


import cobra.util.solver as sutil
from pprint import pprint

try:
   from cplex import Cplex
except:
   pass 


def round_sig(x, sig=2):
  if x==0:
    value=0
  else:
     value=round(x, sig-int(math.floor(math.log10(abs(x))))-1)
  return value


from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

#from cobra.manipulation.delete import  remove_genes

def list_to_str(input_list,separator="; "):
    if(isinstance(input_list,str)):input_list=eval(input_list)
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


#Formerly remove_innactive
def remove_blocked_reactions(model,remove=True,fva=None,reaction_id_remove=None,min_flux=1e-8):
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
        if abs(fva2remove[reaction_id]["maximum"])<min_flux and abs(fva2remove[reaction_id]["minimum"])<min_flux:
            #print fva2remove[reaction_id]
            reaction_to_remove.append(model.reactions.get_by_id(reaction_id))
    print(len(reaction_to_remove), "to remove")
    if not remove:
       return  fva, [x.id for x in reaction_to_remove]
  
  else:
     reaction_to_remove=[model.reactions.get_by_id(x) for x in reaction_id_remove]
     original_fva={}
  model.remove_reactions(reaction_to_remove) 
  #for reaction in reaction_to_remove:
  #      reaction.remove_from_model()
  #  
  #model.reactions.get_by_id("EX_hdcea(e)").lower_bound=0 #Set the 
  #Remove empty genes and reactions
  genes_to_remove=[]
  for gene in model.genes:
      if len(gene.reactions)==0:
        #print(gene)
        genes_to_remove.append(gene)
  """for gene in genes_to_remove:
        try:
          gene.remove_from_model()  
        except:
           print("Gene "+ gene.id+" could not be removed")
  """
  if(len(genes_to_remove)>0):
     cobra.manipulation.delete.remove_genes(model,genes_to_remove)            
  """metabolites_to_remove=[]
  for metabolite in model.metabolites:
       if len(metabolite.reactions)==0:
          metabolites_to_remove.append(metabolite) 
  for metabolite in metabolites_to_remove:
        metabolite.remove_from_model()"""
  cobra.manipulation.delete.prune_unused_metabolites(model)      
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

#This keeps the output of the function as dict to keep compatibility with old functions
#Should be deprectaed at some point
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



def sampling(model,n=100,processes=6,objective=None,starts=1,return_matrix=False,return_dataframe=False,method="optgp",thinning=100):
    print(method, thinning)
    reaction_ids=[x.id for x in model.reactions]
    if objective!=None:
        print(model.reactions.get_by_id(objective).lower_bound)
    flux_dict_list=[]
    for i in range(0,starts):
       result_matrix = sample(model, n,processes=processes,method=method,thinning=thinning).to_numpy() #Valid methods are optgp and achr. Process is only used in optgp. Thinning (“Thinning” means only recording samples every n iterations) is only used in both
       result_matrix=np.asmatrix(result_matrix)
       if return_matrix or return_dataframe:
            if i==0:
               aggregated_results=result_matrix
            else:
               aggregated_results=np.vstack((aggregated_results,result_matrix)) 
       else:
         for row in result_matrix:
          flux_dict={}
          for n_flux,flux in enumerate(row):
            flux_dict[reaction_ids[n_flux]]=flux
          flux_dict_list.append(flux_dict)
          if objective!=None:
             print(flux_dict[objective])
    if return_dataframe:
        aggregated_results=pd.DataFrame(aggregated_results)
        aggregated_results.columns=reaction_ids
        aggregated_results['sample_n'] = ["sample_"+str(x) for x in range(n)]
        aggregated_results.set_index("sample_n",inplace=True)
        return  aggregated_results       
    elif return_matrix:
       return np.transpose(aggregated_results), reaction_ids
    else:
       return flux_dict_list



    


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

"""
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
"""

def get_equation(model,reaction_id,include_compartment=False,get_compartment_from_met_id=False):
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



def get_metabolites_by_name(model,metabolite_name,compartment=None):
    met_list=[]
    for met in model.metabolites:
        if met.name.lower()==metabolite_name.lower():
           if compartment==None or  compartment==met.compartment:
              met_list.append(met)
    return met_list



def is_transport_reaction(reaction, transport_subsystem="Transport", metabolite_to_test=None, strict=False):
    """
    Determines if a reaction is a transport reaction.
    If strict=True, all metabolites must be transported (i.e., for every metabolite, there is a matching one in a different compartment with the same formula or name).
    Metabolite_to_test can be provided to check if only a specific metabolite is being transported.
    """
    is_transporter = False
    if strict:
        # Check that all metabolites are transported
        for reaction_metabolite1 in reaction.metabolites:
            transported = False
            for reaction_metabolite2 in reaction.metabolites:
                if reaction_metabolite1.id != reaction_metabolite2.id and reaction_metabolite1.compartment != reaction_metabolite2.compartment:
                    if reaction_metabolite1.formula == reaction_metabolite2.formula or reaction_metabolite1.name == reaction_metabolite2.name:
                        transported = True
                        break
            if not transported:
                return False
        return True
    else:
        for reaction_metabolite1 in reaction.metabolites:
            if metabolite_to_test is not None:
                if reaction_metabolite1.id != metabolite_to_test.id:
                    if reaction_metabolite1.compartment != metabolite_to_test.compartment:
                        if (reaction_metabolite1.formula == metabolite_to_test.formula or
                            reaction_metabolite1.name == metabolite_to_test.name):
                            is_transporter = True
            else:
                for reaction_metabolite2 in reaction.metabolites:
                    if reaction_metabolite1.id != reaction_metabolite2.id:
                        if reaction_metabolite1.compartment != reaction_metabolite2.compartment:
                            if (reaction_metabolite1.formula == reaction_metabolite2.formula or
                                reaction_metabolite1.name == reaction_metabolite2.name):
                                is_transporter = True
                                break
            if is_transporter:
                break
        # Use subsystem for any we might have missed unless we only care for a specific metabolite
        if metabolite_to_test is not None and transport_subsystem is not None:
            if reaction.subsystem not in (None, "") and reaction.subsystem.lower() == transport_subsystem.lower():
                is_transporter = True
        return is_transporter
        
        
def is_boundary_reaction(reaction,boundary_subsystem="Exchange/demand reactions",patterns_to_omit=["usage_prot"]):
    is_boundary=False
    if(len(reaction.metabolites)==1):
      is_boundary=True  
    #Use subsystem for any we might have missed
    if reaction.subsystem.lower()==boundary_subsystem.lower():
      is_boundary=True
    for pattern in patterns_to_omit:
        if pattern in reaction.id:
           is_boundary=False
           break
    return(is_boundary) 



def met_explorer(model,met_name,met_compartment,exclude_transporters=True,solve=True,hide0=True,corereactions=[],fva=pd.DataFrame(),pfba_flag=True,flux_solution=pd.Series(),reaction_expression_dict={}):
    metabolite=get_metabolites_by_name(model,met_name,met_compartment)[0]
    formula=metabolite.formula
    print(metabolite.id, metabolite.name, metabolite.formula)
    sol={}
    if len(flux_solution==0):
        sol= flux_solution
    elif solve:
       if pfba_flag:
         sol=pfba(model).fluxes
       else:
         sol=model.optimize().fluxes
    
    for reaction in metabolite.reactions:
       if reaction.id in model.reactions:
        gene_expression=reaction_expression_dict.get(reaction.id)
        if gene_expression==None:
           gene_expression="" 
        is_transporter=False
        if exclude_transporters:
           is_transporter=is_transport_reaction(reaction,metabolite_to_test=metabolite)
        if not is_transporter :
           if reaction.id in fva.index:
                 fva_out=str(fva.loc[reaction.id]["minimum"])+" "+str(fva.loc[reaction.id]["maximum"])
           else:
                 fva_out=""
           if solve:
              if not hide0 or  abs(sol[reaction.id])>1e-6:
                 print(round(sol[reaction.id],5),reaction.id, get_equation(model,reaction.id,True), reaction.id in corereactions, fva_out,gene_expression)
           else:
              print(reaction.id, get_equation(model,reaction.id,True),fva_out,gene_expression, reaction.bounds,sol.get(reaction.id))




def load_constraints(model,constraint_filename,copy_model=False,remove_inactive=False,precision=None):
    if(copy_model):
       model=model.copy()
    constraint_df=pd.read_csv(constraint_filename)
    blocked_reactions=[]
    for row in constraint_df.itertuples(index=True):
        reaction_id=row.reaction_id
        if(reaction_id in model.reactions):
           lower_bound=row.lower_bound
           upper_bound=row.upper_bound
           if precision!=None:
              lower_bound=round_sig(lower_bound,precision)
              upper_bound=round_sig(upper_bound,precision)
           objective_coefficient=row.objective_coefficient
           print(reaction_id,lower_bound,upper_bound,objective_coefficient,get_equation(model,reaction_id,include_compartment=True,get_compartment_from_met_id=False))
           reaction=model.reactions.get_by_id(reaction_id) 
           reaction.lower_bound=lower_bound
           reaction.upper_bound=upper_bound
           if(lower_bound==0 and upper_bound==0):
             blocked_reactions.append(reaction.id)  
           #Only change if its different
           if reaction.objective_coefficient!=objective_coefficient:
              reaction.objective_coefficient=objective_coefficient
    if remove_inactive:
       model.remove_reactions(blocked_reactions) #it gives a warning but nothing we can do about it
       print(len(blocked_reactions),"reactions removed from model")           
    return(model)


#Write Model Functions 

def write_model_to_csv(model,out_prefix=None,include_metabolites=True,include_fva=False,fva_fraction_optimum=1,fva=pd.DataFrame(),include_pfba=True,pfba_solution=pd.Series()):
    if include_pfba:
           if len(pfba_solution)==0:
              print("Running pFBA")  
              pfba_solution=pfba(model).fluxes 
    if include_fva:
           if len(fva)==0:
              print("Running FVA") 
              fva=cobra_flux_variability_analysis(model=model,fraction_of_optimum=fva_fraction_optimum,processes=4,pfba_factor=None,loopless=False)
    out_data_reactions=[]
    for reaction in model.reactions:
        #Get reaction type
        if is_boundary_reaction(reaction,boundary_subsystem="Exchange/demand reactions"):
           reaction_type="Boundary" 
        elif is_transport_reaction(reaction,transport_subsystem="Transport",metabolite_to_test=None):
           reaction_type="Transporter" 
        else:
           reaction_type="Reaction" 
        reaction_dict={}
        reaction_dict["reaction_id"]=reaction.id
        reaction_dict["name"]=reaction.name
        reaction_dict["subsystem"]=reaction.subsystem
        reaction_dict["reaction_type"]=reaction_type
        reaction_dict["stoichiometry_ids"]=reaction.reaction
        reaction_dict["stoichiometry_names"]=get_equation(model,reaction.id,include_compartment=True,get_compartment_from_met_id=False)
        reaction_dict["lower_bound"]=reaction.lower_bound
        reaction_dict["upper_bound"]=reaction.upper_bound
        reaction_dict["objective_coefficient"]=reaction.objective_coefficient
        reaction_dict["gene_reaction_rule"]=reaction.gene_reaction_rule
        if(include_pfba):
           reaction_dict["pfba_solution"]=pfba_solution[reaction.id]
        if(include_fva):
          reaction_dict["fva_minimum"]=fva.loc[reaction.id]["minimum"]
          reaction_dict["fva_maximum"]=fva.loc[reaction.id]["maximum"]
        out_data_reactions.append(reaction_dict)
    out_df_reactions = pd.DataFrame(out_data_reactions)
    if(out_prefix!=None):
      out_df_reactions.to_csv(out_prefix+"_reactions.csv",index=False)
    #Add metabolite description    
    if(include_metabolites):
       out_metabolite_data=[]
       for metabolite in model.metabolites:
           metabolite_dict={}
           metabolite_dict["metabolite_id"]=metabolite.id
           metabolite_dict["name"]=metabolite.name
           metabolite_dict["compartment"]=metabolite.compartment
           metabolite_dict["formula"]=metabolite.formula
           metabolite_dict["charge"]=metabolite.charge
           metabolite_dict["n_reactions"]=len(metabolite.reactions)
           if include_pfba:
              total_flux=0
              for reaction in metabolite.reactions:
                  total_flux=+abs(pfba_solution[reaction.id])
              metabolite_dict["metabolite_total_flux_pfba"]=total_flux
           out_metabolite_data.append(metabolite_dict) 
       out_df_metabolites = pd.DataFrame(out_metabolite_data)
       if(out_prefix!=None):
         out_df_metabolites.to_csv(out_prefix+"_metabolites.csv",index=False)
       return(out_df_reactions,out_df_metabolites)
    else:     
      return(out_df_reactions)    


def get_exchange_reaction_from_metabolite_name(model,metabolite_name,compartment="e"):
    met_object=get_metabolites_by_name(model,metabolite_name,compartment)
    if len(met_object)!=1: raise Exception("Cannot Find Metabolite: "+metabolite_name)   
    met_object=met_object[0]
    exchange_reaction=[]
    for reaction in met_object.reactions:
        if len(reaction.metabolites)==1 and reaction.metabolites[met_object]==-1:
           exchange_reaction.append(reaction)
    if len(exchange_reaction)!=1: raise Exception("Cannot Find Single Exchange Reactions for "+metabolite_name)
    return(exchange_reaction[0])

#Copied from qmta.py
def relax_constraints(model,factor,max_flux):
    for reaction in model.reactions:
        reaction.lower_bound=round_sig(max(min(0,model.reactions.get_by_id(reaction.id).lower_bound*factor),-max_flux),2)
        reaction.upper_bound=round_sig(min(max(0,model.reactions.get_by_id(reaction.id).upper_bound*factor),max_flux),2)
        #Ensure all reactions can potentially carry a minimum a flux
        if reaction.lower_bound<0:
           reaction.lower_bound=min(reaction.lower_bound,-1e-5*factor)
        if reaction.upper_bound>0:
           reaction.upper_bound=max(reaction.upper_bound,1e-5*factor) 


########Function to updade multiple bounds as effciently as possible


def update_reaction_bounds(
    model,
    lower_bounds=None,
    upper_bounds=None,
    copy_model=False,
    use_cplex_direct=False, #Note that when using this option True the bounds appearing in reaction.bounds will not be updated use caution when using this option
    verbose=False
):
    """
    Update reaction bounds for models     
    Can use either optlang (solver-agnostic) or direct CPLEX interface (much faster).
    For reversible reactions split into forward and reverse variables at the solver level:
    - net_flux = forward - reverse
    - Setting bounds requires updating both forward and reverse upper bounds. This is done automatically by cobr
    
    Parameters
    ----------
    model : cobra.Model
        The metabolic model
    lower_bounds : dict or pd.Series, optional
        Dictionary mapping reaction IDs to new lower bounds
    upper_bounds : dict or pd.Series, optional
        Dictionary mapping reaction IDs to new upper bounds
    copy_model: bool
        Create and return a copy of the model
    use_cplex_direct : bool
        If True, use direct CPLEX interface (faster for large batches).
        If False, use optlang (solver-agnostic but slower).
        Default is False.
    verbose : bool
        Print progress
        
    Returns
    -------
    cobra.Model
        Model with updated bounds
        
    Examples
    --------
    # Update using optlang (works with any solver)
    update_reaction_bounds(model, 
                          upper_bounds={'PGI': 100, 'PFK': 50},
                          use_cplex_direct=False)
    
    # Update using CPLEX direct (much faster, requires CPLEX solver)
    model.solver = 'cplex'
    update_reaction_bounds(model,
                          lower_bounds={'PGI': -100},
                          upper_bounds={'PGI': 100, 'PFK': 50},
                          use_cplex_direct=True)
    """
    # Convert Series to dict if needed
    if isinstance(lower_bounds, pd.Series):
        lower_bounds = lower_bounds.to_dict()
    if isinstance(upper_bounds, pd.Series):
        upper_bounds = upper_bounds.to_dict()
    
    if lower_bounds is None:
        lower_bounds = {}
    if upper_bounds is None:
        upper_bounds = {}
    
    # Verify we have something to do
    if len(lower_bounds) == 0 and len(upper_bounds) == 0:
        if verbose:
            print("No bounds to update")
        return model
    if copy_model:
       model=model.copy() 
    # Route to appropriate method
    if use_cplex_direct:
        return _update_reaction_bounds_cplex(model, lower_bounds, upper_bounds, verbose)
    else:
        return _update_reaction_bounds_optlang(model, lower_bounds, upper_bounds, verbose)


def _update_reaction_bounds_optlang(model, lower_bounds, upper_bounds, verbose):
    """
    Update reaction bounds using optlang (solver-agnostic).
    Handles forward/reverse variable splitting automatically.
    """
    # Get all reactions that need updating
    all_rxn_ids = set(lower_bounds.keys()) | set(upper_bounds.keys())
    
    updated_count = 0
    skipped_count = 0
    
    for rxn_id in all_rxn_ids:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            
            # Get new bounds (use current if not specified)
            new_lb = lower_bounds.get(rxn_id, rxn.lower_bound)
            new_ub = upper_bounds.get(rxn_id, rxn.upper_bound)
            
            # Update bounds using cobra's built-in method
            # This automatically handles forward/reverse variable updates
            rxn.bounds = (new_lb, new_ub)
            
            updated_count += 1
            
            if verbose:
                print(f"  {rxn_id}: [{new_lb}, {new_ub}]")        
        except KeyError:
            print(f"Warning: Reaction {rxn_id} not found, skipping")
            skipped_count += 1
    
    print(f"Updated bounds for {updated_count} reactions using default cobrapy bindings")
    if skipped_count > 0:
       print(f"Skipped {skipped_count} reactions (not found in model)")
    
    return model


def _update_reaction_bounds_cplex(model, lower_bounds, upper_bounds, verbose):
    """
    Update reaction bounds using direct CPLEX interface (much faster).
    Handles forward/reverse variable splitting.
    """
    # Verify CPLEX solver
    current_solver = sutil.interface_to_str(model.problem)
    if 'cplex' not in current_solver.lower():
        raise ValueError(
            f"Model is using {current_solver} solver, but use_cplex_direct=True requires CPLEX.\n"
            f"Either:\n"
            f"  1. Set model.solver = 'cplex' before calling this function, or\n"
            f"  2. Use use_cplex_direct=False to use optlang (works with any solver)"
        )
    
    # Access CPLEX problem directly
    lp = model.solver.problem
    
    #if not isinstance(lp, Cplex):
    #    raise RuntimeError(f"Expected CPLEX problem but got {type(lp)}")
    
    # Get all reactions that need updating
    all_rxn_ids = set(lower_bounds.keys()) | set(upper_bounds.keys())
    
    # Build lists for batch updates
    forward_lb_updates = []
    forward_ub_updates = []
    reverse_lb_updates = []
    reverse_ub_updates = []
    
    updated_count = 0
    skipped_count = 0
    
    for rxn_id in all_rxn_ids:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
        except KeyError:
            if verbose:
                print(f"Warning: Reaction {rxn_id} not found, skipping")
            skipped_count += 1
            continue
        
        # Get new bounds (use current if not specified)
        new_lb = lower_bounds.get(rxn_id, rxn.lower_bound)
        new_ub = upper_bounds.get(rxn_id, rxn.upper_bound)
        
        # Get forward and reverse variable indices
        forward_var = rxn.forward_variable
        reverse_var = rxn.reverse_variable
        
        try:
            forward_idx = lp.variables.get_indices(forward_var.name)
            reverse_idx = lp.variables.get_indices(reverse_var.name)
        except:
            print(f"Warning: Could not find forward/reverse variables for {rxn_id}")
            skipped_count += 1
            continue
        
        # Calculate new bounds for forward and reverse variables
        # For net_flux = forward - reverse to be in [lb, ub]:
        
        # Forward variable bounds: [max(0, lb), max(0, ub)]
        new_forward_lb = max(0, new_lb)
        new_forward_ub = max(0, new_ub)
        
        # Reverse variable bounds: [max(0, -ub), max(0, -lb)]
        new_reverse_lb = max(0, -new_ub)
        new_reverse_ub = max(0, -new_lb)
        
        #Add to the list
        forward_lb_updates.append((forward_idx, new_forward_lb))
        forward_ub_updates.append((forward_idx, new_forward_ub))
        reverse_lb_updates.append((reverse_idx, new_reverse_lb))
        reverse_ub_updates.append((reverse_idx, new_reverse_ub))
        
        updated_count += 1
        
        if verbose:
            print(f"  {rxn_id}: [{new_lb}, {new_ub}] -> "
                  f"forward=[{new_forward_lb}, {new_forward_ub}], "
                  f"reverse=[{new_reverse_lb}, {new_reverse_ub}]")    
    # Apply updates using CPLEX batch operations (THIS IS THE KEY SPEEDUP)
    lp.variables.set_lower_bounds(forward_lb_updates)
    lp.variables.set_upper_bounds(forward_ub_updates)
    lp.variables.set_lower_bounds(reverse_lb_updates)
    lp.variables.set_upper_bounds(reverse_ub_updates)    
    total_updates = len(forward_lb_updates) + len(forward_ub_updates) + len(reverse_lb_updates) + len(reverse_ub_updates)
    print(f"Updated bounds for {updated_count} reactions using CPLEX direct interface. \nWARNING: Values in the cobrapy model object have not been updated.")
    if verbose:
              print(f"  ({total_updates} total variable bound updates: "
              f"{len(forward_lb_updates)} forward_lb, {len(forward_ub_updates)} forward_ub, "
              f"{len(reverse_lb_updates)} reverse_lb, {len(reverse_ub_updates)} reverse_ub)")
              if skipped_count > 0:
                 print(f"Skipped {skipped_count} reactions (not found in model)")    
    return model

####
