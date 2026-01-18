
#Functions for working with enzyme constrained models

import pandas as pd
import numpy as np
import cobra
import re
import copy

from typing import TYPE_CHECKING, Dict, List

if TYPE_CHECKING:
    from cobra import Model, Metabolite, Reaction


from .cobra_functions import round_sig, get_equation

from cobra.core import Metabolite, Model, Reaction
from cobra.util import solver as sutil



def get_ec_expanded_reaction_mapping(model,reverse_reaction_pattern="_REV",isoenzyme_reaction_pattern="_EXP_\d+",patterns_to_ommit=["^usage_prot_"],verbose=True):
    #The function aims to find all the forward and reverse reactions expanded reactions
    #For example for a reaction A + B <=> C catalyzed by two isoenzymes E1 and E2
    #We will have the following reactions in the model
    # R_EXP_1: A + B + E1 --> C + E1
    # R_REV_EXP_1: C + E1 --> A + B + E1
    # R_EXP_2: A + B + E2 --> C + E2
    # R_REV_EXP_2: C + E2 --> A + B + E2
    #The function will return a mapping like this
    # {"R":{"forward_reactions":["R_EXP_1","R_EXP_2"],"reverse_reactions":["R_REV_EXP_1","R_REV_EXP_2"]},
    
    #Regex for reactions corresponding to isoenzymes e.g. EXP_1, EXP_2
    isoenzyme_reaction_regex = re.compile(isoenzyme_reaction_pattern)
    if verbose:
       print("Isoenzyme reaction pattern regex: "+isoenzyme_reaction_pattern)

    #Regex for reverse reactions
    rev_regex = re.compile(reverse_reaction_pattern)
    if verbose:
       print("Reverse reaction pattern regex: "+reverse_reaction_pattern)
    #Regex for patterns to ommit. Merge them into a single regex
    if len(patterns_to_ommit)>0:
       ommit_pattern_regex = re.compile("|".join(patterns_to_ommit))
       if verbose:
          print("Ommitting reactions matching patterns: "+str(patterns_to_ommit))
    else:
       ommit_pattern_regex = None
    
    #Initialize mapping dictionary
    mapping_dict={}
    #Initialize list of all reverse reactions
    reverse_reactions=[]
    
    
    for reaction in model.reactions:
        base_reaction_id=reaction.id
        reverse_reaction_flag=False
        #Skip reactions matching ommit patterns
        if ommit_pattern_regex is not None:
           if ommit_pattern_regex.search(base_reaction_id):
               continue
        #Get Base Reaction ID by removing isoenzyme and reverse reaction patterns        
        if isoenzyme_reaction_regex.search(base_reaction_id):
            base_reaction_id=isoenzyme_reaction_regex.sub("",base_reaction_id)
        if rev_regex.search(base_reaction_id):
            base_reaction_id=rev_regex.sub("",base_reaction_id)
            reverse_reaction_flag=True
        #Initialize entry in mapping dictionary if not already present
        if base_reaction_id not in mapping_dict:
            mapping_dict[base_reaction_id]={"forward_reactions":[],"reverse_reactions":[]}
        #Add reaction to appropriate list
        if reverse_reaction_flag:
            mapping_dict[base_reaction_id]["reverse_reactions"].append(reaction.id)
            reverse_reactions.append(reaction.id)
        else:
            mapping_dict[base_reaction_id]["forward_reactions"].append(reaction.id)    
    
    return mapping_dict, reverse_reactions


#Remove reactions without net flux
def remove_blocked_reactions_ec_model(model,min_flux=1e-8,expanded_reaction_mapping_dict=None,protein_metabolite_prefix="prot",net_flux_fva=None,test_individual_reactions=True,fva_processes=None,fva_solver_tolerance_feasibility=None,fva_solver_tolerance_optimality=None,verbose=False):
    from .netflux_variability import flux_variability_analysis_net_flux
    print("Removing blocked reactions from EC model with min_flux="+str(min_flux))
    print("Model has "+str(len(model.reactions))+" reactions and "+str(len(model.metabolites))+" metabolites before removing blocked reactions")
    #Start by removing reactions that are already blocked
    blocked_reactions=[x for x in model.reactions if x.bounds==(0,0)]
    print(str(len(blocked_reactions))+" reactions are already blocked")
    model.remove_reactions(blocked_reactions)
    cobra.manipulation.delete.prune_unused_metabolites(model) 
    print("Model has "+str(len(model.reactions))+" reactions and "+str(len(model.metabolites))+" metabolites after removing reactions with 0,0 bounds")
    if expanded_reaction_mapping_dict is None:
           if verbose:
              print("Building default expanded_reaction_mapping_dict with reverse_reaction_pattern=_REV, isoenzyme_reaction_pattern=_EXP_\\d+ and patterns_to_ommit=[^usage_prot_]")
           expanded_reaction_mapping_dict, _ = get_ec_expanded_reaction_mapping(
                  model, 
                  reverse_reaction_pattern="_REV", 
                  isoenzyme_reaction_pattern="_EXP_\\d+",
                  patterns_to_ommit=["^usage_prot_"],
                  verbose=False)
           if test_individual_reactions:
              #You might still want to test individual reactions
              if  verbose:
                  print("Adding individual reactions to the list")
              for reaction in  model.reactions:
                  expanded_reaction_mapping_dict[reaction.id+"__individual"]={"forward_reactions":[reaction.id],"reverse_reactions": []}
    if net_flux_fva==None:
      net_flux_fva=flux_variability_analysis_net_flux( model=model,
         expanded_reaction_mapping_dict=expanded_reaction_mapping_dict,
         flux_list= list(expanded_reaction_mapping_dict.keys()),
         fraction_of_optimum= 0,
         processes= fva_processes,
         solver_tolerance_feasibility= fva_solver_tolerance_feasibility,
         solver_tolerance_optimality= fva_solver_tolerance_optimality)
    else:
       print ("Reusing Net FluxFva")
    #Add reactions without net flux to the list to remove
    reactions_to_remove=[]        
    for net_flux_id in net_flux_fva.index:
        if abs(net_flux_fva.loc[net_flux_id, "maximum"])<min_flux and abs(net_flux_fva.loc[net_flux_id, "minimum"])<min_flux:
           #Get the individual reactions in the net flux
           individual_rids=expanded_reaction_mapping_dict[net_flux_id]["forward_reactions"]+expanded_reaction_mapping_dict[net_flux_id]["reverse_reactions"]
           reactions_to_remove+=individual_rids
           if verbose:
              print(net_flux_id+" "+str(net_flux_fva.loc[net_flux_id, "minimum"])+" "+str(net_flux_fva.loc[net_flux_id, "maximum"])+" "+str(individual_rids))  
    #Make sure there are no duplicate elements 
    reactions_to_remove=list(set(reactions_to_remove))
    print(str(len(reactions_to_remove))+" reactions to remove")
    reaction_objects_to_remove=[model.reactions.get_by_id(x) for x in reactions_to_remove]
    model.remove_reactions(reaction_objects_to_remove)
    #Some usage reactions that could technically be active might have been blocked when reactions without net flux were removed
    orphan_proteins=[]
    for metabolite in model.metabolites:
      if len(metabolite.reactions)<2 and protein_metabolite_prefix in metabolite.id:
         orphan_proteins.append(metabolite)
    if len(orphan_proteins)>0:
       print("Removing "+str(len(orphan_proteins))+" orphan enzymes")
       model.remove_metabolites(orphan_proteins,destructive=True)

    #Remove empty genes
    genes_to_remove=[]
    for gene in model.genes:
      if len(gene.reactions)==0:
        #print(gene)
        genes_to_remove.append(gene)
    if(len(genes_to_remove)>0):
     cobra.manipulation.delete.remove_genes(model,genes_to_remove)            
    #Remove orphan metabolites
    cobra.manipulation.delete.prune_unused_metabolites(model) 


    print("Model has "+str(len(model.reactions))+" reactions and "+str(len(model.metabolites))+" metabolites after removing blocked reactions")      
    return net_flux_fva, reactions_to_remove



"""Add reporter metabolites and reactions to measure net fluxes."""
def add_net_flux_reporter_reactions(
    model: "Model",
    expanded_reaction_mapping_dict: Dict[str, Dict[str, List[str]]],
    reporter_reaction_prefix: str ="netflux__",
    copy_model: bool = True,
    only_add_in_multireaction_fluxes: bool = True,
) -> "Model":
    """Add reporter metabolites and reactions to measure net fluxes.
    
    For each net flux defined in expanded_reaction_mapping_dict, this function:
    1. Creates a reporter metabolite (e.g., "reporter_met_PGI_net")
    2. Adds the metabolite to forward reactions (produced with stoich +1)
    3. Adds the metabolite to reverse reactions (consumed with stoich -1)
    4. Creates a reporter reaction that consumes the metabolite (unbounded)
    
    The flux through each reporter reaction equals the net flux.
    
    WARNING: This function modifies the model structure. It should NOT be used
    if other options are available like get_net_fluxes_from_ec_model or
    flux_variability_analysis_net_flux, which don't require model modification.
    
    Parameters
    ----------
    model : cobra.Model
        The model to add reporter reactions to.
    expanded_reaction_mapping_dict : dict
        A dictionary mapping flux IDs to forward and reverse reactions.
        Format:
        {
           "flux_id1": {
               "forward_reactions": ["reaction_id1", "reaction_id2"],
               "reverse_reactions": ["reaction_id3"]
           },
           "flux_id2": {
               "forward_reactions": ["reaction_id4"],
               "reverse_reactions": []
           },
        }
    copy_model : bool, optional
        Whether to copy the model before modification (default True).
        
    Returns
    -------
    cobra.Model
        The modified model with reporter metabolites and reactions added.
        
    Examples
    --------
    >>> mapping = {
    ...     "PGI_net": {
    ...         "forward_reactions": ["PGI_forward"],
    ...         "reverse_reactions": ["PGI_reverse"]
    ...     },
    ...     "PFK": {
    ...         "forward_reactions": ["PFK"],
    ...         "reverse_reactions": []
    ...     }
    ... }
    >>> model_with_reporters = add_net_flux_reporter_reactions(model, mapping)
    >>> 
    >>> # Now you can optimize and read net fluxes directly
    >>> solution = model_with_reporters.optimize()
    >>> pgi_net_flux = solution.fluxes["reporter_rxn_PGI_net"]
    >>> pfk_net_flux = solution.fluxes["reporter_rxn_PFK"]
    """
    if copy_model:
        model = model.copy()
    
    # Check if model already has reporter reactions
    #existing_reporters = [rxn.id for rxn in model.reactions if rxn.id.startswith(reporter_reaction_prefix)]
    existing_reporters = [met.id for met in model.metabolites if met.id.startswith("reporter_met_")]
    if existing_reporters:
        raise ValueError(
            f"Model already contains reporter metabolites: {existing_reporters}. "
            "Cannot add net flux reporters to a model that already has them."
        )
    reporter_rxns_list=[]
    for net_flux_id in expanded_reaction_mapping_dict.keys():
        forward_reactions = expanded_reaction_mapping_dict[net_flux_id]["forward_reactions"]
        reverse_reactions = expanded_reaction_mapping_dict[net_flux_id]["reverse_reactions"]
        n_reactions=len(forward_reactions)+len(reverse_reactions)
        if only_add_in_multireaction_fluxes and n_reactions==1:
           continue 
        # Create reporter metabolite
        reporter_met_id = f"reporter_met_{net_flux_id}"
        reporter_met = Metabolite(reporter_met_id)
        reporter_met.name = f"Reporter metabolite for {net_flux_id}"
        reporter_met.compartment = "c"  # Default to cytosol
        
        # Add reporter metabolite to forward reactions (produced)
        reporter_ub=0 #Starting Value, will be increase by other reactions.
        reporter_lb=0 #Starting Value, will be increase by other reactions. 
        for forward_rid in forward_reactions:
            rxn = model.reactions.get_by_id(forward_rid)
            rxn.add_metabolites({reporter_met: 1.0})
            reporter_ub+=max(rxn.upper_bound,0)
            reporter_lb+=min(rxn.lower_bound,0) #if a reaction has not been split it will still have a lower bound

        
        # Add reporter metabolite to reverse reactions (consumed)
        for reverse_rid in reverse_reactions:
            rxn = model.reactions.get_by_id(reverse_rid)
            rxn.add_metabolites({reporter_met: -1.0})
            reporter_lb-=max(rxn.upper_bound,0)
        
        # Create reporter reaction that consumes the metabolite
        reporter_rxn_id = reporter_reaction_prefix + net_flux_id
        reporter_rxn = Reaction(reporter_rxn_id)
        reporter_rxn.name = f"Reporter reaction for {net_flux_id}"
        reporter_rxn.lower_bound = reporter_lb  # Sum of forward reactions bounds with some extra marging
        reporter_rxn.upper_bound = reporter_ub   # Sum of reverse reactions bounds with some extra margin
        reporter_rxn.add_metabolites({reporter_met: -1.0})
        reporter_rxns_list.append(reporter_rxn)
    # Add reporter reactions to model
    model.add_reactions(reporter_rxns_list)
    
    return model





def get_net_fluxes_from_ec_model(model,fluxes,output_flux_breakdown=False,ec_expanded_reaction_mapping_dict=None,reverse_reaction_pattern="_REV",isoenzyme_reaction_pattern="_EXP_\d+"):
    #Given a dictionary or pd.series of fluxes for a model with forward and reverse reactions, return a dictionary with net fluxes
    if isinstance(fluxes, pd.Series):
       fluxes = fluxes.to_dict()
    if ec_expanded_reaction_mapping_dict is None:
         ec_expanded_reaction_mapping_dict, reverse_reactions = get_ec_expanded_reaction_mapping(model,reverse_reaction_pattern=reverse_reaction_pattern,isoenzyme_reaction_pattern=isoenzyme_reaction_pattern  )
         

    net_fluxes=[]
    for net_flux_id in ec_expanded_reaction_mapping_dict:
        flux_str=""
        forward_flux=0
        reverse_flux=0
        for forward_r_id in ec_expanded_reaction_mapping_dict[net_flux_id]["forward_reactions"]:
            forward_flux+=fluxes.get(forward_r_id,0)
            flux_str+=forward_r_id+"="+str(round_sig(fluxes.get(forward_r_id,0),3))+";"
        for reverse_r_id in ec_expanded_reaction_mapping_dict[net_flux_id]["reverse_reactions"]:
            reverse_flux+=fluxes.get(reverse_r_id,0)
            flux_str+=reverse_r_id+"="+str(round_sig(fluxes.get(reverse_r_id,0),3))+";"
        net_fluxes.append({"base_rid":net_flux_id,"flux_breakdown":flux_str,"forward_flux":forward_flux,"reverse_flux":reverse_flux,"net_flux":forward_flux - reverse_flux})
    net_fluxes=pd.DataFrame(net_fluxes)
    net_fluxes.set_index('base_rid',inplace=True)
    if not output_flux_breakdown:
       #Drop flux breakdown column
       net_fluxes=net_fluxes.drop(columns=["flux_breakdown"])
    return net_fluxes

"""
def protein_usage_pfba(ec_model,enzyme_kcat_scaling_factor_dict={},gene_weight_dict={},enzyme_list=[],verbose=False,pfba_fraction_of_optimum=1/0.999):
    #This will minimize the total protein usage in the model
    #Input must be an ec model with protein usage reactions
    #If enzyme list is empty all enzymes in the model will be used
    # enzyme_kcat_scaling_factor_dict is a dictionary with the kcat scaling factor for each enzyme used when building the ec model
    #gene_weights can be used to provide additional weights to genes (for instance based on expression). High values means that the enzyme will be less likely to be used
    #genes are taken from the gene associated to the protein usage reaction
    #Note that because our objectibe is a minimization , fraction of optimum must be >1
    model_pfba=ec_model.copy()
    if len(enzyme_list)==0:
           enzyme_list=[x.id.replace("usage_prot_","") for x in model_pfba.reactions.query("usage_prot_")]
    if verbose:
           print("N Minimized Enzymes: "+str(len(enzyme_list)))
           print("Setting Objectives")
    objective_dict={}
    for enzyme in enzyme_list:
            reaction=model_pfba.reactions.get_by_id("usage_prot_"+enzyme)
            #Scaled_enzyme_usage is -1/scaled_kcat*E=-1/(kcat/enzyme_kcat_scaling_factor)*E=-1*enzyme_kcat_scaling_factor/kcat*E  
            # Hence real enzyme usage is Scaled_enzyme_usage/enzyme_kcat_scaling_factor 
            #Reactions wth higher kcat will have a lower objective coefficient 
            kcat_scaling_factor=enzyme_kcat_scaling_factor_dict.get(enzyme,1) #If not found assume 1
            genes=list(reaction.genes)
            if len(genes)!=1:
                raise Exception("Wrong number of genes in "+reaction.id)
            gene_weight=gene_weight_dict.get(genes[0].id,1) #If not found assume 1
            #print(reaction.id,genes[0].id,gene_weight,gene_weight_dict.get(genes[0].id,None))
            #reaction.objective_coefficient=round_sig(-1/kcat_scaling_factor,2) #Minimize total protein usage scaled by kcat
            #We want to minimize the real usage which is scaled by the kcat scaling factor
            #Because usage reactions are defined as negative (prot_A0A0U1RQ18 <--) the objective coefficient must be positive (maximize)           
            objective_dict[reaction]=gene_weight/kcat_scaling_factor #Minimize total protein usage scaled kcat scaling factor and gene weight 
        #Lets scale the objective to avoid very small coefficients
    median_coefficient=np.median([objective_dict[x] for x in objective_dict])
    if verbose:
           print("Median objective coefficient before scaling "+str(median_coefficient))
    for reaction in objective_dict:
            objective_dict[reaction]=min(max(round_sig(objective_dict[reaction]/median_coefficient,4),1e-6),1e6) #Avoid very small or very large coefficients
            #print(reaction.id,objective_dict[reaction])      
    if verbose:
           print("Minimum objective coefficient after scaling "+str(np.min([objective_dict[x] for x in objective_dict])))
           print("Maximum objective coefficient after scaling "+str(np.max([objective_dict[x] for x in objective_dict])))
           #print("Running pFBA")
    model_pfba.objective=objective_dict
    print(model_pfba.optimize().status)
    pfba_solution=cobra.flux_analysis.pfba(model_pfba,fraction_of_optimum=pfba_fraction_of_optimum)
    if verbose:
             print("pFBA status "+pfba_solution.status) 
             print(pfba_solution)
    return pfba_solution, objective_dict
"""
#Replaces protein usage pfba which had some issues
def minimize_protein_usage_fba(ec_model,enzyme_kcat_scaling_factor_dict={},gene_weight_dict={},enzyme_list=[],base_minimization_coefficient=0.1,verbose=False,fraction_of_optimum=1):
    #This will minimize the total protein usage in the model
    #Input must be an ec model with protein usage reactions
    #If enzyme list is empty all enzymes in the model will be used
    # enzyme_kcat_scaling_factor_dict is a dictionary with the kcat scaling factor for each enzyme used when building the ec model
    #gene_weights can be used to provide additional weights to genes (for instance based on expression). High values means that the enzyme will be less likely to be used
    #genes are taken from the gene associated to the protein usage reaction
    #Note that because our objectibe is a minimization , fraction of optimum must be >1
    model_min_enzyme_usage=ec_model.copy()
    #Set objetcive 
    sutil.fix_objective_as_constraint(model_min_enzyme_usage, fraction=fraction_of_optimum)
    
    
    if len(enzyme_list)==0:
           enzyme_list=[x.id.replace("usage_prot_","") for x in model_min_enzyme_usage.reactions.query("usage_prot_")]
    if verbose:
           print("N Minimized Enzymes: "+str(len(enzyme_list)))
           print("Base Minimization Weight:"+str(base_minimization_coefficient))
           print("Setting Objectives")
    objective_dict={}
    for enzyme in enzyme_list:
            reaction=model_min_enzyme_usage.reactions.get_by_id("usage_prot_"+enzyme)
            #Scaled_enzyme_usage is -1/scaled_kcat*E=-1/(kcat/enzyme_kcat_scaling_factor)*E=-1*enzyme_kcat_scaling_factor/kcat*E  
            # Hence real enzyme usage is Scaled_enzyme_usage/enzyme_kcat_scaling_factor 
            #Reactions wth higher kcat will have a lower objective coefficient 
            kcat_scaling_factor=enzyme_kcat_scaling_factor_dict.get(enzyme,1) #If not found assume 1
            genes=list(reaction.genes)
            if len(genes)!=1:
                raise Exception("Wrong number of genes in "+reaction.id)
            gene_weight=gene_weight_dict.get(genes[0].id,1) #If not found assume 1
            #print(reaction.id,genes[0].id,gene_weight,gene_weight_dict.get(genes[0].id,None))
            #reaction.objective_coefficient=round_sig(-1/kcat_scaling_factor,2) #Minimize total protein usage scaled by kcat
            #We want to minimize the real usage which is scaled by the kcat scaling factor
            #Because usage reactions are defined as negative (prot_A0A0U1RQ18 <--) the objective coefficient must be positive (maximize)           
            #objective_dict[reaction]=base_minimization_coefficient+gene_weight/kcat_scaling_factor #Minimize total protein usage scaled kcat scaling factor and gene weight 
            objective_dict[reaction.reverse_variable]=base_minimization_coefficient+gene_weight/kcat_scaling_factor #Minimize total protein usage scaled kcat scaling factor and gene weight 

        #Lets scale the objective to avoid very small coefficients
    median_coefficient=np.median([objective_dict[x] for x in objective_dict])
    if verbose:
           print("Median objective coefficient before scaling "+str(median_coefficient))
    for reaction in objective_dict:
            #if we are mininimizing the reverase reaction it must be negative
            objective_dict[reaction]=-1*min(max(round_sig(objective_dict[reaction]/median_coefficient,3),1e-3),1000) #Avoid very small or very large coefficients
            #print(reaction.id,objective_dict[reaction])      
    if verbose:
           print("Minimum objective coefficient after scaling "+str(np.min([objective_dict[x] for x in objective_dict])))
           print("Maximum objective coefficient after scaling "+str(np.max([objective_dict[x] for x in objective_dict])))
           #print("Running pFBA")
    #model_min_enzyme_usage.objective=objective_dict
    model_min_enzyme_usage.solver.objective.set_linear_coefficients(objective_dict)
    sol=model_min_enzyme_usage.optimize()
    if verbose:
             print("Sol status "+sol.status) 
             print(sol)
    return sol, model_min_enzyme_usage


def get_enzyme_usage_dataframe(model,fluxes,enzyme_kcat_scaling_factor_dict,gene_expression_df,gene_id_column="Gene name",gene_expression_column="tpm_normoxia_mean",only_nonzero_enzymes=True):
   out_data=[]
   reactions_to_evaluate=[x for x in model.reactions if "usage_prot_" in x.id]
   for reaction in reactions_to_evaluate:
         enzyme=reaction.id.replace("usage_prot_","")
         enzyme_metabolite=model.metabolites.get_by_id("prot_"+enzyme)
         scaled_enzyme_usage_flux=-1*fluxes[reaction.id] #Reaction is structure like this prot_A0A0U1RQ18 <-- -1 v so flux is negative
         #Scaled_enzyme_usage is -1/scaled_kcat*E=-1/(kcat/enzyme_kcat_scaling_factor)*E=-1*enzyme_kcat_scaling_factor/kcat*E  
         # Hence real enzyme usage is Scaled_enzyme_usage/enzyme_kcat_scaling_factor 
         real_enzyme_usage_flux=scaled_enzyme_usage_flux/enzyme_kcat_scaling_factor_dict.get(enzyme)
         #Get the max enzyme usage flux possible
         scaled_max_enzyme_usage_flux= -1*reaction.lower_bound #Reaction is structure like this prot_A0A0U1RQ18 <-- -1 v so flux is negative
         real_max_enzyme_usage_flux=scaled_max_enzyme_usage_flux/enzyme_kcat_scaling_factor_dict.get(enzyme)
         #Get Spare capacity
         scaled_spare_capacity=scaled_max_enzyme_usage_flux-scaled_enzyme_usage_flux
         real_spare_capacity=real_max_enzyme_usage_flux-real_enzyme_usage_flux
         spare_capacity_fraction= scaled_spare_capacity/scaled_max_enzyme_usage_flux if scaled_max_enzyme_usage_flux>0 else 0 
         #Add information about the reactions using the enzyme
         enzyme_reactions_with_non_zero_flux = [enzyme_reaction for enzyme_reaction in enzyme_metabolite.reactions if (abs(fluxes[enzyme_reaction.id]) > 1e-7 and enzyme_reaction.id != reaction.id)]
         enzyme_reactions_str_list = []
         for enzyme_reaction in enzyme_reactions_with_non_zero_flux:
             equation_str = get_equation(model, enzyme_reaction.id)
             flux_value = round_sig(fluxes[enzyme_reaction.id], 3)
             enzyme_reactions_str_list.append(enzyme_reaction.id + "(" + str(flux_value) + "): " + equation_str)
         enzyme_reactions_str = ";".join(enzyme_reactions_str_list) if len(enzyme_reactions_str_list) > 0 else ""
         row = {
             "rid": reaction.id,
             "enzyme": enzyme,
             "scaled_enzyme_usage": scaled_enzyme_usage_flux,
             "available_enzyme_scaled": scaled_max_enzyme_usage_flux,
             "spare_capacity_scaled": scaled_spare_capacity,
             "enzyme_kcat_scaling_factor": enzyme_kcat_scaling_factor_dict.get(enzyme),
             "enzyme_usage": real_enzyme_usage_flux,
             "available_enzyme": real_max_enzyme_usage_flux,
             "spare_capacity": real_spare_capacity,
             "relative_spare_capacity": spare_capacity_fraction,
             "enzyme_reactions": enzyme_reactions_str,
             "gene_id": reaction.gene_reaction_rule
         }
         out_data.append(row)
   enzyme_usage_df=pd.DataFrame(out_data)
   if only_nonzero_enzymes:
      enzyme_usage_df=enzyme_usage_df[enzyme_usage_df["enzyme_usage"]>0]
   #enzyme_usage_df=enzyme_usage_df.sort_values('enzyme_usage_flux', ascending=False)
   enzyme_usage_df = pd.merge(enzyme_usage_df,gene_expression_df[[gene_id_column,gene_expression_column]],
      how="left",left_on="gene_id", # column in left df
      right_on=gene_id_column # column in right df
   )
   enzyme_usage_df["ratio_enzyme_usage_to_expression"]=enzyme_usage_df["enzyme_usage"]/(enzyme_usage_df[gene_expression_column]) 
   max_ratio=enzyme_usage_df["ratio_enzyme_usage_to_expression"].max()
   quantile_ratio_95=enzyme_usage_df["ratio_enzyme_usage_to_expression"].quantile(0.95)
   print("Max ratio enzyme usage to expression: "+str(max_ratio))
   print("95th percentile ratio enzyme usage to expression: "+str(quantile_ratio_95))
   enzyme_usage_df=enzyme_usage_df.sort_values('ratio_enzyme_usage_to_expression', ascending=False)
   return(enzyme_usage_df, max_ratio,quantile_ratio_95)

def set_enzyme_usage_bounds_from_gene_expression(model,gene_expression_dict,enzyme_kcat_scaling_factor_dict,enzyme_to_gene_expression_factor=1,reactions_to_omit=[],proteins_to_omit=[]):
    #Set enzyme usage bounds based on gene expression
    #Gene expression is converted to enzyme usage by multiplying by enzyme_to_gene_expression_factor
    #enzyme_kcat_scaling_factor_dict is a dictionary with the kcat scaling factor for each enzyme used when building the ec model
    missing_genes=[]
    #From the reactions to omit get the enzymes to omit
    proteins_to_omit=copy.deepcopy(proteins_to_omit)
    for rid in reactions_to_omit:
        if rid in model.reactions:
           reaction=model.reactions.get_by_id(rid)
           reaction_proteins=[x.id.replace("prot_","") for x in reaction.metabolites if x.id.startswith("prot_")]
           proteins_to_omit+=reaction_proteins
        else:
           print("Reaction "+rid+" to omit not in model") 

    for reaction in model.reactions.query("usage_prot_"):
        enzyme=reaction.id.replace("usage_prot_","")
        if enzyme in proteins_to_omit:
            print("Skipping enzyme "+enzyme+" as it is in the omit list")
            continue
        genes=list(reaction.genes)
        if len(genes)!=1:
            raise Exception("Wrong number of genes in "+reaction.id)
        gene=genes[0].id
        gene_expression=gene_expression_dict.get(gene,None)
        if gene_expression is None or gene_expression<=0:
           missing_genes.append(gene)
           print("No gene expression found for gene "+gene+" associated to enzyme "+enzyme+" so skipping")
           continue
        #Get the max enzyme usage flux possible
        max_enzyme_usage=enzyme_to_gene_expression_factor*gene_expression
        scaled_max_enzyme_usage=max_enzyme_usage*enzyme_kcat_scaling_factor_dict.get(enzyme,1)
        reaction.lower_bound=-1*scaled_max_enzyme_usage #Reaction is structure like this prot_A0A0U1RQ18 <--  so flux is negative
    return missing_genes

def find_lowest_feasible_enzyme_expression_factor(
	model,
	gene_expression_dict,
	enzyme_kcat_scaling_factor_dict,
	min_factor=0,
	initial_ratio_estimate=1,
	tol=1e-6,
    reactions_to_omit=[],
    proteins_to_omit=[],
	verbose=True
):
	"""
	Iteratively find the lowest enzyme_to_gene_expression_factor that gives a feasible solution.
	Returns the lowest feasible factor and the corresponding solution.
	"""
	low = min_factor
	high = initial_ratio_estimate
	best_factor = None
	best_solution = None

	while high - low > tol:
		mid = (low + high) / 2
		test_model = model.copy()
		set_enzyme_usage_bounds_from_gene_expression(
			test_model,
			gene_expression_dict,
			enzyme_kcat_scaling_factor_dict=enzyme_kcat_scaling_factor_dict,
			enzyme_to_gene_expression_factor=mid,reactions_to_omit=reactions_to_omit,proteins_to_omit=proteins_to_omit
		)
		solution = test_model.optimize()
		if verbose:
			print(f"Testing factor: {mid:.6g}, status: {solution.status}")
		if solution.status == "optimal":
			best_factor = mid
			best_solution = solution
			high = mid
		else:
			low = mid

	if best_factor is not None:
		print(f"Lowest feasible enzyme_to_gene_expression_factor: {best_factor}")
		#Run pfba to get the flux distribution
		test_model = model.copy()	
		set_enzyme_usage_bounds_from_gene_expression(
			test_model,
			gene_expression_dict,
			enzyme_kcat_scaling_factor_dict=enzyme_kcat_scaling_factor_dict,
			enzyme_to_gene_expression_factor=best_factor,reactions_to_omit=reactions_to_omit,proteins_to_omit=proteins_to_omit
		)
		best_solution = cobra.flux_analysis.pfba(test_model)
	return best_factor, best_solution
