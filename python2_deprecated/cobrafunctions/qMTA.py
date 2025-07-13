# -*- coding: utf-8 -*-
import cobra
import math
import numpy as np
import sys
import operator
#sys.path.append("/home/carles/colon/code/")

import copy
from warnings import warn

import scipy
from scipy.sparse import dok_matrix
try:
    from sympy import Basic, Number
except:
    class Basic:
        pass

metabolite_ex_dict={}

from cplex import Cplex, SparsePair,SparseTriple
from cplex.exceptions import CplexError
from six import iteritems, string_types
from six.moves import zip
import json
import cobra
from cobra import Reaction, Metabolite
import scipy
from write_spreadsheet import write_spreadsheet
from read_spreadsheets import read_spreadsheets
from cobra_functions import *
import cobra_functions


from cobra.core.solution import LegacySolution

# solver specific parameters for Cplex
parameter_defaults = {'objective_sense': 'maximize',
                      'tolerance_optimality': 1e-9,
                      'tolerance_feasibility': 1e-9,
                      'tolerance_integer': 1e-9,
                      'lp_method': 1,
                      'tolerance_barrier': 1e-9,
                      'verbose': False,
                      'qpmethod': 1}

parameter_mappings = {'lp_method': 'lpmethod',
                      'lp_parallel': 'threads',
                      'qpmethod':'qpmethod',
                      'threads': 'threads',
                      'objective_sense': 'objective_sense',
                      'time_limit': 'timelimit',
                      'iteration_limit': 'simplex.limits.iterations',
                      'tolerance_barrier': 'barrier.convergetol',
                      'tolerance_feasibility': 'simplex.tolerances.feasibility',
                      'tolerance_markowitz': 'simplex.tolerances.markowitz',
                      'tolerance_optimality': 'simplex.tolerances.optimality',
                      'tolerance_integer': 'mip.tolerances.integrality',
                      'MIP_gap_abs': 'mip.tolerances.absmipgap',
                      'MIP_gap': 'mip.tolerances.mipgap'}

variable_kind_dict ={'continuous': "C", 'integer': "I"}
#variable_kind_dict = {'continuous': Cplex.variables.type.continuous, 'integer': Cplex.variables.type.integer} #{'continuous': "C", 'integer': "I"}
status_dict = {'MIP_infeasible': 'infeasible',
               'integer optimal solution': 'optimal',
               'MIP_optimal': 'optimal',
               'MIP_optimal_tolerance': 'optimal',
               'MIP_unbounded':  'unbounded',
               'infeasible': 'infeasible',
               'integer infeasible': 'infeasible',
               'optimal': 'optimal',
               'optimal_tolerance': 'optimal',
               'unbounded': 'unbounded',
               'integer optimal, tolerance': 'optimal',
               'time limit exceeded': 'time_limit'}

def round_sig(x, sig=2):
  if x==0:
    value=0
  else:
     value=round(x, sig-int(math.floor(math.log10(abs(x))))-1)
  return value


def relax_constraints(model,factor,max_flux):
    for reaction in model.reactions:
        reaction.lower_bound=round_sig(max(min(0,model.reactions.get_by_id(reaction.id).lower_bound*factor),-max_flux),2)
        reaction.upper_bound=round_sig(min(max(0,model.reactions.get_by_id(reaction.id).upper_bound*factor),max_flux),2)
        #Ensure all reactions can potentially carry a minimum a flux
        if reaction.lower_bound<0:
           reaction.lower_bound=min(reaction.lower_bound,-1e-5*factor)
        if reaction.upper_bound>0:
           reaction.upper_bound=max(reaction.upper_bound,1e-5*factor) 


def remove_nearZeroVar(rows,uniqueCut=10 ,freqCut=95.0/10.0,rows_to_omit=["dummy"]):
    #Adapted from the nearZeroVar from the R caret package
    #freqCut	the cutoff for the ratio of the most common value to the second most common value
    #uniqueCut the cutoff for the percentage of distinct values out of the number of total samples
    out_rows=[]
    header=[x for x in rows[0] if x not in rows_to_omit]
    rowlen=len(header)-1 #Excluding row names
    print rowlen
    for row in rows:
        #print row
        row_name=row[0]
        if row_name in rows_to_omit:
           out_rows.append(row)
        row_name=row_name.replace("[","_").replace("]","_")
        row_lite=row[1:len(row)] #Remove rowname
        #print len(row_lite)
        #print row_name,row_lite
        #Get unique values and their frequency
        count_dict = {i:row_lite.count(i) for i in row_lite}
        unique_values_percentatge=float(len(count_dict))/float(rowlen)*100
        sorted_tuple = sorted(count_dict.items(), key=operator.itemgetter(1), reverse=True)
        if len(count_dict)>1:
           freq_first_to_second=float(sorted_tuple[0][1])/float(sorted_tuple[1][1])
        else:
           freq_first_to_second=9999999 
        out_str=str(row_name)+" "+str(unique_values_percentatge),str(freq_first_to_second)
        #print  out_str
        if unique_values_percentatge<uniqueCut and freq_first_to_second>freqCut:
           pass
           #print  out_str
        else:
           row_lite.insert(0,row_name) #Add corrected rowname
           out_rows.append(row)
    return out_rows



def get_status(lp):
    status = lp.solution.get_status_string().lower()
    return status_dict.get(status, status)


def get_objective_value(lp):
    return lp.solution.get_objective_value()



def cplex_format_solution(lp, cobra_model, **kwargs):
    from cplex import Cplex, SparsePair
    from cplex.exceptions import CplexError
    status = lp.solution.get_status_string().lower()
    print status
    if status in ('optimal', 'time_limit', 'non-optimal',"integer optimal solution","integer optimal, tolerance"):
        objective_value = lp.solution.get_objective_value()
        x_dict = dict(zip(lp.variables.get_names(),
                          lp.solution.get_values()))
        """for r in cobra_model.reactions:
            if (r.id not in x_dict) and ("_"+r.id in x_dict):
               x_dict[r.id]=x_dict["_"+r.id]
               del x_dict["_"+r.id]
               print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        """
        
        x = lp.solution.get_values()
        # MIP's don't have duals
        if lp.get_problem_type() in (Cplex.problem_type.MIQP, Cplex.problem_type.MILP):
            y = y_dict = None
        else:
            y_dict = dict(zip(lp.linear_constraints.get_names(),lp.solution.get_dual_values()))
            y = lp.solution.get_dual_values()
    else:
        x = y = x_dict = y_dict = objective_value = None
    return LegacySolution(objective_value, x=x, x_dict=x_dict, status=status,y=y, y_dict=y_dict)


def set_parameter(lp, parameter_name, parameter_value):
    if parameter_name == 'objective_sense':
        parameter_value = getattr(lp.objective.sense, parameter_value)
        lp.objective.set_sense(parameter_value)
        return
    elif parameter_name == 'the_problem':
        warn('option the_problem removed')
        return
    elif parameter_name == 'verbose':
        return
        """
        if parameter_value:
            lp.set_log_stream(sys.stdout)
            lp.set_results_stream(sys.stdout)
            lp.set_warning_stream(sys.stderr)
            # If the value passed in is True, it shold be 1. MIP display can
            # be as high as 5, but the others only go up to 2.
            value = int(parameter_value)
            set_parameter(lp, 'mip.display', value)
            set_parameter(lp, 'simplex.display', min(value, 2))
            set_parameter(lp, 'barrier.display', min(value, 2))
        else:
            lp.set_log_stream(None)
            lp.set_results_stream(None)
            lp.set_warning_stream(None)
            set_parameter(lp, 'mip.display', 0)
            set_parameter(lp, 'simplex.display', 0)
            set_parameter(lp, 'barrier.display', 0)
        """
    try:
        cplex_name = parameter_mappings.get(parameter_name, parameter_name)
        cplex_value = parameter_value
        print cplex_name, cplex_value
        param = lp.parameters
        for i in cplex_name.split("."):
            param = getattr(param, i)
        if isinstance(cplex_value, string_types) and \
                hasattr(param.values, cplex_value):
            cplex_value = getattr(param.values, cplex_value)
        param.set(cplex_value)
    except (CplexError, AttributeError) as e:
        raise ValueError("Failed to set %s to %s: %s" %
                         (parameter_name, str(parameter_value), repr(e)))


def _float(value):
    if isinstance(value, Basic) and not isinstance(value, Number):
        return 0.
    else:
        return float(value)


def read_gene_data(fname,model,log2_str="log2FoldChange",log2_factor=1,padj_str="padj",p_th=0.25,log2fc_th=0,gene_str="NCBI.gene.ID",p_weight_formula="-1*math.log(p_value,10)",sheet_dict={},ignore_p_value=False):
    up_genes=[]
    down_genes=[]
    log2fold_change_dict={}
    p_value_dict={}
    milp_weight_dict={}
    milp_weight_list=[] #For normalization
    if sheet_dict=={}:
       sheet_dict=read_spreadsheets(fname)
    sheet=sheet_dict[sheet_dict.keys()[0]]           
    gene_n=0 #If gene is not defined we will assume is the first one
    for n_row,row in enumerate(sheet): 
      if n_row==0: #Header
         for n, element in enumerate(row):
           if element==None:
               continue   
           if log2_str==element:
               log2_n=n
           elif padj_str==element:
                p_adj_n=n
           elif gene_str==element:
                gene_n=n
         continue   
         
      gene_id= str(row[gene_n])
      log2fc=row[log2_n]
      if not ignore_p_value:
         p_value=row[p_adj_n]
      else:
         p_value=0 
      if p_value in ("NA","",None):
         continue
      log2fc=float(log2fc)*log2_factor
      p_value=float(p_value)
      #print  gene_id, log2fc, p_value   
      if gene_id in model.genes:
         #print row
         if abs(log2fc)>log2fc_th and p_value<p_th:
            if log2fc>=0:
               up_genes.append(gene_id) 
            elif log2fc<0:
               down_genes.append(gene_id)
            log2fold_change_dict[gene_id]=log2fc
            p_value_dict[gene_id]=p_value
            milp_weight=eval(p_weight_formula)
            milp_weight_dict[gene_id]=milp_weight
            milp_weight_list.append(milp_weight)    
    
    milp_mean=np.mean([milp_weight_dict[x] for x in milp_weight_dict])
    milp_weight_normalized_dict={x:round(milp_weight_dict[x]/milp_mean,3) for x in milp_weight_dict}
    return  up_genes, down_genes, log2fold_change_dict,   p_value_dict ,  milp_weight_dict, milp_weight_normalized_dict



def build_weight_dicts(signficant_gene_list,cobra_model,max_reactions4gene=10,gene_weight=0.5,non_gene_met_reaction_weight=0.5,gene_weight_dict={},non_gene_met_reaction_weight_dict={},fold_change_dict={},vref_dict={},max_fold_change=99999999999,min_flux4weight=1e-6,normalize_by_scale_genes=True,normalize_by_scale_unchanged_reactions=True,genes_log2_to_lineal=True,precision=6,min_flux_fold_change=1e-9):
    gene_met_reactions=set()
    genes_to_omit=[]
    quadaratic_dict={}
    coefficient_dict={}
    gene_target_dict={}
    reactions2omit=[]
    for gene_id in signficant_gene_list:
       gene=cobra_model.genes.get_by_id(gene_id) 
       if len(gene.reactions)>max_reactions4gene:
          genes_to_omit.append(gene_id)
    
    
    signficant_gene_list_corrected=[x for x in signficant_gene_list if x not in genes_to_omit]
    ################## Signficant weight
    for gene_id in signficant_gene_list_corrected:
        if not gene_id  in cobra_model.genes:
            continue
        fold_change=max(min(max_fold_change,fold_change_dict[gene_id]),-max_fold_change)
        if genes_log2_to_lineal:
                fold_change=pow(2,fold_change)
        for reaction in cobra_model.genes.get_by_id(gene_id).reactions:
            rid=reaction.id
            gene_met_reactions.add(rid)
            vref=vref_dict[rid]
            #Todo add correction to genes catalyzing many reactions
            target_vref=round(vref*fold_change,7)
            weight=float(gene_weight)
            if normalize_by_scale_genes==True:
               flux_factor= max(abs(pow(vref-target_vref,2)),min_flux4weight)
               #Minimal fold change
               flux_factor_min= max(abs(pow(vref*0.999,2)),min_flux4weight)
               #print flux_factor, flux_factor_min
               weight/= max(flux_factor,0)
            """if normalize_by_scale_genes=="debug":
			   weight/= max(abs(target_vref),min_flux4weight) """
            if gene_id in gene_weight_dict:
               weight*=gene_weight_dict[gene_id]
            """weight=round_sig(weight,precision)"""
            if abs(vref)<min_flux_fold_change:
               weight=0 
            if rid not in coefficient_dict or rid not in quadaratic_dict:
               coefficient_dict[rid]=0
               quadaratic_dict[rid]=0
            if gene_id not in gene_target_dict:
               gene_target_dict[gene_id]={}
            coefficient_dict[rid]+=-2*target_vref*weight#-2ab*weight
            quadaratic_dict[rid]+=weight*2 #b2 Multipy by 2
            gene_target_dict[gene_id][rid]=round_sig(target_vref,4)
            if weight==0:
               gene_target_dict[gene_id][rid]=0  
    
    non_gene_met_reactions=[x for x in cobra_model.reactions if x.id not in  gene_met_reactions]
    for reaction in non_gene_met_reactions:
        genes_id=[x.id for x in reaction.genes]
        if len(set(genes_id).intersection(genes_to_omit))>0:
           reactions2omit.append(reaction.id) 
        rid=reaction.id
        vref=vref_dict[rid]
        weight=non_gene_met_reaction_weight
        if normalize_by_scale_unchanged_reactions==True:
           #weight/= max(abs(pow(vref,2)),min_flux4weight)
           weight/= max(abs(vref),min_flux4weight) #This works better than squared
        """if normalize_by_scale_unchanged_reactions=="debug":
           weight/= max(abs(pow(vref,2)),1e-8)"""
           #weight/= max(abs(vref),min_flux4weight) #This works better than squared
        #if normalize_by_scale_unchanged_reactions:   
        if rid in non_gene_met_reaction_weight_dict:
           weight*=non_gene_met_reaction_weight_dict[rid]
        """if weight!=0:
           weight=round_sig(weight,precision)""" 
        coefficient_dict[rid]=-2*vref*weight#-2ab*weight
        quadaratic_dict[rid]=weight*2 #b2/2 Multipy by 2 
    return quadaratic_dict, coefficient_dict, gene_target_dict, signficant_gene_list_corrected,reactions2omit




def create_rMTA_problem_quadratic(cobra_model,quadaratic_dict={}, coefficient_dict={},out_name="qrMTA.lp",**kwargs ):
    """Solver-specific method for constructing a solver problem from   
    milp_weight 0.25 and quadratic_component_weight 0.5 is equal to alpha 0.5 ( 0.5 Quadratic +0.5/2*yf+0.5/2*yr) 
    if quadratic_component_weight_dict or milp_weight_dict is defined, specific reatcions will be modified as milp_weight_dict[x]*milp_weight and quadratic_component_weight*quadratic_component_weight_dict[x]
    if absolute_value is used (which forces absolute value on Vref-epsilon constraint, correct for negative reactions should not be necessary)
    """
    #print "Process parameter defaults"
    
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = parameter_defaults.copy()
        the_parameters.update(kwargs)
        print the_parameters 
    if 'relax_b' in the_parameters:
        relax_b = the_parameters.pop("relax_b")
        warn('need to reimplement relax_b')
        relax_b = False
    else:
        relax_b = False
    
    # Begin problem creation
    print "Create problem"
    lp = Cplex()
    
    
    for k, v in iteritems(the_parameters):
        set_parameter(lp, k, v)
    
    
    objective_coefficients=[]
    for x in cobra_model.reactions:
        if x.id in coefficient_dict:
           obj_coefficent=coefficient_dict[x.id]
        else:
           obj_coefficent=0
        objective_coefficients.append(obj_coefficent)
    """objective_coefficients = [float(x.objective_coefficient)
                              for x in cobra_model.reactions]"""
    lower_bounds = [_float(x.lower_bound) for x in cobra_model.reactions]
    upper_bounds = [_float(x.upper_bound) for x in cobra_model.reactions]
    variable_names = cobra_model.reactions.list_attr("id")
    variable_kinds = [variable_kind_dict[x.variable_kind] for x
                      in cobra_model.reactions]
    #Constraints
    lp.variables.add(obj=objective_coefficients, lb=lower_bounds, ub=upper_bounds, names=variable_names, types=variable_kinds)
    constraint_sense = []; constraint_names = []; constraint_limits = []
    for x in cobra_model.metabolites:
        if "MOMA_wt_" in x.id:
            continue
        constraint_sense.append(x._constraint_sense)
        constraint_names.append(x.id)
        constraint_limits.append(float(x._bound))
     
    the_linear_expressions = []
    # NOTE: This won't work with metabolites that aren't in any reaction
    for the_metabolite in cobra_model.metabolites:
        if "MOMA_wt_" in the_metabolite.id:
            continue
        variable_list = []
        coefficient_list = []
        for the_reaction in the_metabolite._reaction:
            variable_list.append(str(the_reaction.id))
            coefficient_list.append(_float(the_reaction._metabolites[the_metabolite]))
        the_linear_expressions.append(SparsePair(ind=variable_list,val=coefficient_list))
    # Set objective to quadratic program
    for n,x in enumerate(cobra_model.reactions):
        if x.id in quadaratic_dict:
           coef=quadaratic_dict[x.id]
        else:
            coef=0
        lp.objective.set_quadratic_coefficients(n, n, coef)
    """if quadratic_component is not None:
        set_quadratic_objective(lp, quadratic_component)"""
    #return   the_linear_expressions ,  constraint_limits, constraint_sense, constraint_names
    if relax_b:
        lp.linear_constraints.add(lin_expr=the_linear_expressions,rhs=constraint_limits, range_values=list(range_values), senses=constraint_sense, names=constraint_names)
                                  
    else:
        lp.linear_constraints.add(lin_expr=the_linear_expressions,rhs=constraint_limits, senses=constraint_sense, names=constraint_names)
    
    # Set the problem type as cplex doesn't appear to do this correctly
    lp.set_problem_type(Cplex.problem_type.QP)
    lp.objective.set_sense(1) #Minimize
    #lp.parameters.emphasis.numerical.set(1) #default for quadratic
    set_parameter(lp,"verbose",5)
    if out_name not in ("",None) and False:
       pass 
       #lp.write(out_name)
    return(lp)


def process_mta_results(analyzed_model,vref_dict={},vres_dict={},gene_target_dict={},up_genes=[],down_genes=[],p_value_dict={},key="MTA_results",omit_reactions=False,reactions_2_omit=[],reaction_pathway_dict={},use_only_first_pathway=True,signficant_met_dict={},signficant_kpc_dict={},signficant_genes_only=False):
  if reaction_pathway_dict in ("",{},None):
     omit_pathways=True
     reaction_pathway_dict={} #So that it remains empty in the next iteration
  else:
     omit_pathways=False 
  reactions_without_pathway=[x.id for x in analyzed_model.reactions if (reaction_pathway_dict.get(x.id)==[] or reaction_pathway_dict.get(x.id)==None)]
  for x in reactions_without_pathway:
    if "SINK" in x:
       reaction_pathway_dict[x]=["metabolomics"] 
    elif "EX_" in x:
       reaction_pathway_dict[x]=["Exchange_reactions"] 
    else:
       reaction_pathway_dict[x]=["no pathway assigned"]
  #print reactions_without_pathway
  variation_dict_dict={}
  reaction_dict_dict={}
  """reactions_2_omit=reactions2omit 
  omit_reactions=False
  analyzed_model=base_model
  use_only_first_pathway=False
  key="prova"
  """
  if True: ####This creates the dict
    variation_dict={}
    reaction_dict={}
    reaction_dict_dict[key]=reaction_dict
    variation_dict_dict[key]=variation_dict
    reaction_gene_dict={}
    """if "reverse" in key:
                factor=-1
    else:
                factor=1"""
    #factor =1 #This factor is to compensate if rxn_fbs has not been corrected before
    #print key, factor
    limited_output_flag=max([len(x.genes) for x in analyzed_model.reactions])==1 #In case we are working with reactions
    for reaction in analyzed_model.reactions:
          ###For display purposes if it is the reverse variation (which is usually the variation we are after) flip upregulated/downeragulated
          #break
          omited_flag=False
          reaction_id=reaction.id
          reaction_gene_dict[reaction_id]={"genes_up":[],"genes_down":[]}
          reaction_genes_up=[]#[str(x.id)+"("+gene_target_dict[x.id][reaction.id]+"/"+str(round(p_value_dict[x.id],3))+")" for x in reaction.genes if x.id in up_genes]
          for x in reaction.genes:
                if x.id in up_genes and not omited_flag:
                   try: 
                      if limited_output_flag:
                         reaction_str= str(gene_target_dict[x.id][reaction.id])
                      else:    
                          reaction_str=str(x.id)+"("+str(gene_target_dict[x.id][reaction.id])+"/"+str(round(p_value_dict[x.id],3))+")"
                      reaction_genes_up.append(reaction_str)
                      reaction_gene_dict[reaction_id]["genes_up"].append(x.id)
                   except:
                      pass     
          reaction_genes_down=[]
          for x in reaction.genes:
                if x.id in down_genes and not omited_flag:
                   try:
                     if limited_output_flag:
                         reaction_str= str(gene_target_dict[x.id][reaction.id])             
                     else:
                         reaction_str=str(x.id)+"("+str(gene_target_dict[x.id][reaction.id])+"/"+str(round(p_value_dict[x.id],3))+")"
                     reaction_genes_down.append(reaction_str)
                     reaction_gene_dict[reaction_id]["genes_down"].append(x.id)
                   except:
                     pass                     
          signficant_genes_flag=len(reaction_genes_up+reaction_genes_down)>0 or not signficant_genes_only
          for n,system in enumerate(reaction_pathway_dict[reaction_id]):
            if use_only_first_pathway and n>0:
                  break
            if system not in variation_dict:
               variation_dict[system]={"n_increase":0, "n_decrease":0,"increased":[],"decreased":[],"total_increase":0,"total_increase_genes":0,"total_decrease_genes":0,"total_decrease":0,"omitted":[],"n_omitted":0,"total_increase_normalized":0,"total_decrease_normalized":0,"normalized_n_reaction_score_increased":0,"normalized_n_reaction_score_decreased":0,"total_flux_vref":0,"total_flux_vres":0,"log2fc":0,"reaction_fc":[],"reaction_log2fc":[],"genes_up":[],"genes_down":[]} 
            vref=vref_dict[reaction_id]
            vres=vres_dict[reaction_id]
            variation_dict[system]["genes_up"]+=reaction_gene_dict[reaction_id]["genes_up"]
            variation_dict[system]["genes_down"]+=reaction_gene_dict[reaction_id]["genes_down"]
            if signficant_genes_flag and not (reaction_id in reactions_2_omit and omit_reactions):
               variation_dict[system]["total_flux_vref"]+=abs(vref)
               variation_dict[system]["total_flux_vres"]+=abs(vres)
            #variation_dict[system]["log2fc"]+=abs(vres)
            """if vref==0:
               continue""" 
            if reaction_id in reactions_2_omit and omit_reactions:
               #print "omited " +reaction_id
               variation_dict[system]["n_omitted"]+=1
               variation_dict[system]["omitted"].append(reaction_id)
               omited_flag=True
               diff="omited"
               #continue
            else:
               diff=abs(vres-vref)
               omited_flag=False 
               #continue            
            #reaction_genes_down=[str(x.id)+"("+str(round(p_value_dict[x.id],3))+")" for x in reaction.genes if x.id in down_genes]
            if (abs(vres)-abs(vref))>=1e-6 and not omited_flag:
               variation_dict[system]["n_increase"]+=1
               #variation_dict[system]["normalized_n_reaction_score_increased"]+=normalized_n_reaction_score
               variation_dict[system]["increased"].append(reaction_id)
               if signficant_genes_flag:
                  variation_dict[system]["total_increase"]+= (abs(vres)-abs(vref))
               #variation_dict[system]["total_increase_normalized"]+= abs(normalized_diff)
               """if success=="consistent":
                   variation_dict[system]["total_increase_genes"]+=1"""
            elif abs(vres)-abs(vref)<=-1e-6 and not omited_flag:
               variation_dict[system]["n_decrease"]+=1
               #variation_dict[system]["normalized_n_reaction_score_decreased"]+=normalized_n_reaction_score
               variation_dict[system]["decreased"].append(reaction_id)
               if signficant_genes_flag:
                  variation_dict[system]["total_decrease"]+= abs(abs(vres)-abs(vref))
               #variation_dict[system]["total_decrease_normalized"]+= abs(normalized_diff)
               """if success=="consistent":
                   variation_dict[system]["total_decrease_genes"]+=1"""
            if abs(vref)>1e-5:
               reaction_fc=(max(abs(vres),1e-5)/abs(vref))
               #print vres, vref, reaction_fc
               reaction_log2fc=math.log(reaction_fc,2)
               if signficant_genes_flag:
                  variation_dict[system]["reaction_fc"].append(reaction_fc)
                  variation_dict[system]["reaction_log2fc"].append(reaction_log2fc)
            else:
               reaction_log2fc=""
            try:
               reaction_dict[reaction_id]=[vref,vres,diff,reaction_log2fc,float((reaction_genes_up+reaction_genes_down)[0])] 
            except:
               reaction_dict[reaction_id]=[vref,vres,diff,reaction_log2fc,",".join(reaction_genes_up),",".join(reaction_genes_down)]
  ###Metabolomics: 
  for met in signficant_met_dict:
      reaction=signficant_met_dict[met]["reaction"]
      p_value=str(round(signficant_met_dict[met]['p_adj'],6))
      lfc=str(signficant_met_dict[met]['logfc'])
      #print signficant_met_dict
      #print met,signficant_met_dict[met]
      target_vres=str(round(signficant_met_dict[met]["target_vref"],7))
      metabolomics_str="logfc:"+lfc+" p_adj:"+str(p_value)+" target_vres:"+target_vres
      reaction_dict[reaction][-2]=metabolomics_str
  for met in signficant_kpc_dict:
      reaction=signficant_kpc_dict[met]["reaction"]
      p_value=str(round(signficant_kpc_dict[met]['p_adj'],6))
      lfc=str(signficant_kpc_dict[met]['logfc'])
      #print signficant_kpc_dict
      print met
      target_vres=str(round(signficant_kpc_dict[met]["target_vref"],7))
      metabolomics_str="logfc:"+lfc+" p_adj:"+str(p_value)+" target_vres:"+target_vres
      reaction_dict[reaction][-2]=metabolomics_str      
  output_sheet={}
  model=analyzed_model
  header_2=["rid","name","reaction","genes","vref","vres","variation","log2fc","reaction_genes_up","reaction_genes_down"]
  header1=["system","n_omited","omited","n_genes_up","n_genes_down","total_flux_vref","total_flux_vres","n_increase","n_decrease","increased","decreased","total_increase","total_decrease","total_variation","total_fold_change","Log2_fold_change","Average log2FC"]
  for key in variation_dict_dict:
    variation_dict= variation_dict_dict[key]
    reaction_dict=reaction_dict_dict[key]
    output_sheet[key+"_a"]=[header1]
    output_sheet[key+"_b"]=[]
    
    for system in variation_dict:
      if system=="no pathway assigned":# or (variation_dict[system]["total_increase"]-variation_dict[system]["total_decrease"])==0:
         continue
      fc=variation_dict[system]["total_flux_vres"]/max(variation_dict[system]["total_flux_vref"],1e-6)
      if fc!=0:
         logfc=math.log(fc,2)
      else:
         logfc=0
      if len(variation_dict[system]["reaction_log2fc"])>0:
         mean_log2=np.mean(variation_dict[system]["reaction_log2fc"])
      else:
         mean_log2=0    
      row=[system,variation_dict[system]["n_omitted"],str(variation_dict[system]["omitted"]),len(set(variation_dict[system]["genes_up"])),len(set(variation_dict[system]["genes_down"])),variation_dict[system]["total_flux_vref"],variation_dict[system]["total_flux_vres"],variation_dict[system]["n_increase"],variation_dict[system]["n_decrease"],str(variation_dict[system]["increased"]),str(variation_dict[system]["decreased"]),variation_dict[system]["total_increase"],variation_dict[system]["total_decrease"],variation_dict[system]["total_increase"]-variation_dict[system]["total_decrease"],fc,logfc,mean_log2]
      output_sheet[key+"_a"].append(row)
      output_sheet[key+"_b"].append(header1)
      output_sheet[key+"_b"].append(row)
      output_sheet[key+"_b"].append(["","increased"]+header_2)
      for reaction_id in variation_dict[system]["increased"]:
        """if reaction_id in reaction_signficance_dict:
           lfc=reaction_signficance_dict[reaction_id][n_lfc]
           p_value=reaction_signficance_dict[reaction_id][n_pvalue]
        else:
           lfc=""
           p_value=""
        """
        reaction=model.reactions.get_by_id(reaction_id)
        reaction_row=["","",reaction.id,reaction.name,reaction.reaction,reaction.gene_reaction_rule]+reaction_dict[reaction_id]
        output_sheet[key+"_b"].append(reaction_row)
      output_sheet[key+"_b"].append(["","decreased"]+header_2)
      for reaction_id in variation_dict[system]["decreased"]:
          """
          if reaction_id in reaction_signficance_dict:
             lfc=reaction_signficance_dict[reaction_id][n_lfc]
             p_value=reaction_signficance_dict[reaction_id][n_pvalue]
          else:
             lfc=""
             p_value="" 
          """
          reaction=model.reactions.get_by_id(reaction_id)
          reaction_row=["","",reaction.id,reaction.name,reaction.reaction,reaction.gene_reaction_rule]+reaction_dict[reaction_id]
          output_sheet[key+"_b"].append(reaction_row)
    
    #reaction_row=[reaction.id,reaction.name,reaction.reaction,reaction.gene_reaction_rule]+reaction_dict[reaction_id]
    #output_sheet[key+"all_reactions"].append(reaction_row)
    if omit_pathways:
        output_sheet={} 
    output_sheet[key+"_all_reactions"]=[["Reaction ID","Reaction Name","Reaction id","Reaction","Vref","Vres","Variation","Log2FC","Target flux"]]
    for reaction in analyzed_model.reactions:
        reaction_row=[reaction.id,reaction.name,reaction.reaction,reaction.gene_reaction_rule]
        if reaction.id in reaction_dict:
           reaction_row+=reaction_dict[reaction.id]  
        output_sheet[key+"_all_reactions"].append(reaction_row)
    
    return reaction_dict_dict, variation_dict_dict , output_sheet



def run_qMTA(target_model,reference_model,gene_fname,vref_dict={},gene_parameters={},gene_weight=0.5,unchanged_reaction_weight=0.5,reaction_pathway_dict={},key="",coef_precision=7,max_reactions4gene=10,use_only_first_pathway=False,output_omit_reactions_with_more_than_max_genes=False,output_signficant_genes_only=False,normalize_by_scale_genes=True,min_flux4weight=1e-6,normalize_by_scale_unchanged_reactions=True,differential_expression_sheet_dict={},min_flux_fold_change=1e-9,qpmethod=1,n_threads=0,sample_name="",debug_prefix="",non_gene_met_reaction_weight_dict={},detailed_output=True):
    ref_gene_parameters={"log2_str":"log2FoldChange","log2_factor":1,"padj_str":"padj","p_th":0.25,"log2fc_th":0,"gene_str":"NCBI.gene.ID","p_weight_formula":"-1*math.log(p_value,10)","ignore_p_value":False}
    ref_gene_parameters.update(gene_parameters)
    log2_str=ref_gene_parameters["log2_str"]
    log2_factor=ref_gene_parameters["log2_factor"]
    padj_str=ref_gene_parameters["padj_str"]
    p_th=ref_gene_parameters["p_th"]
    log2fc_th=ref_gene_parameters["log2fc_th"]
    gene_str=ref_gene_parameters["gene_str"]
    p_weight_formula=ref_gene_parameters["p_weight_formula"]
    ignore_p_value=ref_gene_parameters["ignore_p_value"]
     
    up_genes, down_genes, log2fold_change_dict,   p_value_dict ,  gene_weight_dict, gene_weight_normalized_dict,=read_gene_data(fname=gene_fname,model=reference_model,log2_str=log2_str,log2_factor=log2_factor,padj_str=padj_str,p_th=p_th,log2fc_th=log2fc_th,gene_str=gene_str,p_weight_formula=p_weight_formula,sheet_dict=differential_expression_sheet_dict,ignore_p_value=ignore_p_value)
    
    quadaratic_dict, coefficient_dict, gene_target_dict, signficant_gene_list_corrected,reactions2omit=build_weight_dicts(up_genes+down_genes,reference_model,max_reactions4gene=max_reactions4gene,gene_weight=gene_weight,non_gene_met_reaction_weight=unchanged_reaction_weight,gene_weight_dict=gene_weight_normalized_dict,non_gene_met_reaction_weight_dict=non_gene_met_reaction_weight_dict,fold_change_dict=log2fold_change_dict,vref_dict=vref_dict,max_fold_change=99999999999,normalize_by_scale_genes=normalize_by_scale_genes,normalize_by_scale_unchanged_reactions=normalize_by_scale_unchanged_reactions,min_flux4weight=min_flux4weight,genes_log2_to_lineal=True,precision=coef_precision,min_flux_fold_change=min_flux_fold_change)
    problem=create_rMTA_problem_quadratic(target_model,quadaratic_dict=quadaratic_dict, coefficient_dict=coefficient_dict,out_name="rMTA.lp",qpmethod=qpmethod)
    problem.parameters.timelimit.set(300)
    problem.parameters.threads.set(n_threads)
    problem.parameters.emphasis.numerical.set(1)
    #problem.parameters.emphasis.memory.set(True) This will reduce memory usage but also lower performance
    #problem.parameters.workmem.set(3038) #after this value it starts to compress data
    
    try:
      problem.solve()
      status = problem.solution.get_status_string().lower()
      if status!="optimal":
          #raw_input("Press Enter to continue...") 
          problem.parameters.qpmethod.set(1)
          try:
             print "Trying qpmethod 1"
             problem.solve()
             status = problem.solution.get_status_string().lower()
             if status!="optimal":
                raise Exception('solver error')
          except:
             #Reset parameters and test all methods
             problem.parameters.reset()
             problem.parameters.timelimit.set(150)
             problem.parameters.threads.set(n_threads)
             for emphasis_numerical in (0,1):
               for method in  (2,3,4,5,6,4,1):
                 try:
                     print "Trying qpmethod",method,"emphasis",emphasis_numerical
                     problem.parameters.emphasis.numerical.set(emphasis_numerical)
                     problem.parameters.qpmethod.set(method)
                     problem.solve()
                     status = problem.solution.get_status_string().lower()
                     print "qpmethod",method,"emphasis",emphasis_numerical,status
                     if status=="optimal": break
                 except:
                     print "solver error"
               if status=="optimal": break
      if  status!="optimal": #Tunning parameters should not do anything not already done, but we leave it just in case
          #raw_input("Press Enter to continue...")   
          problem.parameters.qpmethod.set(1) 
          print "tuning parameters"
          try:
             problem.parameters.tune_problem()
             problem.parameters.write_file("parameters.txt")
             problem.solve()
             status = problem.solution.get_status_string().lower()
             print "tunning",status
             if status!="optimal":
                raise Exception('solver error')
          except:
             print "Swicthing to default qpMethod to attempt to get a solution even if its not fully optimal" 
             problem.parameters.qpmethod.set(0)
             problem.solve()
             """if status!="optimal":
                problem.write("debug.lp") 
                f=open(debug_prefix+"debug.txt","a")
                f.write(sample_name+"\n")
                f.close()"""
                
      solution=cplex_format_solution(problem, target_model) 
      for reaction in reference_model.reactions:
        if reaction.id not in solution.x_dict:
           solution.x_dict[reaction.id]=0 
    except: raise Exception('solver error')
    
    if detailed_output:       
       reaction_dict_dict, variation_dict_dict , output_sheet=process_mta_results(reference_model,vref_dict=vref_dict,vres_dict=solution.x_dict,gene_target_dict=gene_target_dict,up_genes=up_genes,down_genes=down_genes,p_value_dict=p_value_dict,key=key,omit_reactions=output_omit_reactions_with_more_than_max_genes,reactions_2_omit=reactions2omit,reaction_pathway_dict=reaction_pathway_dict,use_only_first_pathway=use_only_first_pathway,signficant_genes_only=output_signficant_genes_only)
    else:
       reaction_dict_dict={}
       variation_dict_dict={}
       output_sheet={}
    return output_sheet, solution.x_dict, reaction_dict_dict, variation_dict_dict,up_genes, down_genes, log2fold_change_dict,   p_value_dict ,  gene_weight_dict








