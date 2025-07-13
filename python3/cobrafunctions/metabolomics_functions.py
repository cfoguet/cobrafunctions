# -*- coding: utf-8 -*-
import cobra
import math
import numpy as np
import sys
import operator
import pandas as pd

import copy
from warnings import warn

import scipy

from cobra import Reaction, Metabolite
from .write_spreadsheet import write_spreadsheet
from .read_spreadsheets import read_spreadsheets
from .cobra_functions import *
from . import cobra_functions


#######Metabolomics functions

def add_sink_reactions_with_multiple_compartments(model,stat_dict,metabolite_name_compartment_dict={},lb=None,ub=None,condition="Control",precision=7,factor=1,compartments=[]):
    #sink_met_dict={}
    met_sink_dict= {}
    rejected_list=[]
    for met in stat_dict:
      if met in metabolite_name_compartment_dict:
        metabolite_id=metabolite_name_compartment_dict[met]['met_id'] #Metabolite ID without compartment tag
        rid='SINK_'+metabolite_id
        met_sink_dict[met]=rid   
        if rid not in model.reactions:
          reporter_metabolite=Metabolite("reporter_"+met)
          for compartment in metabolite_name_compartment_dict[met]["compartment"]:
               met_id=metabolite_id+"_"+compartment
               if met_id in model.metabolites:
                 sink_rid='SINK_'+met_id
                 print(met_id, compartment, metabolite_id, rid, met)
                 sink_reaction=Reaction(sink_rid)
                 sink_reaction.name = 'sink of'+met_id
                 sink_reaction.subsystem = 'sink'
                 sink_reaction.add_metabolites({reporter_metabolite: 1.0,model.metabolites.get_by_id(met_id):-1})
                 model.add_reactions([sink_reaction])
          #Full reaction
          reaction = Reaction(rid)
          #met_sink_dict[rid]=met
          reaction.name = 'sink of'+metabolite_id
          reaction.subsystem = 'sink'
          reaction.add_metabolites({reporter_metabolite: -1.0})   
          model.add_reactions([reaction])
        else:
           reaction=model.reactions.get_by_id(rid)
        if lb!=None:
          reaction.lower_bound=round(stat_dict[met][condition][lb]*factor,precision)
        if ub!=None: 
          reaction.upper_bound=round(stat_dict[met][condition][ub]*factor,precision)
          
      else:
       rejected_list.append(met)
    print("done")   
    return met_sink_dict, rejected_list

#formerly read_biocrates_metabolomics_data
def read_metabolomics_data(metabolomics_file,max_NA_fraction=0.5,lower_metabolite_names=True):
    #Metaboanlyst like, Assume samples in rows and the fisrt row is samples name followes by condition
    data_sheet=read_spreadsheets(metabolomics_file)
    data_dict={}
    n_met_dict={}
    for sheet in data_sheet:
      na_counter={}  
      max_na=max_NA_fraction*(len(data_sheet[sheet])-1) 
      for n_row,row in enumerate(data_sheet[sheet]):
       print(row)
       condition=row[1]
       sample=row[0]
       for n_col,col in enumerate(row):
           if n_col>1:
              if n_row==0:
                 if lower_metabolite_names:
                    col=col.lower() 
                 n_met_dict[n_col]=col
                 print(n_met_dict)
              else:
                 if condition=="":
                    continue
                 met=n_met_dict[n_col] 
                 #try: # to check if float
                 try:
                   value=float(col)
                 except:
                     #"""
                     #print na_counter
                     if met not in na_counter:
                        na_counter[met]=1
                     else:
                        na_counter[met]+=1
                     #"""  
                     #continue
                     value="Na"
                 if value=="Na":
                    continue 
                 if met not in data_dict:
                    data_dict[met]={}
                 #print data_dict
                 if condition not in data_dict[met]:
                    data_dict[met][condition]=[]
                 data_dict[met][condition].append(value)
                 #print data_dict
      #Remove NA
      for met in na_counter:
        n_na=na_counter[met]
        if n_na>=max_na and met in data_dict:
           del(data_dict[met]) 
    stat_dict={}
    for met in data_dict:
           if lower_metabolite_names:
              met=met.lower()
           stat_dict[met]={}
           for condition in data_dict[met]:
               stat_dict[met][condition]={}
               values=data_dict[met][condition]
               mean=np.mean(values)
               std=np.std(values)
               stat_dict[met][condition]["mean"]=mean
               stat_dict[met][condition]["std"]=std
               stat_dict[met][condition]["se"]=std/float(len(values))
               stat_dict[met][condition]["values"]=values
    return stat_dict


def statistical_difference_metabolomics(stat_dict,cond1="target",cond2="control",convert_to_log=True, p_adj_th=0.05,met_list=None,p_weight_formula="-1*math.log(p_value,10)",normalize_p_weight=True,met_sink_dict={},log2_factor_met=1):
    #log2_factor_met is used to do the reverse MTA
    if met_list==None:
       met_list=list(stat_dict.keys()) 
    p_value_list=[]
    signficant_met_dict={}
    for met in sorted(met_list):
        if met not in stat_dict:
           continue
        if cond1 not in stat_dict[met]:
           continue 
        if cond2 not in stat_dict[met]:
           continue 
        values1=stat_dict[met][cond1]["values"]
        values2=stat_dict[met][cond2]["values"]
        mean1=stat_dict[met][cond1]["mean"]
        mean2=stat_dict[met][cond2]["mean"] 
        if convert_to_log:
           comparisson=cond1+"/"+cond2 
           values1=[math.log(x,2) for x in values1] 
           values2=[math.log(x,2) for x in values2]
        else:
           comparisson=cond1+"-"+cond2  
        result=scipy.stats.ttest_ind(values1, values2, axis=None, equal_var=True, nan_policy='propagate',alternative='two-sided')
        p_value=result[1]
        p_value_list.append(p_value)
        fc=mean1/mean2
        if convert_to_log:
           lfc= math.log(fc,2)*log2_factor_met
        else:
            lfc=None
        if log2_factor_met==1:
           local_dict={"p_value":p_value,"fc":fc,"logfc":lfc,"mean_source":mean2,"mean_target":mean1}
        elif log2_factor_met==-1:
           local_dict={"p_value":p_value,"fc":fc,"logfc":lfc,"mean_source":mean1,"mean_target":mean2}
        else:
            raise Exception('log2_factor_met should be 1 or -1')
        if "statistics" not in stat_dict:
            stat_dict[met]["statistics"]={}
        stat_dict[met]["statistics"][comparisson]=local_dict
    ##Correct p value
    corrected_p_values=scipy.stats.false_discovery_control(p_value_list, axis=0, method='bh') #list(r_stats.p_adjust(FloatVector(p_value_list), method = 'BH'))
    n_counter=0
    for met in sorted(met_list):
         if met not in stat_dict:
            continue
         if "statistics" not in stat_dict[met]:
             continue 
         stat_dict[met]["statistics"][comparisson]["p_adj"]=corrected_p_values[n_counter]
         p_value=corrected_p_values[n_counter]
         values1=stat_dict[met][cond1]["values"] #To be used in formulas
         values2=stat_dict[met][cond2]["values"] #To be used in formulas        
         sd1=np.std(values1) #To be used in formulas
         sd2=np.std(values2) #To be used in formulas
         stat_dict[met]["statistics"][comparisson]["p_weight"]=eval(p_weight_formula)
         if  corrected_p_values[n_counter]<p_adj_th:
             signficant_met_dict[met]=stat_dict[met]["statistics"][comparisson]
             signficant_met_dict[met]["reaction"]=met_sink_dict.get(met)
         n_counter+=1
    
    if normalize_p_weight:
       p_factor=p_mean=np.mean([ signficant_met_dict[met]["p_weight"] for met in signficant_met_dict])
    else:
       p_factor=p_mean=1.0 
    p_mean_normalized_dict={met:round_sig(signficant_met_dict[met]["p_weight"]/p_mean,3) for met in signficant_met_dict}
    for met in signficant_met_dict:
        signficant_met_dict[met]["p_weight_norm"]=p_mean_normalized_dict[met]
    return signficant_met_dict


def modify_model_for_seahorse_recon3d(model,remove=False):
   #Shared Token metabolites
   print("CYOOm_token_atpsyn" not in model.metabolites and "ATPS4mi_token" not in model.metabolites and "NADH2_u10mi_token" not in model.metabolites)
   if "CYOOm_token_atpsyn" not in model.metabolites and "ATPS4mi_token" not in model.metabolites and "NADH2_u10mi_token" not in model.metabolites: #Other checks could be added but in theory those should be enought
      CYOOm_token_atpsyn=Metabolite("CYOOm_token_atpsyn"); CYOOm_token_atpsyn.name="CYOOm_token_atpsyn"
      ATPS4mi_token=Metabolite("ATPS4mi_token"); ATPS4mi_token.name="ATPS4mi_token"
      NADH2_u10mi_token=Metabolite("NADH2_u10mi_token"); NADH2_u10mi_token.name="NADH2_u10mi_token"
      fadhdh_like_token=Metabolite("fadhdh_like_token"); fadhdh_like_token.name="fadhdh_like_token"
      sulfox_like_token=Metabolite("sulfox_like_token"); sulfox_like_token.name="sulfox_like_token"
      
      #Nadhdh
      model.reactions.NADH2_u10mi.add_metabolites({NADH2_u10mi_token:1}) #	5.0 h_m + nadh_m + q10_m --> 4.0 h_i + nad_m + q10h2_m + NADH2_u10mi_token
      RGROUP_non_atpsyn_NADH2_u10mi=reaction_from_string(model,"RGROUP_non_atpsyn_NADH2_u10mi","NADH2_u10mi_token ->",bounds=None,gpr="")
      RGROUP_atpsyn_NADH2_u10mi=reaction_from_string(model,"RGROUP_atpsyn_NADH2_u10mi","NADH2_u10mi_token -> 0.5 CYOOm_token_atpsyn + 2.5 ATPS4mi_token ",bounds=None,gpr="")
      
      FADH_Like_list=[]
      for reaction in model.metabolites.get_by_id("q10_m").reactions:
          if reaction.metabolites[model.metabolites.get_by_id("q10_m")]==-1:
             if reaction.metabolites.get(model.metabolites.get_by_id("q10h2_m"))==1:
                 if model.metabolites.get_by_id("h_i") not in reaction.metabolites:
                    print(reaction.id, reaction.reaction) 
                    FADH_Like_list.append(reaction)
      
      for reaction in FADH_Like_list:
          reaction.add_metabolites({fadhdh_like_token:1}) #  etfrd_m + q10_m --> etfox_m + q10h2_m + ETFQO_token
      
      RGROUP_non_atpsyn_fadhdh_like=reaction_from_string(model,"RGROUP_non_atpsyn_fadhdh_like","fadhdh_like_token ->",bounds=None,gpr="") 
      RGROUP_atpsyn_fadhdh_like=reaction_from_string(model,"RGROUP_atpsyn_fadhdh_like","fadhdh_like_token -> 0.5 CYOOm_token_atpsyn + 1.5 ATPS4mi_token ",bounds=None,gpr="")
      ###sulfox Like
      SULFOX_Like_list=[]
      for reaction in model.metabolites.get_by_id("ficytC_m").reactions:
          if reaction.metabolites[model.metabolites.get_by_id("ficytC_m")]==-2:
             if reaction.metabolites[model.metabolites.get_by_id("focytC_m")]==2:
                 if model.metabolites.get_by_id("h_i") not in reaction.metabolites:
                    print(reaction.id, reaction.reaction) 
                    SULFOX_Like_list.append(reaction)
      
      for reaction in SULFOX_Like_list:
          reaction.add_metabolites({sulfox_like_token:1}) #  etfrd_m + q10_m --> etfox_m + q10h2_m + ETFQO_token
      
      RGROUP_non_atpsyn_sulfox_like=reaction_from_string(model,"RGROUP_non_atpsyn_sulfox_like","sulfox_like_token ->",bounds=None,gpr="") 
      RGROUP_atpsyn_sulfox_like=reaction_from_string(model,"RGROUP_atpsyn_sulfox_like","sulfox_like_token -> 0.5 CYOOm_token_atpsyn + 0.5 ATPS4mi_token ",bounds=None,gpr="")
      #ATP synthase
      model.reactions.ATPS4mi.add_metabolites({ATPS4mi_token:-1})
      #O2 linked to ATPsynthase
      RGROUP_CYOOm_atpsyn=reaction_from_string(model,"RGROUP_CYOOm_atpsyn","CYOOm_token_atpsyn -> ",bounds=None,gpr="")
      RGROUP_CYOOm3i_CYOOm2i=define_reaction_group(model,{"CYOOm3i":1,"CYOOm2i":1},group_reaction_id="RGROUP_CYOOm3i_CYOOm2i",lower_bound=None,upper_bound=None,objective_coefficient=0)
      #identfy proton transport reactions: 
      proton_transport_dict={}
      for reaction in model.metabolites.get_by_id("h_e").reactions:
          if model.metabolites.get_by_id("h_c") in reaction.metabolites: #It is a transport reaction
             factor=reaction.metabolites[model.metabolites.get_by_id("h_e")]
             print(reaction.reaction, factor)
             proton_transport_dict[reaction.id]=factor
      
      RGROUP_h_transport=define_reaction_group(model,proton_transport_dict,group_reaction_id="RGROUP_h_transport",lower_bound=None,upper_bound=None,objective_coefficient=0)
   if remove:
       #Mets
       CYOOm_token_atpsyn=model.metabolites.CYOOm_token_atpsyn
       ATPS4mi_token=model.metabolites.ATPS4mi_token
       NADH2_u10mi_token=model.metabolites.NADH2_u10mi_token
       fadhdh_like_token=model.metabolites.fadhdh_like_token
       sulfox_like_token=model.metabolites.sulfox_like_token
       mRGROUP_CYOOm3i_CYOOm2i=model.metabolites.mRGROUP_CYOOm3i_CYOOm2i
       mRGROUP_h_transport=model.metabolites.mRGROUP_h_transport
       
       mets2remove=[CYOOm_token_atpsyn,ATPS4mi_token,NADH2_u10mi_token,fadhdh_like_token,sulfox_like_token,mRGROUP_CYOOm3i_CYOOm2i,mRGROUP_h_transport]
       print("mets")
       for met in mets2remove:
           met.remove_from_model(destructive=False)
       #Reactions
       reaction_id=["RGROUP_non_atpsyn_NADH2_u10mi","RGROUP_atpsyn_NADH2_u10mi","RGROUP_non_atpsyn_fadhdh_like","RGROUP_atpsyn_fadhdh_like","RGROUP_non_atpsyn_sulfox_like","RGROUP_atpsyn_sulfox_like","RGROUP_CYOOm_atpsyn","RGROUP_h_transport","RGROUP_CYOOm3i_CYOOm2i"]
       print("reactions")
       for rid in  reaction_id:
           reaction=model.reactions.get_by_id(rid)
           reaction.remove_from_model()  
