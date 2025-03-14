import re
import copy
#import Tkinter
#import tkFileDialog
import csv
from openpyxl import load_workbook
import math

import math

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


def read_flux_constraints(metabolic_model,file_name=None,create_copies=True,delimiter=",",boundary_precision=1e-6):
    precision=int(-1*(math.log10(boundary_precision)))
    """if file_name==None:
      tk=Tkinter.Tk()
      tk.withdraw()
      loaded_file = tkFileDialog.askopenfile(title='Choose a constraint file',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')]) 
      file_name=loaded_file.name
      tk.destroy()
      if len(file_name)==0:
         return model,ratio_dict"""
    if create_copies:
       model=copy.deepcopy(metabolic_model)
       ratio_dict={} 
    else:
       model=metabolic_model
    ratio_re=re.compile("^(.+)/(.+)")
    row_list=[]
    if ".xlsx"in file_name:
        wb = load_workbook(file_name, read_only=True,data_only=True)
        ws=wb.active
        for n,row_xlsx in enumerate(ws.rows):
            row=[]
            for cell in row_xlsx:   
                row.append(cell.value)
            row_list.append(row)
    else:
        csv_file=open(file_name)
        csv_reader=csv.reader(csv_file,delimiter=delimiter)
        row_list=list(csv_reader)
        csv_file.close()
        print row_list  
    for n,row in enumerate(row_list):
        row_len=len(row)
        if row_len==0:
           continue
        if row[0]!=None and row[0]!="":
           if "#" in row[0]:
               continue
           contraint_id=row[0].replace(" ","") #Remove empty spaces
           if contraint_id in model.reactions:
              reaction_id=str(contraint_id)
              reaction=model.reactions.get_by_id(reaction_id)
              if row_len>1:
                 if row[1]!=None and row[1]!="":
                    try:
                      reaction.lower_bound=round(float(row[1]),precision)
                      print reaction_id+": lower bound set to "+ str(reaction.lower_bound)
                    except:
                       print reaction_id+": lower bound not valid"
              if row_len>2:
                 if row[2]!=None and row[2]!="":
                    try:
                      reaction.upper_bound=round(float(row[2]),precision)
                      print reaction_id+": upper bound set to "+ str(reaction.upper_bound)
                    except:
                       print reaction_id+": upper bound not valid"
              if row_len>3:
                 if row[3]!=None and row[3]!="":
                    try:
                      reaction.objective_coefficient=round(float(row[3]),precision)  
                      print reaction_id+": objective coefficient set to "+ str(reaction.objective_coefficient)
                    except:
                      print reaction_id+": objective coefficient not valid"
                 #print model.optimize()
    
    return model   

