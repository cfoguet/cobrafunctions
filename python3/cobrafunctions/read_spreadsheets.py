import csv
from openpyxl import load_workbook
import re
#import Tkinter
#import tkFileDialog
import copy


def read_spreadsheets(file_names=None,csv_delimiter=',',more_than_1=True,tkinter_title="Chose a file"):
    """if file_names==None:
      tk=Tkinter.Tk()
      tk.withdraw()
      if more_than_1==True:
         loaded_files = tkFileDialog.askopenfiles(title=tkinter_title,filetypes=[("csv","*.csv"),("xlsx","*.xlsx"),('All files','*.*')]) 
         file_names=[x.name for x in loaded_files]
      else:
         loaded_file = tkFileDialog.askopenfile(title=tkinter_title,filetypes=[("csv","*.csv"),("xlsx","*.xlsx"),('All files','*.*')])  
         file_names=[loaded_file.name]
      tk.destroy()
      print file_names"""  
    if not isinstance(file_names,list):
       file_names=[file_names]
    condition_rows_dict={}
    for file_name in file_names:
        file_type=file_name[-6:] #Take some extra caracters just in case
        if any(x in file_type for x in [".xlsx",".xlsm",".xltx",".xltm"]):
           print(file_type)
           wb = load_workbook(file_name, read_only=True,data_only=True)
           for ws in wb.worksheets:
                  condition=ws.title
                  #
                  if "Sheet" in condition:
                     condition=condition.replace("Sheet","Label_")
                  elif "Hoja" in condition:
                     condition=condition.replace("Hoja","Label_") 
                  condition_rows_dict[condition]=[]
                  for xlsx_row in ws.rows:
                      row=[]
                      for cell in xlsx_row:
                          row.append(cell.value)
                      condition_rows_dict[condition].append(row)
        else:
           csv_file=open(file_name)
           csv_reader=csv.reader(csv_file,delimiter=csv_delimiter)
           row_list=list(csv_reader)
           condition=(file_name.split(".")[0]).split("/")[-1] #Get the name of the file minus the extension and the path
           condition_rows_dict[condition]=row_list
           csv_file.close() 
    return condition_rows_dict



               
