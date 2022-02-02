import string
import openpyxl 
import csv

alphabet=list(string.ascii_uppercase)
alphabet2=[]
for letter1 in alphabet:
    for letter2 in alphabet:
        alphabet2.append(letter1+letter2) 

alphabet+=alphabet2

def write_spreadsheet(file_name,sheet_row_data_dict,sheet_order=None,force_extension=False):
    if sheet_order==None or sheet_order==[]:
       sheet_order=sorted(sheet_row_data_dict)
    print(file_name)
    if ".xlsx" in file_name:
       wb = openpyxl.Workbook()
       for n_sheet,sheet_id in enumerate(sheet_order):
           sheet_title = (str(n_sheet)+"_"+sheet_id[:27] + '..') if len(sheet_id) > 31 else sheet_id
           if n_sheet==0:
              sheet = wb.active
              sheet.title=sheet_title
           else:
                          
              sheet=wb.create_sheet(title=sheet_title)
           for n_row, row in enumerate(sheet_row_data_dict[sheet_id]):
               for n_col,data in enumerate(sheet_row_data_dict[sheet_id][n_row]):
                    try:
                      sheet[alphabet[n_col]+str(n_row+1)]=data
                    except:
                      sheet[alphabet[n_col]+str(n_row+1)]="encoding_error" 
       wb.save(file_name)
       return
    else:
        if ".txt" not in file_name.lower() and ".csv" not in file_name.lower() and force_extension: 
           file_name+=".csv" 
        csvFile = open(file_name, 'w')
        outputWriter = csv.writer(csvFile)
        for n_sheet,sheet_id in enumerate(sheet_order):
            if len(sheet_order)>1:
               outputWriter.writerow(['###'+sheet_id+"###"])
            for n_row, row in enumerate(sheet_row_data_dict[sheet_id]):
                outputWriter.writerow(row)
            if len(sheet_order)>1: 
                outputWriter.writerow(["######"])
        csvFile.close()
        return       
 
#write_spreadsheet("essborram.csv",sheet_row_data_dict={"hola":[[1,2],[3,4]],"adeu":[["a","b"],["c","d"]]},sheet_order=["hola","adeu"])

def cbind_spreadsheets(file_names=[],out_name="",csv_delimiter=',',more_than_1=True,tkinter_title="Chose a file",add_header=True,default_cols=50):
    ###Get dimensions
    sheet_rows_dict={}
    n_rows=0
    n_cols=0
    for file_name in sorted(file_names):
        condition_rows_dict=read_spreadsheets(file_names=file_name,csv_delimiter=csv_delimiter,more_than_1=more_than_1,tkinter_title="Chose a file")
        for sheet in condition_rows_dict:
            sheet_rows_dict[file_name+sheet]=condition_rows_dict[sheet]
            local_nrows=len(sheet)
            sheet_rows=condition_rows_dict[sheet]
            n_rows=max(len(sheet_rows),n_rows)
            n_cols+=max([len(x) for x in sheet_rows])
    rows=[]
    print(n_rows,n_cols)
    #print sheet_rows_dict
    for n in range(0,n_rows):
        row=[]
        for sheet in sorted(sheet_rows_dict):
            try:
                if n==0:
                   mod_row=[str(x)+"__"+sheet for x in sheet_rows_dict[sheet][n]]
                   row+=mod_row
                else:
                    row+=sheet_rows_dict[sheet][n]
                    
            except:
                row+=[""]*len(sheet_rows_dict[sheet][0])
            #print sheet_rows_dict[sheet][n]
            """if len(sheet_rows_dict[sheet])<n:
               row+=sheet_rows_dict[sheet][n]
            else:
               row+=[""]*len(sheet_rows_dict[sheet][0])
               print "case 2"""
        
        rows.append(row)
    #print row,n_rows,n_cols    
    #print rows
    write_spreadsheet(out_name,{"merged":rows},sheet_order=None,force_extension=False)
    
