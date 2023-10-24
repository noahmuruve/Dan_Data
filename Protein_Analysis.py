# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 08:01:52 2023

@author: noahm
"""

##############################################################################
directory='C:/Users/noahm/Desktop/projects/Dan_Data/'  
##############################################################################

import os
import pandas as pd
import math as m
import numpy as np
import scipy.stats as ss

class file_info:
    def __init__(self, file_name, protein_list):
        self.file_name = file_name
        self.protein_list = protein_list

class protein:
    def __init__ (self, protein_ID, name, symbol, uniprot, sample_list, calculated_values_list):
        self.protein_ID = protein_ID
        self.name=name
        self.symbol=symbol
        self.uniprot=uniprot
        self.sample_list = sample_list
        self.calculated_values_list = calculated_values_list

class sample:
    def __init__ (self, patient_ID, sex, age, treatment, time, matrix, value, fold_diff_time):
        self.patient_ID = patient_ID
        self.sex = sex
        self.age = age
        self.treatment = treatment
        self.time = time
        self.matrix = matrix
        self.value = value
        self.fold_diff_time = fold_diff_time
        
class calculated_values:
    def __init__ (self,
                sex = "",
                age = 0,
                treatment = "",
                time = "",
                matrix = "",
                
                values = [],
                fold_diff_time_list = [],
                
                mean = 0,
                std = 0,
                ttest = 0,
                mw = 0,
                
                fold_diff_treatment_mean = False,
                fold_diff_time_mean = False,
                fold_diff_time_std = False,
                fold_diff_time_ttest = False,
                fold_diff_time_mw = False,
                
                fold_diff_treatment_time_mean = False):

        self.sex = sex
        self.age = age
        self.treatment = treatment
        self.time = time
        self.matrix = matrix
        
        self.values = values
        self.fold_diff_time_list = fold_diff_time_list
        
        self.mean = mean
        self.std = std
        self.ttest = ttest
        self.mw = mw
        
        self.fold_diff_treatment_mean = fold_diff_treatment_mean
        self.fold_diff_time_mean = fold_diff_time_mean
        self.fold_diff_time_std = fold_diff_time_std
        self.fold_diff_time_ttest = fold_diff_time_ttest
        self.fold_diff_time_mw = fold_diff_time_mw
        
        self.fold_diff_treatment_time_mean = fold_diff_treatment_time_mean

class order_ob:
    def __init__(self, typ, value, options):
        self.typ = typ
        self.value = value
        self.options = options
        
# read each file in the directory #
def Read_Data(direc, datasub):
    error_log = []
    # --------------------------------------------------------------------
    ## read key ##
    key_file="key/AB002 Subjects - Demographics.csv"
    key=pd.read_csv(open(direc + key_file, 'r'),header=None)
    
    # --------------------------------------------------------------------
    ## generate protein_list for file ##
    protein_list = []
    first_file = True
    for file in os.listdir(direc + datasub):
        raw_file = pd.read_csv(open(direc + datasub + '/' + file, 'r'))
        
        # --------------------------------------------------------------------
        ## organize raw_file for easier column searching ##
        # find title row for single file #
        i=0
        while pd.isnull(raw_file.iat[i,0]):
            i=i+1
        sample_title_row=i
        
        # fix column titles #
        j=0
        for title in raw_file.columns:
            if title.startswith("Unnamed:"):
                swap = { title : raw_file.iat[sample_title_row,j] }
                raw_file.rename(columns=swap, inplace=True)
                
            j=j+1
        # --------------------------------------------------------------------
        ## organize raw_file for easier row searching ##
        # Find title column for single file
        j = 0
        while pd.isnull(raw_file.iat[0,j]):
            j=j+1
        sample_title_column=j
        
        # create title list #
        title_list = []
        i=0
        for index, row in raw_file.iterrows():
            if pd.isnull(raw_file.iat[i,sample_title_column]) and i != sample_title_row:
                title = row.SampleId
            elif i != sample_title_row:
                title = raw_file.iat[i,sample_title_column]
            title_list.append(title)
            i=i+1
            
        raw_file.index = title_list
        
        # --------------------------------------------------------------------
        # Populate Protein Info #
        for i in range(sample_title_column+1, len(raw_file.columns)):
            for j in range(sample_title_row+1, len(raw_file)):
                
                # if first time looking at protein generate general protein data #
                if j == sample_title_row+1 and first_file:
                    protein_list.append(protein(
                    raw_file.columns[i],
                    raw_file.loc['TargetFullName'][i],
                    raw_file.loc['EntrezGeneSymbol'][i],
                    raw_file.loc['UniProt'][i],
                    [],
                    [])
                    )
                
                    if not protein_list[-1].protein_ID or not protein_list[-1].name or not protein_list[-1].symbol or not protein_list[-1].uniprot:
                        error_msg = "Error, missing protein info for column " + str(i)
                        if error_msg not in error_log:
                            error_log.append(error_msg)
                
                # Get Sample Patient ID #
                patient_ID = raw_file.SampleId[j].split()[0].strip()
                if len(patient_ID) == 0:
                    error_msg = "Error, empty patient ID for: " + str(raw_file.SampleId[j].strip())
                    if error_msg not in error_log:
                        error_log.append(error_msg)
                
                # Get Sample Time #
                temp_time = raw_file.SampleId[j].split()[1:]
                timevalid=False
                time = "NA"
                for item in temp_time:
                    if item.startswith("D"):
                        day=item[1:]
                        timevalid=True
                        if float(day) == 1:
                            time = "Day_0" #zero
                        elif float(day) <= 3:
                            time = "Day_3" #3
                        elif float(day) > 3:
                            time = "EOT" #EOT
                        else:
                            error_msg = "Error, sample time not valid" + str(patient_ID) + str(temp_time)
                            if error_msg not in error_log:
                                error_log.append(error_msg)
                            
                    elif item == "EOT":
                        time = "EOT"
                        timevalid=True
                
                if not timevalid:
                    error_msg = "Error, sample time not valid for" + str(patient_ID) + str(temp_time)
                    if error_msg not in error_log:
                        error_log.append(error_msg)
                
                # look through key for patient info #
                subfound=False
                for sub in key.iterrows():
                    if patient_ID==sub[1][0]:
                        sex=sub[1][1].strip()
                        age=sub[1][2] #this one is int so does not need strip()
                        treatment=sub[1][3].strip()
                        subfound=True
                        break
                    
                # if not found means its a control or something # 
                if not subfound:
                    sex="NA"
                    age="NA"
                    treatment="NA"
                    error_msg = "sample ID not found on key:" + str(patient_ID)
                    if error_msg not in error_log:
                        error_log.append(error_msg)
                    
                if first_file:
                    protein_idx = -1
                else:
                    protein_idx = i-sample_title_column-1
                    assert raw_file.columns[i] == protein_list[protein_idx].protein_ID
                        
                # create protein with all the values from sample #    
                protein_list[protein_idx].sample_list.append(sample(
                           patient_ID,
                           sex,
                           age,
                           treatment,
                           time,
                           raw_file.SampleMatrix[j],
                           float(raw_file.iat[j,i]),
                           False))
        
        first_file = False
        
    # --------------------------------------------------------------------
    # now that we've fully populated all sample lists, we can reliably #
    # calculate fold_diffs #
    for prot in protein_list:
        for samp in prot.sample_list:
            
            # time fold diffs are stored in Day_3 and EOT so find samp that #
            # will have fold diff #
            if samp.time != "NA" and samp.time != "Day_0":
                
                # find the sample from the patient on day_0 for this protein #
                for ref_samp in prot.sample_list:
                    if ref_samp.time == "Day_0" and ref_samp.patient_ID == samp.patient_ID:
                        if ref_samp.value:    
                            samp.fold_diff_time = samp.value/ref_samp.value
                        break
                
    return protein_list, error_log

def Parse_Data(protein_list, sex_flag, age_flag, age_cutoff, treatment_flag, time_flag, matrix_flag):
    # start by generating the calculated_value_list per protein #
    for prot in protein_list:
        
        # find which calculated value bin the sample fits in #
        for samp in prot.sample_list:
            found = False
            
            for cv in prot.calculated_values_list:
               
                if sex_flag and samp.sex != cv.sex:
                    continue
                
                if age_flag and samp.age != cv.age: #todo: will need to update for age bin
                    continue
                
                if treatment_flag and samp.treatment != cv.treatment:
                     continue
                 
                if time_flag and samp.time != cv.time:
                     continue
    
                if matrix_flag and samp.matrix != cv.matrix:
                     continue
                
                # if we got here then we matched #
                found = True
                cv.values.append(samp.value)
                if samp.fold_diff_time:
                    cv.fold_diff_time_list.append(samp.fold_diff_time)
                break
                    
            # if not found generate a new bin #
            if not found:
                prot.calculated_values_list.append(calculated_values(samp.sex,
                                            samp.age,
                                            samp.treatment,
                                            samp.time,
                                            samp.matrix,
                                            [samp.value]))
                if samp.fold_diff_time:
                    prot.calculated_values_list[-1].fold_diff_time_list = [samp.fold_diff_time]
        
        # ones we've finished collecting values for a single protein, we can calculate values per bin #
        for cv in prot.calculated_values_list:
            cv.mean = np.mean(cv.values)
            cv.std = np.std(cv.values)
            if cv.fold_diff_time_list:
                cv.fold_diff_time_mean = np.mean(cv.fold_diff_time_list)
                cv.fold_diff_time_std = np.std(cv.fold_diff_time_list)
            
        # Once values per treatment type (drug vs placebo) have been found, 
        # multi calc value variables can be calculated. For treatment fold #
        # diffs (drug/placebo) have to find reference calc values. #
        # start by finding a calc value that will have fold diff populated. #
        # drug/placebo fold diffs are stored in the drug calc value #
        for cv in prot.calculated_values_list:
            if cv.time != "NA" and cv.treatment == "LSALT":
                
                # find the calc value that is the placebo at the same time for this protein #
                for ref_cv in prot.calculated_values_list:
                    if ref_cv.treatment == "P" and ref_cv.time == cv.time:
                        
                        if ref_cv.mean:
                            cv.fold_diff_treatment_mean = cv.mean/ref_cv.mean
                        cv.ttest = ss.ttest_ind(cv.values, ref_cv.values, equal_var=False).pvalue
                        cv.mw = ss.mannwhitneyu(cv.values, ref_cv.values).pvalue
                        
                        if cv.time !="Day_0":
                            cv.fold_diff_time_ttest = ss.ttest_ind(cv.fold_diff_time_list, ref_cv.fold_diff_time_list, equal_var=False).pvalue
                            cv.fold_diff_time_mw = ss.mannwhitneyu(cv.fold_diff_time_list, ref_cv.fold_diff_time_list).pvalue
                            
                            cv.fold_diff_treatment_time_mean = cv.fold_diff_time_mean/ref_cv.fold_diff_time_mean
                            
def create_array(protein_list, sex_flag, age_flag, age_cutoff, treatment_flag, time_flag, matrix_flag, order):
    # function assumes time_flag and treatment_flag are true #
    assert time_flag
    assert treatment_flag
    
    # associate category with what order they should be sorted in #
    # todo: probably will want to shift this up so that options aren't hard coded #
    order_list = [order_ob("sex",order[0], ["M", "F"]),
                  order_ob("age",order[1],[69]),
                  order_ob("treatment",order[2], ["LSALT","P"]),
                  order_ob("time",order[3], ["Day_0","Day_3", "EOT"]),
                  order_ob("matrix",order[4], ["EDTA Plasma", "Serum"])]

    order_list.sort(reverse = True, key=lambda order_ob: order_ob.value)
    
    easy_array = []
    
    for prot in protein_list:
        
        # start protein row #
        row = [prot.protein_ID, prot.name, prot.symbol, prot.uniprot]
        
        # generate baseline list of active options #
        active_options = []
        for each in order_list:
            active_options.append(each.options[0])
        
        # generate a protein row #
        if order_list[0].value:
            order_cv(order_list, 0, active_options, prot, row)
    
        # add row to array #
        easy_array.append(row)
    
    easy_array = np.array(easy_array)    
    return easy_array

def order_cv(order_list, idx, active_options, prot, row):
    
    # look through options at current depth (idx) or order list #
    for each in order_list[idx].options:
        
        # update active option array #
        active_options[idx] = each
        
        # check if we need to go further down the list #
        if idx+1<len(order_list) and order_list[idx+1].value:
            order_cv(order_list,idx+1,active_options,prot,row)
        
        else:
            # find the calculated calues object that matches active options #
            for cv in prot.calculated_values_list:
                match = True
                i=0
                # have to walk throuh order_list because active_options is just a list of strings #
                for option in order_list:
                    if option.value:
                        if option.typ == "sex":
                            if cv.sex != active_options[i]:
                                match = False
                                break
                            i=i+1
                            continue
                        
                        if option.typ == "age":
                            if cv.age != active_options[i]:
                                match = False
                                break
                            i=i+1
                            continue
                        
                        if option.typ == "treatment":
                            if cv.treatment != active_options[i]:
                                match = False
                                break
                            i=i+1
                            continue
                        
                        if option.typ == "time":
                            if cv.time != active_options[i]:
                                match = False
                                break
                            i=i+1
                            continue
                        
                        if option.typ == "matrix":
                            if cv.matrix != active_options[i]:
                                match = False
                                break
                            i=i+1
                            continue
                
                if match:
                    row.append(cv.mean)
                    row.append(cv.std)
                    if cv.ttest:
                        row.append(cv.ttest)
                        row.append(cv.mw)
                    
                    if cv.fold_diff_treatment_mean:
                        row.append(cv.fold_diff_treatment_mean)
                        
                    if cv.fold_diff_time_mean:
                        row.append(cv.fold_diff_time_mean)
                        row.append(cv.fold_diff_time_std)
                        
                    if cv.fold_diff_time_ttest:
                        row.append(cv.fold_diff_time_ttest)
                        row.append(cv.fold_diff_time_mw)
                    
                    if cv.fold_diff_treatment_mean and cv.fold_diff_time_mean:
                        row.append(cv.fold_diff_treatment_time_mean)

def check(direc,easy_array):
    test_array=np.array(pd.read_csv(open(direc + 'test_results/array_results.csv', 'r'),header=None))
    error_list=[]
    for i in range(0,len(test_array)):
        for j in range(0,len(test_array[0])):
            test_pass=True
            try:
                if abs(float(easy_array[i][j]) - float(test_array[i][j]))>0.0001 or (m.isnan(float(test_array[i][j])) == False and m.isnan(float(easy_array[i][j]) - float(test_array[i][j]))):
                    test_pass=False
            except:
                if easy_array[i][j] != test_array[i][j]:
                    test_pass=False
            
            if not test_pass:    
                error_list.append("error at easy_array i = " + str(i)+", j = " + str(j) + " expected " + str(test_array[i][j]) + " got " + str(easy_array[i][j]))
    
    if len(test_array) != len(easy_array) or len(test_array[0]) != len(easy_array[0]):
        error_list.append("array size missmatch")
    
    if len(error_list)==0:
        error_list.append("all pass")
    
    return error_list

def main(sex, age, treatment, time, matrix, order):
    protein_list,error_log = Read_Data(directory,'data')
    Parse_Data(protein_list,sex,age,0,treatment,time,matrix)
    easy_array = create_array(protein_list,sex,age,0,treatment,time,matrix,order)
    
    return easy_array,error_log

def test(sex, age, treatment, time, matrix, order):
    protein_list,error_log = Read_Data(directory,'test')
    Parse_Data(protein_list,sex,age,0,treatment,time,matrix)
    easy_array = create_array(protein_list,sex,age,0,treatment,time,matrix,order)
    errors=check(directory, easy_array)
    return errors, easy_array, error_log


print("enter Ture or False for sex, age, treatment, time, matrix. Then enter order array (higher number is lower down on the excel sheet).\n")
# e,a, el = test(False, False, True, True, False, [0,0,1,2,0])
easy_array, error_log = main(False, False, True, True, False, [0,0,1,2,0])    