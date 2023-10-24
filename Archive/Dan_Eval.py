import os
import pandas as pd
import math as m
import numpy as np
import scipy

##############################################################################
directory='C:/Users/noahm/Desktop/projects/Dan_Data/'  
##############################################################################

class fileinfo:
    def __init__(self,patient_list, protein_list):
        self.patient_list=patient_list
        self.protein_list=protein_list
        
class patient:
    def __init__(self,ID,sample_list,sex,age,drug,fold_diff_plus_array):
        self.ID=ID
        self.sample_list=sample_list
        self.sex=sex
        self.age=age
        self.drug=drug
        self.fold_diff_plus_array = fold_diff_plus_array
        
class sample:
    def __init__(self,SampleId = "", SampleTime = "", SampleMatrix = "", Proteins=[]):
        self.SampleId = SampleId
        self.SampleTime = SampleTime
        self.SampleMatrix =SampleMatrix
        self.Proteins = Proteins

class protein_sample:
    def __init__(self,ID,data,name, symbol, uniprot,value):
        self.ID=ID
        self.data=data
        self.name=name
        self.symbol=symbol
        self.uniprot=uniprot
        self.value = value

# class protein_sample:
#     def __init__(self, time, values):
#         self.time=time
#         self.values=values

class patient_fold_diff_plus:
    def __init__(self, protein_ID, test_day, fold_diff, ttest, mw):
        self.protein_ID = protein_ID
        self.test_day = test_day
        self.fold_diff = fold_diff
        self.ttest = ttest
        self.mw = mw

class protein_interval:
    def __init__ (self, ):
        self.ID = ID
        self.sex = sex
        self.age = age
        self.drug = drug
        self.time = time
        self.matrix = matrix
        self.name=name
        self.symbol=symbol
        self.uniprot=uniprot
        self.values=values
class protein_data:
    def __init__(self, 

                 meanP = 0,
                 meanD = 0,
                 meanMP = 0,
                 meanMD = 0,
                 meanFP = 0,
                 meanFD = 0,
                 
                 sumP = 0,
                 sumD = 0,
                 sumMP = 0,
                 sumMD = 0,
                 sumFP = 0,
                 sumFD = 0,
                 
                 numP = 0,
                 numD = 0,
                 numMP = 0,
                 numMD = 0,
                 numFP = 0,
                 numFD = 0,
                 
                 fold_diff = 0,
                 fold_diff_M = 0,
                 fold_diff_F = 0,
                 
                 stdD = 0,
                 stdP = 0,
                 stdMP = 0,
                 stdMD = 0,
                 stdFP = 0,
                 stdFD = 0,
                 
                 ttestD = 0,
                 ttestP = 0,
                 ttestMD = 0,
                 ttestMP = 0,
                 ttestFP = 0,
                 ttestFD = 0,
                 
                 mwD = 0,
                 mwP = 0,
                 mwMD = 0,
                 mwMP = 0,
                 mwFP = 0,
                 mwFD = 0):
        
        
        self.meanP = meanP
        self.meanD = meanD
        self.meanMP = meanMP
        self.meanMD = meanMD
        self.meanFP = meanFP
        self.meanFD = meanFD
        
        self.sumP = sumP
        self.sumD = sumD
        self.sumMP = sumMP
        self.sumMD = sumMD
        self.sumFP = sumFP
        self.sumFD = sumFD
       
        self.numP = numP
        self.numD = numD
        self.numMP = numMP
        self.numMD = numMD
        self.numFP = numFP
        self.numFD = numFD
        
        self.fold_diff = fold_diff
        self.fold_diff_M = fold_diff_M
        self.fold_diff_F = fold_diff_F
        
        self.stdD = stdD
        self.stdP = stdP
        self.stdMP = stdMP
        self.stdMD = stdMD
        self.stdFP = stdFP
        self.stdFD = stdFD
        
        self.ttestD = ttestD
        self.ttestP = ttestP
        self.ttestMP = ttestMP
        self.ttestMD = ttestMD
        self.ttestFP = ttestFP
        self.ttestFD = ttestFD
        
        self.mwD = mwD
        self.mwP = mwP
        self.mwMP = mwMP
        self.mwMD = mwMD
        self.mwFP = mwFP
        self.mwFD = mwFD
# read each file in the directory #
def Read_Data(direc, datasub):
    # --------------------------------------------------------------------
    ## read key ##
    
    key_file="key/AB002 Subjects - Demographics.csv"
    key=pd.read_csv(open(direc + key_file, 'r'),header=None)
    
    # --------------------------------------------------------------------
    ## generate patient_list for file ##
    file_list = []
    for file in os.listdir(direc + datasub):
        patient_list = []
        protein_list = []
        raw_file = pd.read_csv(open(direc + datasub + '/' + file, 'r'))
        
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
        
        # populate patient data #
        for i in range(sample_title_row+1,len(raw_file)):
            
            # seperate patient_ID from sample time stamp #
            temp = raw_file.SampleId[i].split()
            raw_file.SampleId[i]=sample(temp[0].strip(),temp[1:])
            patient_ID=raw_file.SampleId[i].SampleId
            timevalid=False
            for item in raw_file.SampleId[i].SampleTime:
                if item.startswith("D"):
                    day=item[1:]
                    timevalid=True
                    if float(day) < 3:
                        raw_file.SampleId[i].SampleTime = "Day_0" #zero
                    elif float(day) == 3:
                        raw_file.SampleId[i].SampleTime = "Day_3" #3
                    elif float(day) > 3:
                        raw_file.SampleId[i].SampleTime = "EOT" #EOT
                    else:
                        print("Error, sample time not valid",patient_ID, raw_file.SampleId[i].SampleTime)
                        
                elif item == "EOT":
                    raw_file.SampleId[i].SampleTime = "EOT"
                    timevalid=True
            
            if not timevalid:
                print("Error, sample time not valid for",patient_ID, raw_file.SampleId[i].SampleTime)
                
            # add sample to existing patient or add patient #
            added = False
            for subject in patient_list:
                if subject.ID == patient_ID:
                    subject.sample_list.append(raw_file.loc[i])
                    added=True
                    break
            if not added:
                
                # if creating a new patient, find their age, sex, and drug #
                # based on key #
                subfound=False
                for sub in key.iterrows():
                    if patient_ID==sub[1][0]:
                        s=sub[1][1].strip()
                        a=sub[1][2] #this one is int so does not need strip()
                        d=sub[1][3].strip()
                        subfound=True
                        break
                    
                # if not found means its a control or something # 
                if not subfound:
                    s="NA"
                    a="NA"
                    d="NA"
                    
                patient_list.append(patient(patient_ID,[raw_file.loc[i]],s,a,d))
                
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
                title = row.SampleId.SampleId
                for each in row.SampleId.SampleTime:
                    title= title + ' ' + each
            elif i != sample_title_row:
                title = raw_file.iat[i,sample_title_column]
            title_list.append(title)
            i=i+1
        
        raw_file.index = title_list
        
        # --------------------------------------------------------------------
        # generate protein List #  
        for i in range (sample_title_column+1, len(raw_file.columns)):
            protein_list.append(protein(raw_file.columns[i],
                            [protein_sample("Day_0", protein_data([],[],[],[],[],[],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),
                             protein_sample("Day_3",protein_data([],[],[],[],[],[],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),
                             protein_sample("EOT",protein_data([],[],[],[],[],[],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))],
                            raw_file.loc['TargetFullName'][i],
                            raw_file.loc['EntrezGeneSymbol'][i],
                            raw_file.loc['UniProt'][i])) 
        
        ## Additional test for key association ##
        for each in patient_list:
            if each.age == "NA":
                print("sample ID not found on key:",each.ID)
    
        
        file_list.append(fileinfo(patient_list, raw_file, protein_list))
        
    return file_list

def Parse_Data(files):
    
    # get down to the sample level #
    for f in files:
        for patient in f.patient_list:
            for sample in patient.sample_list:
                
                                        
                if sample.SampleId.SampleTime == "Day_0":
                    ii=0
                elif sample.SampleId.SampleTime == "Day_3":
                    ii=1
                else:
                    ii=2
                
                i=0
                j=0
                after_SeqId = False
                for each in sample:

                    # we now know we are looking at actual data. see elif #
                    # below #
                    if after_SeqId:
                        
                        # log2() all data #
                        #sample[i]= m.log2(float(each))
                        sample[i]=float(each)
                                                
                        if patient.sex != "NA":
                            # add data to protein data summary #
                            if patient.sex == "M":
                                if patient.drug == "P":
                                    f.proteins[j].data[ii].values.sumMP+=sample[i]
                                    f.proteins[j].data[ii].values.numMP+=1
                                    f.proteins[j].data[ii].values.valuesP.append(sample[i])
                                    f.proteins[j].data[ii].values.valuesMP.append(sample[i])
                                else:
                                    f.proteins[j].data[ii].values.sumMD+=sample[i]
                                    f.proteins[j].data[ii].values.numMD+=1
                                    f.proteins[j].data[ii].values.valuesD.append(sample[i])
                                    f.proteins[j].data[ii].values.valuesMD.append(sample[i])
                            elif patient.sex =="F":
                                if patient.drug =="P":    
                                    f.proteins[j].data[ii].values.sumFP+=sample[i]
                                    f.proteins[j].data[ii].values.numFP+=1
                                    f.proteins[j].data[ii].values.valuesP.append(sample[i])
                                    f.proteins[j].data[ii].values.valuesFP.append(sample[i])
                                else:
                                    f.proteins[j].data[ii].values.sumFD+=sample[i]
                                    f.proteins[j].data[ii].values.numFD+=1
                                    f.proteins[j].data[ii].values.valuesD.append(sample[i])
                                    f.proteins[j].data[ii].values.valuesFD.append(sample[i])
                        j=j+1
                            
                    # This makes sure only actual data is looked at and #
                    # preamble is ignored #
                    elif sample.index[i] == 'SeqId':
                        after_SeqId = True
                    i=i+1


# unsused function. Purpose: if user wants to separate data by file #    
def Seperate_Files(files):
    for f in files:
        for protein in f.proteins:
            for ii in range(0,len(protein.data)):
                # calculate standard deviation #
                protein.data[ii].values.stdD = np.std(protein.data[ii].values.valuesD)
                protein.data[ii].values.stdP = np.std(protein.data[ii].values.valuesP)
                
                # calculate gendered standard deviation #
                protein.data[ii].values.stdMD = np.std(protein.data[ii].values.valuesMD)
                protein.data[ii].values.stdMP = np.std(protein.data[ii].values.valuesMP)
                protein.data[ii].values.stdFD = np.std(protein.data[ii].values.valuesFD)
                protein.data[ii].values.stdFP = np.std(protein.data[ii].values.valuesFP)
                
                #calculate means #
                if protein.data[ii].values.numMP or protein.data[ii].values.numFP :
                    protein.data[ii].values.meanP=(protein.data[ii].values.sumMP+protein.data[ii].values.sumFP)/(protein.data[ii].values.numMP+protein.data[ii].values.numFP)
                if protein.data[ii].values.numMD or protein.data[ii].values.numFD :
                    protein.data[ii].values.meanD=(protein.data[ii].values.sumMD+protein.data[ii].values.sumFD)/(protein.data[ii].values.numMD+protein.data[ii].values.numFD)                
                
                # calculate gendered means #
                if protein.data[ii].values.numMP:
                    protein.data[ii].values.meanMP=protein.data[ii].values.sumMP/protein.data[ii].values.numMP
                if protein.data[ii].values.numMD:
                    protein.data[ii].values.meanMD=protein.data[ii].values.sumMD/protein.data[ii].values.numMD
                if protein.data[ii].values.numFP:
                    protein.data[ii].values.meanFP=protein.data[ii].values.sumFP/protein.data[ii].values.numFP
                if protein.data[ii].values.numFD:
                    protein.data[ii].values.meanFD=protein.data[ii].values.sumFD/protein.data[ii].values.numFD
                
                # Calculate fold Difference #
                if protein.data[ii].values.meanP:
                    protein.data[ii].values.fold_diff= protein.data[ii].values.meanD/protein.data[ii].values.meanP
                
                # Calculate gendered fold Difference #
                if protein.data[ii].values.meanMP:
                    protein.data[ii].values.fold_diff_M = protein.data[ii].values.meanMD/protein.data[ii].values.meanMP
                if protein.data[ii].values.meanFP:
                    protein.data[ii].values.fold_diff_F = protein.data[ii].values.meanFD/protein.data[ii].values.meanFP
        
def Combine_Data(files):
    
    # go through each protein and match it to proteins from other files #
    multifile_proteins=[]
    for f in files:
        for prot in f.proteins:
            
            added = False
            for mf_prot in multifile_proteins:
                if mf_prot.ID == prot.ID:
                    for ii in range(0,len(prot.data)):
                        
                        # if found combine data #
                        mf_prot.data[ii].values.valuesD += prot.data[ii].values.valuesD
                        mf_prot.data[ii].values.valuesP += prot.data[ii].values.valuesP
                        
                        mf_prot.data[ii].values.valuesMD += prot.data[ii].values.valuesMD
                        mf_prot.data[ii].values.valuesMP += prot.data[ii].values.valuesMP
                        mf_prot.data[ii].values.valuesFD += prot.data[ii].values.valuesFD
                        mf_prot.data[ii].values.valuesFP += prot.data[ii].values.valuesFP
                        
                        mf_prot.data[ii].values.sumP += prot.data[ii].values.sumP
                        mf_prot.data[ii].values.sumD += prot.data[ii].values.sumD
                        
                        mf_prot.data[ii].values.sumMP += prot.data[ii].values.sumMP
                        mf_prot.data[ii].values.sumMD += prot.data[ii].values.sumMD
                        mf_prot.data[ii].values.sumFP += prot.data[ii].values.sumFP
                        mf_prot.data[ii].values.sumFD += prot.data[ii].values.sumFD
                        
                        mf_prot.data[ii].values.numMP += prot.data[ii].values.numMP
                        mf_prot.data[ii].values.numMD += prot.data[ii].values.numMD
                        mf_prot.data[ii].values.numFP += prot.data[ii].values.numFP
                        mf_prot.data[ii].values.numFD += prot.data[ii].values.numFD
                    
                    added = True
                    break
                
            # if not found then add new protein to list #
            if added == False:
                multifile_proteins.append(prot)
            

    for protein in multifile_proteins:
        for ii in range(0,len(protein.data)):
            # calculate standard deviation #
            protein.data[ii].values.stdD = np.std(protein.data[ii].values.valuesD)
            protein.data[ii].values.stdP = np.std(protein.data[ii].values.valuesP)
            
            # calculate gendered standard deviation #
            protein.data[ii].values.stdMD = np.std(protein.data[ii].values.valuesMD)
            protein.data[ii].values.stdMP = np.std(protein.data[ii].values.valuesMP)
            protein.data[ii].values.stdFD = np.std(protein.data[ii].values.valuesFD)
            protein.data[ii].values.stdFP = np.std(protein.data[ii].values.valuesFP)
            
            # calculate mean #
            if protein.data[ii].values.numMP or protein.data[ii].values.numFP :
                protein.data[ii].values.meanP=(protein.data[ii].values.sumMP+protein.data[ii].values.sumFP)/(protein.data[ii].values.numMP+protein.data[ii].values.numFP)
            if protein.data[ii].values.numMD or protein.data[ii].values.numFD :
                protein.data[ii].values.meanD=(protein.data[ii].values.sumMD+protein.data[ii].values.sumFD)/(protein.data[ii].values.numMD+protein.data[ii].values.numFD)  
            
            # calculate gendered mean #
            if protein.data[ii].values.numMP:
                protein.data[ii].values.meanMP=protein.data[ii].values.sumMP/protein.data[ii].values.numMP
            if protein.data[ii].values.numMD:
                protein.data[ii].values.meanMD=protein.data[ii].values.sumMD/protein.data[ii].values.numMD
            if protein.data[ii].values.numFP:
                protein.data[ii].values.meanFP=protein.data[ii].values.sumFP/protein.data[ii].values.numFP
            if protein.data[ii].values.numFD:
                protein.data[ii].values.meanFD=protein.data[ii].values.sumFD/protein.data[ii].values.numFD
        
            
            # calculate fold difference #
            if protein.data[ii].values.meanP:
                protein.data[ii].values.fold_diff= protein.data[ii].values.meanD/protein.data[ii].values.meanP
            
            # calculate gendered fold Difference #
            if protein.data[ii].values.meanMP:
                protein.data[ii].values.fold_diff_M = protein.data[ii].values.meanMD/protein.data[ii].values.meanMP
            if protein.data[ii].values.meanFP:
                protein.data[ii].values.fold_diff_F = protein.data[ii].values.meanFD/protein.data[ii].values.meanFP
            
    return multifile_proteins

def create_array(combined_data):
    temp=[]
    for each in combined_data:
        row=[each.ID, each.name]
        for time in each.data:
            row.append(time.values.meanD)
            row.append(time.values.stdD)
            row.append(time.values.meanP)
            row.append(time.values.stdP)
            row.append(time.values.fold_diff)
            
            if time.time != "Day_0":
                MD_to_O = time.values.meanD/each.data[0].values.meanD
                MP_to_O = time.values.meanP/each.data[0].values.meanP
                ZZ = MD_to_O/MP_to_O 
                row.append(MD_to_O )
                row.append(MP_to_O)
                row.append(ZZ)
                
        temp.append(row)
        
    easy_array=np.array(temp)
    return easy_array

def check(direc,easy_array):
    test_array=np.array(pd.read_csv(open(direc + 'test_results/array_results.csv', 'r'),header=None))
    error_list=[]
    for i in range(0,len(test_array)):
        for j in range(0,len(test_array[0])):
            test_pass=True
            try:
                if abs(float(easy_array[i][j]) - float(test_array[i][j]))>0.0001:
                    test_pass=False
            except:
                if easy_array[i][j] != test_array[i][j]:
                    test_pass=False
            
            if not test_pass:    
                error_list.append("error at easy_array i = " + str(i)+", j = " + str(j) + " expected " + str(test_array[i][j]) + " got " + str(easy_array[i][j]))
    
    if len(error_list)==0:
        error_list.append("all pass")
    
    return error_list
    
    
def main():
    files = Read_Data(directory,'data')
    Parse_Data(files)
    combined_data=Combine_Data(files)
    easy_array = create_array(combined_data)
    
    return easy_array

def test():
    files = Read_Data(directory,'test')
    Parse_Data(files)
    combined_data=Combine_Data(files)
    easy_array = create_array(combined_data)
    
    errors=check(directory, easy_array)
    return errors, easy_array

e,a=test()