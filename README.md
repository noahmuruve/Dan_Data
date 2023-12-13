# Dan_Data
A data parser for statistical analysis of protein concentrations for a lab study.

## Protein_Analysis.py
Protein_Analysis.py is the main scripot in this project. It was developed in Spyder and was not setup to build an exe. To use, run the code in the IDE's kernal and the output will be saved as the following variables:
- easy_array: the actual data. The columns titles match that of the test_results excel table
- error_log: a log of any time data was not readible for whatever reason

This code was run through test as can be seen in the test() function directly above main(). test() works identically to main, except that it parses a smaller data set located in the test folder, and also outputs an additional error log which compares the output data (easy_array in main) to a truth array saved in the test_results folder. The inputs to main() and test() are identical. When the script is first run, instructions on how to set variables are printed in the IDE's kernal. 

## array_results.xlsx
A working excel sheet used to generate and validate the data in the test_results folder.

## Folder Info
### Archive
Contains archived code and data.

### Commit _*_Output
Contains the log files output by a specific version of Protein_Analysis.py with associated input data.

### data
Contains input data

### key
Contains the patient->sample ID key.

### test
Contains test input data

### test_results
Contains the expected results of the test run. This array was manually generated and validated.
