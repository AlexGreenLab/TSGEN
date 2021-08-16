
# Toehold Switch Generator

The following code enables the automated design, creation, and selection of toehold switches when provided in an excel input file by the user. Toehold switches provide a means to prevent translation of mRNA. Repression is disabled in the presence of a cognate RNA trigger. 

## Installation

The Toehold Design software runs in MATLAB and uses the NUPACK Software Suite. It is recommended that the user have the Parallel Computing MATLAB Toolbox installed for faster design evaluation.

The software here is compatible with versions 3.1 and 3.2 of NUPACK. Follow the installation instructions at [NUPACK.org](http://www.nupack.org/) to install it. 


## Input Sequence Guidelines

The software assumes the user creates an .csv file in the subfolder "/input" with the requisite sequence information. design name, the wild-type target sequence, the mutated target sequence, and the first 29 bases of the output gene. For the target sequence input, the mutation position should be at 15 bases away from the 5’ end of the transcript. Sequences should be specified in the 5’ to 3’ direction.

By default, the code will search for an input file named "design_input_file_old.csv".

## Running the Software


To run the Toehold Switch Design Software, a user will first need to navigate to the software folder with the main function "toehold_switch_design_run_GCM.m" and provide a csv file with the required sequence content outlined above. 

Once the input sequences have been established, the user will want to run the main function:

```matlab
toehold_switch_design_run_GCM()

```
The function can be run in its default configuration with no inputs. In this case the software will evauluate the sequences in the file "design_input_file_old.csv" in parallel, producing six designs each of the Series B (EAR2_H22) and Series C (TSgen2) toeholds for the targets and temperature conditions specified in the input file. Users can specify a diferrent number of designs to return by calling the function as:

```
toehold_switch_design_run_GCM(n)
```
in which `n` is the number of designs.  

Otherwise, the user can provide five function inputs to modify the run conditions as seen below:
```
toehold_switch_design_run_GCM(num_designs, TSGEN2, EAR2_H22, input_file, PARALLEL)
```
```num_designs``` specifies the number of TSgen2 and EAR_H22 toehold switches to return 

`TSGEN2` should be set to 1 if you want to generate Series C designs (Toehold Switch generation 2, Ma et al., OUP Synbio parameters).  Set to 0 if you do not wish the software to return this type of design.

`EAR2_H22` should be set to 1 if you want to generate Series B designs (Zika 2016 paper). Set to 0 if you do not wish the software to return this type of design.

`input_file` is a text string specifying the csv file with sequences required for designs. The file must still be within 'input' folder.

`PARALLEL` should be set to 1 if you want to run on your local computer using multiple cores. Otherwise, set to 0.
 
If no inputs are provided, it will run with the default values 
listed below. If one input is provided it will modify the num_designs 
value only and keep the other default values. If all five inputs are
listed, it will modify each according to the input guide above. 
 
Default input values:
```matlab
num_designs = 6
TSGEN2 = 1; 
EAR2_H22 = 1; 
input_file = 'design_input_file_old.csv';
PARALLEL = 1; 
```


###############
Contributing: Larger projects often have sections on contributing to their project, in which contribution instructions are outlined. Sometimes, this is a separate file. If you have specific contribution preferences, explain them so that other developers know how to best contribute to your work. To learn more about how to help others contribute, check out the guide for setting guidelines for repository contributors.

###############
Credits: Include a section for credits in order to highlight and link to the authors of your project.

###############	
License: Finally, include a section for the license of your project. For more information on choosing a license, check out GitHub’s licensing guide!

