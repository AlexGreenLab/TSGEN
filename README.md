
# Toehold Switch Generator

The following code enables the automated design, creation, and selection of toehold switches when provided in an excel input file by the user. Toehold switches provide a means to prevent translation of mRNA. Repression is disabled in the presence of a cognate RNA trigger. 

## Installation

The Toehold Design software runs in MATLAB and uses the NUPACK Software Suite. It is recommended that the user have the Parallel Computing MATLAB Toolbox installed for faster design evaluation.

The software here is compatible with versions 3.1 and 3.2 of NUPACK. Follow the installation instructions at [NUPACK.org](http://www.nupack.org/) to install it. 


## Input Sequence Guidelines

The software assumes the user creates an .csv file in the subfolder "/input" with the requisite sequence information. design name, the wild-type target sequence, the mutated target sequence, and the first 29 bases of the output gene. For the target sequence input, the mutation position should be at 15 bases away from the 5’ end of the transcript. Sequences should be specified in the 5’ to 3’ direction.

By default, the code will search for an input file named "design_input_file.csv".

## Running the Software


To run the Toehold Switch Design Software, a user will first need to navigate to the software folder with the main function "toehold_switch_design_run.m" and provide a csv file with the required sequence content outlined above. 

Once the input sequences have been established, the user will want to run the main function:

```matlab
toehold_switch_design_run()

```
The function can be run in its default configuration with no inputs. In this case the software will evauluate the sequences in the file "design_input_file.csv" in parallel, producing six designs each of the Series A (TSgen2) and Series B (EAR2_H22) toeholds for the targets and temperature conditions specified in the input file. Otherwise, users can specify the number of designs to return and the input file by calling the function as:

```
toehold_switch_design_run(n, inputFile)
```
in which `n` is the number of designs and `inputFile` is a string naming the target input file. The input file must be located in an /input folder in the location that the code is running.

The user can specify additional arguments, up to the six options described below, to modify the run conditions:
```
toehold_switch_design_run(num_designs, input_file, options)
```
```num_designs``` specifies the number of toehold switches to return. Set to 6 by default. 

`input_file` is a text string specifying the csv file with sequences required for designs. The file must still be within an 'input' sub-folder from the working directory. Set to 'design_input_file.csv' by default.

`SeriesA` should be set to 1 if you want to generate Series A designs (Toehold Switch generation 2, Green et al., 2014).  Set to 0 if you do not wish the software to return this type of design. Set to 1 by default.

`SeriesB` should be set to 1 if you want to generate Series B designs (Ma et al., 2018 paper). Set to 0 if you do not wish the software to return this type of design. Set to 1 by default

`Parallel` should be set to 1 if you want to run on your local computer using multiple cores. Otherwise, set to 0. Set to 1 by default.

`Antisense` can be set to 1 to create antisense toeholds for the targets provided. Set to 0 by default. 
 
If no inputs are provided, it will run with the default values 
listed below. 
 
Default input values:
```matlab
num_designs = 6
input_file = 'design_input_file.csv';
SeriesA = 1; 
SeriesB = 1; 
PARALLEL = 1; 
ANTISENSE = 0;
```
An example syntax might be:
`toehold_switch_design_run(8,'my_sequences.csv','SeriesB',0,'Antisense',1)`
This returns 8 designs of the SeriesA designs for each input sequence in the 'my_sequences.csv' file, as well as their antisense versions. It returns 0 SeriesB designs for any of the sequences.

###############
Contributing: Larger projects often have sections on contributing to their project, in which contribution instructions are outlined. Sometimes, this is a separate file. If you have specific contribution preferences, explain them so that other developers know how to best contribute to your work. To learn more about how to help others contribute, check out the guide for setting guidelines for repository contributors.

###############
Credits: Include a section for credits in order to highlight and link to the authors of your project.

###############	
License: Finally, include a section for the license of your project. For more information on choosing a license, check out GitHub’s licensing guide!

