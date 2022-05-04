# Toehold Switch Generator

The following code enables the automated design, analysis, and computational screening of toehold switches against target sequences of interest. Toehold switches provide a means to activate translation of a user-specified output protein in the presence of a cognate target RNA. Target RNA sequences are specified in an input file supplied by the user.

## Installation

The Toehold Switch Design Software runs in MATLAB version 2019b and higher and uses the NUPACK Software Suite. It is recommended that the user have the Parallel Computing MATLAB Toolbox installed for faster design evaluation.

The software here is compatible with versions 3.1 and 3.2 of NUPACK. Follow the installation instructions at [NUPACK.org](http://www.nupack.org/) to install it. 


## Input Sequence Guidelines

The software assumes that the user creates a .csv file in the subfolder "/input" with the requisite information, including the target name, the intended operating temperature(s) of the toehold switch, and the name and sequence of the output protein. A sample input file named "design_input_file.csv" is provided.

The target sequence is specified in the file using the columns "Outer Target" and "Inner Target". The "Outer Target" column should contain the full sequence of the target to be detected or the full amplicon generated after amplification and/or transcription. The software can simulate the folding of this outer target sequence when scoring different toehold switch designs. The "Inner Target" column is used to specify the sequence region of the outer target to which the toehold switch will bind. Typically, the inner target will exclude primer binding sites on the outer target to avoid potential interference between the toehold switch and the amplification primers. If the toehold switches will not be used with primers, the inner and outer targets can have the same sequence. Computational time for the design process will increase as target sequences get longer.

The intended operating temperature of the toehold switch is specified in the "Temperature" column. More than one temperature can be specified by enclosing the list of temperatures in square brackets separated by commas (e.g. `[25,30,37]`).

During the design process, the algorithm will use up to the first 29-nt of the output gene "Output Sequence" when designing and scoring the toehold switches.

## Running the Software

To run the Toehold Switch Design Software, a user will first need to navigate to the software folder with the main function "toehold_switch_design_run.m" and provide a .csv file with the required information content outlined above. 

Once the input sequences have been established, the user will want to run the main function:

```matlab
toehold_switch_design_run()

```
The function can be run in its default configuration with no inputs. In this case the software will evauluate the sequences in the file "design_input_file.csv" in parallel, producing six designs each of the Series A and Series B toeholds for the targets and temperature conditions specified in the input file. Otherwise, users can specify the number of designs to return and the input file by calling the function as:

```
toehold_switch_design_run(n, inputFile)
```
in which `n` is the number of designs to generate and `inputFile` is a string naming the target input file. The input file must be located in an /input folder in the location that the code is running.

The user can specify additional arguments, up to the six options described below, to modify the run conditions:
```
toehold_switch_design_run(num_designs, input_file, options)
```
```num_designs``` specifies the number of toehold switches to return. Set to 6 by default. 

`input_file` is a text string specifying the csv file with sequences required for designs. The file must still be within an 'input' sub-folder from the working directory. Set to 'design_input_file.csv' by default.

`SeriesA` should be set to 1 if you want to generate Series A designs (Green et al., Cell 2014). Set to 0 if you do not wish the software to return this type of design. Set to 1 by default.

`SeriesB` should be set to 1 if you want to generate Series B designs (Pardee et al., Cell 2016; Ma et al., OUP Syn. Biol. 2018). Set to 0 if you do not wish the software to return this type of design. Set to 1 by default

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
