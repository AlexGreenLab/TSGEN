%toehold_switch_design_run accepts sequence inputs to create
% corresponding designs for toehold switch sensor devices
%
%function toehold_switch_design_run(num_designs,input_file,options)
% The function allows for specification of the number of returned designs
% and the sequence input file (must be .csv file located in a '/input' 
% subfolder from the working directory), as well as other run parameters. 
% If no inputs are provided, the following defaults are assumed:
% 
%Input defaults:
% num_designs = 6
% input_file = 'design_input_file.csv'
%
%Optional run parameters and their defaults are described below
%
% SeriesA (default value = 1)
%   Defines if the function will produce Series A designs (Toehold Switch 
%   generation 2 from Green et al. 2014). Set to 0 to not produce this 
%   design type. 
%
%
% SeriesB (default value = 1) 
%   Defines if the function will produce Series B designs (EAR2_H22 from 
%   Ma et al. 2018). Set to 0 to not produce this design type. 
% 
% Parallel (default value = 1)
%   Enables processing with multiple cores on local computer. Set to 0 if 
%   you do not want parallelization enabled or do not have the Parallel
%   Processing Toolbox.
%
% Antisense (default value = 0)
%   Set to 1 if you wish to create toehold switch sensors for the antisense
%   sequences of those provided in the input file.
%
% Example syntax:
% toehold_switch_design_run(8,'my_sequences.csv','SeriesB',0,'Antisense',1)
%   This returns 8 of the SeriesA designs for each input sequence and its 
%   antisense, and 0 of the SeriesB designs for either sequence. 
%
function toehold_switch_design_run(num_designs,input_file,options)
    arguments
        num_designs (1,1) {mustBeNumeric} = 6
        input_file (1,:) {mustBeText} = 'design_input_file.csv'
        options.SeriesA (1,1) {mustBeNumeric} = 1
        options.SeriesB (1,1) {mustBeNumeric} = 1
        options.Parallel (1,1) {mustBeNumeric} = 1
        options.Antisense (1,1) {mustBeNumeric} = 0
    end
    
	warning('off');
	addpath('lib'); %add required helper functions to Matlab path
    


	% Design parameters:

    FINAL_NUM_DESIGNS = num_designs;
    RUN_SERIESA = options.SeriesA; 
    RUN_SERIESB = options.SeriesB; 
    IS_PARALLEL = options.Parallel; 
    RUN_ANTISENSE = options.Antisense;
    
	% hairpin_suffix0 should be a 21-nt sequence and lack in-frame stop codon; 
	% this is the common linker21 used for toehold switches
	hairpin_suffix0 = 'AACCUGGCGGCAGCGCAAAAG'; 
    
    fprintf("Loading targets from: %s\n",input_file);
    input_table = loadCSVabsFileRawtoStructureQuotes(sprintf('input/%s',input_file),''); 
	input_table(1,:) = [];
	num_designs1 = FINAL_NUM_DESIGNS;
	num_designs2 = num_designs1*2;

	%Calling main design, scoring, and selection scripts:
    
    if IS_PARALLEL
        fprintf("Evaluating designs in parallel\n");
    end
    
    if RUN_ANTISENSE
        fprintf("Generating antisense designs for targets\n");
    end    
    
	%Generate Series A designs (Toehold Switch generation 2, Green et al., Cell 2014)
	if RUN_SERIESA
        fprintf("Creating %d SeriesA designs\n", FINAL_NUM_DESIGNS);
		sensor_SeriesA_stage1_1; %generate sequences and performance parameters
		sensor_SeriesA_stage2_scoring_1; %score each design
		sensor_SeriesA_stage3_selection_1; %Select FINAL_NUM_DESIGNS for testing
	end

	%Generate Series B style designs (EAR2_H22, Pardee et al. 2016)
	if RUN_SERIESB
        fprintf("Creating %d SeriesB designs", FINAL_NUM_DESIGNS);
		sensor_SeriesB_gen_stage1_1; %generate sequences and performance parameters
		sensor_SeriesB_gen_stage2_scoring_1; %score each design
		sensor_SeriesB_gen_stage3_selection_1; %Select FINAL_NUM_DESIGNS for testing
	end


end
