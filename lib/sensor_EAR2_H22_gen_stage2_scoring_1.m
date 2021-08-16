dots = char(zeros(1,5000)+'.');
bleft = char(zeros(1,5000)+'(');
bright = char(zeros(1,5000)+')');
Ns = char(zeros(1,5000)+'N');

design_type = 'SeriesB';
[~,~,~,] = mkdir(design_type);
RBS_seq = 'AGAGGAGA';
probe_index = 12;
target_index = 13;
discard_window = 2;
num_designs = num_designs2; 

cd(design_type);
for c1 = 1:length(mRNA_names)
    mkdir2(mRNA_names{c1});
    cd(mRNA_names{c1});
    design_data = loadCSVabsFileRawtoStructure('mRNA_RR_results_complete.csv','');
    design_data(1,:) = [];
    design_spec_data = [];
    norm_design_spec_data = [];
    for c2 = 1:size(design_data,1)
        target_seq = dna2rna_B(design_data{c2,target_index});
        target_len = length(target_seq);
        probe_seq  = dna2rna_B(design_data{c2,probe_index});
        index_RBS = findstr(probe_seq,RBS_seq);
        index_linker = findstr(probe_seq,hairpin_suffix0);
        index_AUG = findstr(probe_seq,'AUG');
        start_codon_index = 0;
        for c3 = 1:length(index_AUG)
            if index_AUG(c3) < index_RBS(1)
                continue;
            end
            if mod((index_linker - index_AUG(c3)),3) == 0
                start_codon_index = index_AUG(c3);
                break;
            end 
        end
        for c3 = 1:size(design_data,2)-2
            design_spec_data(c2,c3) = str2num(design_data{c2,c3});
        end
        if start_codon_index == 0
            design_spec_data(c2,1) = 0;
            fprintf('***** EAR2_H22 %s ERROR: Suitable start codon not found in design %d. *****\n',mRNA_names{c1},c2);
            continue;
        end
        aa_seq = nt2aa(probe_seq(start_codon_index:end));
        if ~isempty(findstr(aa_seq,'*'))
            design_spec_data(c2,1) = 0;
            fprintf('EAR2_H22 %s: Design %d of %d has in-frame stop codon. Continuing.\n',mRNA_names{c1},c2,size(design_data,1));
            continue;
        end
    end
    %variable sets:
    %design_spec_data
    % 1:'Position'
    % 2:'Probe Defect'
    % 3:'Probe Min. Defect'
    % 4:'Toehold Defect'
    % 5:'Target Defect',
    % 6:'mRNA Defect',
    % 7:'Active mRNA Defect'
    % 8:'DeltaG Active mRNA'
    % 9: 'Probe/mRNA Binding Defect'
    % 10: 'Probe/mRNA Active Translation Defect'
    % 11:'Binding Fraction'
    
    discard_set = [];
    keep_set = [];
    for c2 = 1:size(design_spec_data,2)
        norm_design_spec_data(:,c2) = design_spec_data(:,c2)./max(design_spec_data(:,c2));
        if design_spec_data(c2,1) == 0
            discard_set(end+1,1) = c2;
        else
            keep_set(end+1,1) = c2;
        end
    end
    
    design_score_original = ((2*norm_design_spec_data(:,3) + 5*norm_design_spec_data(:,4) + 2*norm_design_spec_data(:,6) + 4*norm_design_spec_data(:,7) + (1-design_spec_data(:,11))).*(design_spec_data(:,11) > 0.9)).*(design_spec_data(:,1) ~= 0) + (design_spec_data(:,1) == 0)*1e6;
    [sorted_design_score_original,sorted_indices1] = sort(design_score_original);
    
    design_score = (-93.2*design_spec_data(:,2) - 43.3*design_spec_data(:,7) - 22.1*design_spec_data(:,6) - 9.4*design_spec_data(:,5) + 61.3) - (design_spec_data(:,1) == 0)*1e6;
    
    [sorted_design_score,sorted_indices2] = sort(design_score,'descend');
    combined_data = [];
    labels = {'Score with Multiple Targets','Score with Base Target','Rank','Position','Probe Defect','Probe Min. Defect','Toehold Defect','Target Defect','Base mRNA Defect','Active mRNA Defect','DeltaG Active mRNA','Probe/mRNA Binding Defect','Probe/mRNA Active Translation Defect','Binding Fraction',...
        'Probe Sequence','Target Sequence'};
    end_skip = length(discard_set);
    combined_data = [design_score,design_score_original,design_score_original*0,design_spec_data];
    combined_data_sorted = combined_data(sorted_indices2,:);
    combined_data_sorted(:,3) = [1:size(combined_data_sorted,1)]';
    combined_data_sorted_array = {};
    for c2 = 1:size(combined_data_sorted,1)
        for c3 = 1:size(combined_data_sorted,2)
            combined_data_sorted_array{c2,c3} = num2str(combined_data_sorted(c2,c3));
        end
    end
    seq_data = design_data(:,[probe_index,target_index]);
    seq_data_sorted = seq_data(sorted_indices2,:);
    combined_data_sorted_array = [combined_data_sorted_array,seq_data_sorted];
    combined_data_sorted_array = [labels;combined_data_sorted_array];
    writeCSVFilefromStructure(sprintf('../AAB_%s_%s_all_designs_ranked.csv',design_type,mRNA_names{c1}),combined_data_sorted_array);
    
    discard_window1 = discard_window;
    while 1
        best_peaks = locateNpeaks(min(round(size(design_score,1)/2),num_designs),design_spec_data(:,1),design_score,discard_window1);
        if size(best_peaks,1) < num_designs%4
            discard_window1 = discard_window1 - 1;
            if discard_window1 <= 0
                temp_max = size(design_score,1);%min(round(size(design_score,1)/2),num_designs);
                [temp_sorted,temp_indices] = sort(design_score,'descend');
                best_peaks = [design_spec_data(temp_indices(1:temp_max),1),-temp_sorted(1:temp_max)];
                break;
            end
            fprintf('\tReducing discard window (%d).\n',discard_window1);
        else
            break;
        end
    end
    sorted_indices2 = [];
    for c2 = 1:size(best_peaks,1)
        [~,index] = min(abs(best_peaks(c2,1)-design_spec_data(:,1)));
        sorted_indices2(end+1,1) = index(1);
    end
    combined_data = [];
    labels = {'Score with Multiple Targets','Score with Base Target','Rank','Position','Probe Defect','Probe Min. Defect','Toehold Defect','Target Defect','Base mRNA Defect','Active mRNA Defect','DeltaG Active mRNA','Probe/mRNA Binding Defect','Probe/mRNA Active Translation Defect','Binding Fraction',...
        'Probe Sequence','Target Sequence'};
    
    combined_data = [design_score,design_score_original,design_score_original*0,design_spec_data];
    combined_data_sorted = combined_data(sorted_indices2,:);
    combined_data_sorted(:,3) = [1:size(combined_data_sorted,1)]';
    combined_data_sorted_array = {};
    for c2 = 1:size(combined_data_sorted,1)
        for c3 = 1:size(combined_data_sorted,2)
            combined_data_sorted_array{c2,c3} = num2str(combined_data_sorted(c2,c3));
        end
    end
    seq_data = design_data(:,[probe_index,target_index]);
    seq_data_sorted = seq_data(sorted_indices2,:);
    combined_data_sorted_array = [combined_data_sorted_array,seq_data_sorted];
    combined_data_sorted_array = [labels;combined_data_sorted_array];
    writeCSVFilefromStructure(sprintf('../AAA_%s_%s_top_designs.csv',design_type,mRNA_names{c1}),combined_data_sorted_array);
    cd ../
end
cd ../
