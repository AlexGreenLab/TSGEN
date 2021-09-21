hfig = [];
design_type = 'SeriesA';
target_RNA_table = loadCSVabsFileRawtoStructure([design_type,'/ZZZ_original_sequence_info.csv'],'');

linker21 = hairpin_suffix0;
T7_primer = 'GCGCTAATACGACTCACTATAGGG';

cd(design_type)
sel_designs = {}; %name1, designs1, seq, mean score1, name2, seq, designs2, mean score2
mRNA_groups = {};
for c1 = 1:size(target_RNA_table,1)
    curr_name = target_RNA_table{c1,1};
    min_name = curr_name;
    flag = 0;
    for c2 = 1:length(mRNA_groups)
        if issame2(mRNA_groups{c2},min_name)
            flag = 1;
            break;
        end
    end
    if flag == 0
        mRNA_groups{end+1,1} = min_name;
    end
end

design_comp_data = [];
for c1 = 1:length(mRNA_groups)
    file_list = dir(sprintf('AAA_*_%s_top_designs.csv',mRNA_groups{c1}));
    design_set = {};
    design_scores1 = [];
    design_scores2 = [];
    for c2 = 1:length(file_list)
        file_name = file_list(c2).name;
        if ~isempty(findstr(file_name,'_sense_'))
            IS_SENSE = 1;
        else
            IS_SENSE = 0;
        end
        indices = findstr(file_name,'_');
        target_name = file_name(indices(2)+1:indices(end-1)-1);
        [csv,labels] = loadCSVabsFile(file_name,'');
        seq_set = loadCSVabsFileRawtoStructure(file_name,'');
        seq_set(1,:) = [];
        seq_set = seq_set(1:max(num_designs1,num_designs2),end-1:end);
        design_scores1(end+1,:) = [mean(csv(1:num_designs1,1)),min(csv(1:num_designs1,1)),max(csv(1:num_designs1,1))];
        design_scores2(end+1,:) = [mean(csv(1:num_designs2,1)),min(csv(1:num_designs2,1)),max(csv(1:num_designs2,1))];
        design_set(end+1,:) = {target_name,seq_set};
        design_comp_data = [design_comp_data;csv];
    end
    [sorted_design_scores1_mean,indices_mean1] = sort(design_scores1(:,1));
    [~,indices_min1] = sort(design_scores1(:,2));
    [~,indices_max1] = sort(design_scores1(:,3));
    [~,indices_mean2] = sort(design_scores2(:,1));
    [~,indices_min2] = sort(design_scores2(:,2));
    [~,indices_max2] = sort(design_scores2(:,3));
    
    %sort by mean score of num_designs1
    sel_design_set = design_set(indices_mean1(1:min(2,length(indices_mean1))),:);
    sel_targets = {};
    for c2 = 1:size(sel_design_set,1)
        curr_targ = sel_design_set{c2,1};
        for c3 = 1:size(target_RNA_table,1)
            if issame2(curr_targ,target_RNA_table{c3,1})
                sel_targets(end+1,:) = target_RNA_table(c3,:);
                break;
            end
        end
    end
    if size(sel_targets,1) == 1
        sel_designs(end+1,:) = {sel_targets{1,1},sel_targets{1,2},sel_design_set{1,2},sorted_design_scores1_mean(1)};
    else
        sel_designs(end+1,:) = {sel_targets{1,1},sel_targets{1,2},sel_design_set{1,2},sorted_design_scores1_mean(1),sel_targets{2,1},sel_targets{2,2},sel_design_set{2,2},sorted_design_scores1_mean(2)};
    end
end
    
cd ../

% check to confirm that target region is in the full-length version
for c1 = 1:size(sel_designs,1)
    for c2 = 1:round(size(sel_designs,2)/4)
        curr_seq = dna2rna_B(sel_designs{c1,c2*4-2});
        curr_designs = sel_designs{c1,c2*4-1};
        for c3 = 1:size(curr_designs,1)
            if isempty(findstr(curr_seq,dna2rna_B(curr_designs{c3,2})))
                fprintf('PROBLEM exists.\n');
                return;
            end
        end
    end
end

output_table = {};
target_RNAs = {};
sensor_RNAs = {};
short_target_RNAs = {};

counterx = 1;
for c1 = 1:size(sel_designs,1)
    output_table{end+1,1} = sprintf('%s target:',mRNA_groups{c1});
    output_table{end+1,1} = '';
    for c2 = 1:round(size(sel_designs,2)/4)
        output_table{end+1,1} = sprintf('Target region #%d: Average design score (top six devices) = %0.1f',c2,sel_designs{c1,c2*4});
        output_table{end+1,1} = sprintf('%s RNA sequence',sel_designs{c1,c2*4-3});
        output_table{end,2} = dna2rna_B(sel_designs{c1,c2*4-2});
        target_RNAs(end+1,:) = {sel_designs{c1,c2*4-3},dna2rna_B(sel_designs{c1,c2*4-2})};
        output_table{end+1,1} = '';
        output_table(end+1,1:3) = {'Sensor name','Sensor RNA sequence','Target RNA subsequence'};
        temp_table = sel_designs{c1,c2*4-1};
        subtarget_name = sel_designs{c1,c2*4-3};
        if 1
            temp_table(FINAL_NUM_DESIGNS+1:end,:) = [];
        end
        for c3 = 1:size(temp_table,1)
            sensor_seq = dna2rna_B(concatDNA('GGG',temp_table{c3,1}));
            index = findstr(sensor_seq,linker21);
            sensor_seq = [sensor_seq(1:index-1),linker21];
            output_table(end+1,1:3) = {sprintf('SeriesA_%s_sens%s',subtarget_name,char('A'-1+c3)),sensor_seq,dna2rna_B(temp_table{c3,2})};
            sensor_RNAs(end+1,:) = {sprintf('SeriesA_%s_sens%s',subtarget_name,char('A'-1+c3)),sensor_seq};
            short_target_RNAs(end+1,:) = {sprintf('SeriesA_%s_targ%s',subtarget_name,char('A'-1+c3)),dna2rna_B(concatDNA('',temp_table{c3,2}))};
        end
        output_table{end+1,1} = '';
        counterx = counterx + 1;
    end
end    
output_table(end,:) = [];
mkdir2('final_designs');
cd final_designs;
writeCSVFilefromStructure('SeriesA_final_design_info.csv',output_table);
target_DNAs = {};
for c1 = 1:size(target_RNAs,1)
    curr_seq = rna2dna_B(target_RNAs{c1,2});
    curr_seq = [T7_primer(1:end-3),concatDNA('GGG',curr_seq)];
    curr_seq = d_revcomp(curr_seq);
    target_DNAs(end+1,:) = {sprintf('%s_AS',target_RNAs{c1,1}),curr_seq};
end
sensor_DNAs = {};
for c1 = 1:size(sensor_RNAs,1)
    curr_seq = rna2dna_B(sensor_RNAs{c1,2});
    curr_seq = [T7_primer(1:end-3),concatDNA('GGG',curr_seq)];
    fprintf('%s: %s\n',sensor_RNAs{c1,1},curr_seq);
    sensor_DNAs(end+1,:) = {sensor_RNAs{c1,1},curr_seq};
end
writeCSVFilefromStructure('SeriesA_target_DNA_sequences.csv',target_DNAs);
writeCSVFilefromStructure('SeriesA_sensor_DNA_sequences.csv',sensor_DNAs);
writeCSVFilefromStructure('SeriesA_subtarget_RNA_sequences.csv',short_target_RNAs);
cd ../
