

dots = char(zeros(1,5000)+'.');
bleft = char(zeros(1,5000)+'(');
bright = char(zeros(1,5000)+')');
Ns = char(zeros(1,5000)+'N');
%H56 hairpin sequence
H00 = 'GGGUGAAUGAAUUGUAGGCUUGUUAUAGUUAUGAACAGAGGAGACAUAACAUGAACAAGCCUAACCUGGCGGCAGCGCAAAAGAUGCGUAAAGGAGAAGAACUUUUCACUGG';
%     '..............((((((((((...((((((...........))))))...)))))))))).(((.......(((((.....)))))..))).((((....)))).....'


top_part_1 = H00(3+12+9+1:3+12+9+9+11+9);
top_part_2 = H00(3+12+9+1-3:3+12+9+9+11+9+3);

%new designs: use 30-nt toehold to max that out (ideally), and use same
%unwinding of 6-bp

hairpin_suffix_addon_len = 29;
mRNA_names0 = input_table(:,1);
mRNA_sequences0 = input_table(:,2);
mRNA_inner_sequences0 = input_table(:,3);
T_set0 = input_table(:,4);
output_name0 = input_table(:,5);
output_seq0 = input_table(:,6);


if RUN_ANTISENSE
    for c1 = 1:length(mRNA_names0)
        mRNA_names0(end+1) = strcat(mRNA_names0(c1), '_antisense');
        mRNA_sequences0(end+1) = {seqrcomplement(mRNA_sequences0{c1})};
        mRNA_inner_sequences0(end+1) = {seqrcomplement(mRNA_inner_sequences0{c1})};
        T_set0(end+1) = T_set0(c1);
        output_name0(end+1) = output_name0(c1);
        output_seq0(end+1) = output_seq0(c1);
    end
end


mRNA_names = {};
mRNA_sequences = {};
mRNA_inner_sequences = {};
mRNA_temps = [];
mRNA_hairpin_suffix = {};
for c1 = 1:length(mRNA_names0)
    T_set = eval(T_set0{c1});
    curr_output_name = output_name0{c1};
    curr_output_seq = dna2rna_B(output_seq0{c1}(1:min(end,hairpin_suffix_addon_len)));
    for c2 = 1:length(T_set)
        mRNA_names{end+1,1} = sprintf('%s_%s_T%d',mRNA_names0{c1},curr_output_name,T_set(c2));
        mRNA_sequences{end+1,1} = mRNA_sequences0{c1};
        mRNA_inner_sequences{end+1,1} = mRNA_inner_sequences0{c1};
        mRNA_temps(end+1,1) = T_set(c2);
        mRNA_hairpin_suffix{end+1,1} = [hairpin_suffix0,curr_output_seq];
    end
end


d_dom = 'AAC';
len_a = 3;
len_b = 3;
len_c = 30;
len_d = length(d_dom);
counter = 1;
target_len = sum([len_a,len_b,len_c]);    

design_type = 'SeriesA';
mkdir2(design_type);
writeCSVFilefromStructure([design_type,'/ZZZ_original_sequence_info.csv'],[mRNA_names,mRNA_sequences,mRNA_inner_sequences]);

cd(design_type);
for c1 = 1:length(mRNA_names)
    [~,~,~] = mkdir(mRNA_names{c1});
    cd(mRNA_names{c1});
    mRNA_seq = mRNA_sequences{c1};
    mRNA_inner_seq = mRNA_inner_sequences{c1};
    hairpin_suffix = mRNA_hairpin_suffix{c1};
    T = mRNA_temps(c1);

    % Defines the secondary structue of the toehold switch
    toehold_len = len_c;
    hairpin_struc = [dots(1:toehold_len),bleft(1:9),dots(1:3),bleft(1:6),dots(1:11),bright(1:6),dots(1:3),bright(1:9),dots(1:len_d+len_a+3),dots(1:length(hairpin_suffix))];
    hairpin_struc_min = hairpin_struc;
    hairpin_struc_min(end-length(hairpin_suffix)+1:end) = 'N';
    toehold_struc_min = [];
    toehold_struc_min(1:length(hairpin_struc_min)) = 'N';
    toehold_struc_min(1:toehold_len) = '.';
    matrix_hairpin_struc_min = convertStrucString2PairProbTable(hairpin_struc_min);
    matrix_toehold_struc_min = convertStrucString2PairProbTable(toehold_struc_min);
    len_hairpin_struc_min = sum(sum(matrix_hairpin_struc_min));
    target_struc =  [dots(1:target_len)];
    index_after_duplex = target_len + 1;
    
    mRNA_Ns = Ns(1:length(mRNA_seq));
    pp = computePairProbTableRNAT((mRNA_seq),T);
    output_table = cell(length(mRNA_seq)-target_len+1+1,8+3+2);
    output_table(1,:) = {'Position','Probe Defect','Probe Min. Defect','Toehold Defect','Target Defect','mRNA Defect','Active mRNA Defect','DeltaG Active mRNA','Probe/mRNA Binding Defect','Probe/mRNA Active Translation Defect','Binding Fraction','Probe Sequence','Target Sequence'};
    start_point = findstr(mRNA_inner_seq,mRNA_seq);
    end_point = start_point + length(mRNA_inner_seq) - 1;
    if IS_PARALLEL
        parfor c2 = start_point:end_point-target_len+1
            fprintf('%s: Position %d of %d...\n',mRNA_names{c1},c2,length(mRNA_seq)-target_len+1);
            target_seq = mRNA_seq(c2:c2+target_len-1);
            target_seq_comp = rna_comp_strand_string(target_seq);
            toehold = target_seq_comp(1:toehold_len);
            a_dom = target_seq(1:len_a);
            b_dom = target_seq(len_a+1:len_a+len_b);
            a_dom_star = r_revcomp(a_dom);
            b_dom_star = r_revcomp(b_dom);
            hairpin_seq = ['',toehold,b_dom_star,a_dom_star,top_part_2,a_dom,b_dom,d_dom,a_dom_star,r_revcomp(top_part_2(end-2:end)),hairpin_suffix];
            mRNA_struc = mRNA_Ns;
            mRNA_struc(c2:c2+target_len-1) = '.';
            mRNA_pp_struc = convertStrucString2PairProbTable(mRNA_struc);
            hairpin_defect = checkDefectArbT(hairpin_seq,hairpin_struc,T);
            target_defect = checkDefectArbT(target_seq,target_struc,T);
            mRNA_defect = 1-sum(sum(pp.*mRNA_pp_struc))./target_len;
            pp_hpin = computePairProbTableRNAT(hairpin_seq,T);
            active_mRNA_defect = checkUnpairedT(hairpin_seq(index_after_duplex:end),T); 
            [active_mRNA_deltaG,rna_struc] = computeRNAdeltaGT(hairpin_seq(index_after_duplex:end),T);
            hairpin_defect_min = 1-sum(sum(matrix_hairpin_struc_min.*pp_hpin))./len_hairpin_struc_min;
            toehold_defect_min = 1-sum(sum(matrix_toehold_struc_min.*pp_hpin))./(toehold_len);
            [conc1,conc2,conc_dimer,pair_structure,deltaG,pair_prob_table_set] = computeRNApairStructureDeltaGPairsT(hairpin_seq,mRNA_seq,T);

            probe_mRNA_bind_struc = DUnotation2DotBracket(sprintf('D%d (N%d   N%d) N%d',...
                target_len,length(hairpin_seq)-target_len,...
                c2-1,...
                length(mRNA_seq)-c2-target_len+1));
            probe_mRNA_bind_pp_struc = convertStrucString2PairProbTable(probe_mRNA_bind_struc);
            probe_mRNA_bind_len = sum(sum(probe_mRNA_bind_pp_struc));
            probe_mRNA_active_struc = DUnotation2DotBracket(sprintf('N%d U%d N%d',target_len,length(hairpin_seq) - target_len,length(mRNA_seq)));
            probe_mRNA_active_pp_struc = convertStrucString2PairProbTable(probe_mRNA_active_struc);
            probe_mRNA_active_len = sum(sum(probe_mRNA_active_pp_struc));
            probe_mRNA_pp = pair_prob_table_set{4};
            probe_mRNA_bind_defect   = 1 - sum(sum(probe_mRNA_pp.*probe_mRNA_bind_pp_struc))./probe_mRNA_bind_len;
            probe_mRNA_active_defect = 1 - sum(sum(probe_mRNA_pp.*probe_mRNA_active_pp_struc))./probe_mRNA_active_len;

            output_table(c2+1,:) = {num2str(c2),num2str(hairpin_defect),num2str(hairpin_defect_min),num2str(toehold_defect_min),num2str(target_defect),num2str(mRNA_defect),num2str(active_mRNA_defect),num2str(active_mRNA_deltaG),num2str(probe_mRNA_bind_defect),num2str(probe_mRNA_active_defect),num2str(conc_dimer),hairpin_seq,target_seq};
            counter = counter + 1;
        end
    else
        for c2 = start_point:end_point-target_len+1
            fprintf('%s: Position %d of %d...\n',mRNA_names{c1},c2,length(mRNA_seq)-target_len+1);
            target_seq = mRNA_seq(c2:c2+target_len-1);
            target_seq_comp = rna_comp_strand_string(target_seq);
            toehold = target_seq_comp(1:toehold_len);
            a_dom = target_seq(1:len_a);
            b_dom = target_seq(len_a+1:len_a+len_b);
            a_dom_star = r_revcomp(a_dom);
            b_dom_star = r_revcomp(b_dom);
            hairpin_seq = ['',toehold,b_dom_star,a_dom_star,top_part_2,a_dom,b_dom,d_dom,a_dom_star,r_revcomp(top_part_2(end-2:end)),hairpin_suffix];
            mRNA_struc = mRNA_Ns;
            mRNA_struc(c2:c2+target_len-1) = '.';
            mRNA_pp_struc = convertStrucString2PairProbTable(mRNA_struc);
            hairpin_defect = checkDefectArbT(hairpin_seq,hairpin_struc,T);
            target_defect = checkDefectArbT(target_seq,target_struc,T);
            mRNA_defect = 1-sum(sum(pp.*mRNA_pp_struc))./target_len;
            pp_hpin = computePairProbTableRNAT(hairpin_seq,T);
            active_mRNA_defect = checkUnpairedT(hairpin_seq(index_after_duplex:end),T); 
            [active_mRNA_deltaG,rna_struc] = computeRNAdeltaGT(hairpin_seq(index_after_duplex:end),T);
            hairpin_defect_min = 1-sum(sum(matrix_hairpin_struc_min.*pp_hpin))./len_hairpin_struc_min;
            toehold_defect_min = 1-sum(sum(matrix_toehold_struc_min.*pp_hpin))./(toehold_len);
            [conc1,conc2,conc_dimer,pair_structure,deltaG,pair_prob_table_set] = computeRNApairStructureDeltaGPairsT(hairpin_seq,mRNA_seq,T);

            probe_mRNA_bind_struc = DUnotation2DotBracket(sprintf('D%d (N%d   N%d) N%d',...
                target_len,length(hairpin_seq)-target_len,...
                c2-1,...
                length(mRNA_seq)-c2-target_len+1));
            probe_mRNA_bind_pp_struc = convertStrucString2PairProbTable(probe_mRNA_bind_struc);
            probe_mRNA_bind_len = sum(sum(probe_mRNA_bind_pp_struc));
            probe_mRNA_active_struc = DUnotation2DotBracket(sprintf('N%d U%d N%d',target_len,length(hairpin_seq) - target_len,length(mRNA_seq)));
            probe_mRNA_active_pp_struc = convertStrucString2PairProbTable(probe_mRNA_active_struc);
            probe_mRNA_active_len = sum(sum(probe_mRNA_active_pp_struc));
            probe_mRNA_pp = pair_prob_table_set{4};
            probe_mRNA_bind_defect   = 1 - sum(sum(probe_mRNA_pp.*probe_mRNA_bind_pp_struc))./probe_mRNA_bind_len;
            probe_mRNA_active_defect = 1 - sum(sum(probe_mRNA_pp.*probe_mRNA_active_pp_struc))./probe_mRNA_active_len;

            output_table(c2+1,:) = {num2str(c2),num2str(hairpin_defect),num2str(hairpin_defect_min),num2str(toehold_defect_min),num2str(target_defect),num2str(mRNA_defect),num2str(active_mRNA_defect),num2str(active_mRNA_deltaG),num2str(probe_mRNA_bind_defect),num2str(probe_mRNA_active_defect),num2str(conc_dimer),hairpin_seq,target_seq};
            counter = counter + 1;
        end
    end
    output_table(2:start_point,:) = [];
    flag_index = 0;
    for c2 = 1:size(output_table,1)
        if isempty(output_table{c2,1})
            flag_index = c2;
            break;
        end
    end
    if flag_index
        output_table(flag_index:end,:) = [];
    end
    writeCSVFilefromStructure(sprintf('mRNA_RR_results_complete.csv'),output_table);
    counter = 1;
    cd ..
end
cd ..

