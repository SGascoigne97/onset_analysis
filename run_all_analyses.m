% Run this script for all results and supplementary results
%% Code for results in Incomplete resection of the icEEG seizure onset zone
% is not associated with post-surgical outcomes"
% (S J Gascoigne et al., 2024)

%% Load data and Lausanne atlases
% Set working directory to subquestions folder
addpath(genpath('functions'))
addpath('subquestions/')

load("tables/final_output.mat")
load('roi_info/ATLAS.mat')
final_output = final_output(final_output.outcome ~= 8,:); % ILAE 8 is used 
% to encode unknown outcome - remove such subjects BEFORE completing downstram analysis

n_perm = 1000; % All AUC p-values are computed based on permutation test 
               % with 1000 permutations
out_thresh = 2;%2; % ILAE 1-2 is considered as a 'favourable' outcome whilst 
                % ILAE 3+ is considered as an 'unfavourable' outcome
consensus_thresh = 0.5; % Set threshold for inclusion of regions in 
                        % consensus onset, here we will include regions 
                        % present in at least half of the subject's seizures
most_resec_thresh = 0.5; % Set threshold above which subjects are considered 
                         % as having 'most' of their onset resected

det_meths = ["clo", "imprint", "EI"]; % List onset detection methods you 
                                      % are interested in analysing
                                      

% Set parameters for saving figures
save_fig = 1; % Each figure that is created will be saved 
save_folder = 'figures/paper_figures/Figure 2/'; % Location where subfolders will
                                          % be created to store figures and tables
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end
file_type = "svg"; % Figures will be saved as svgs so they can be pulled 
                   % into illustrator for paper figures

%% Set parameters for analyses
for parc = ["chan", "roi_120", "roi_250"]    
    [chan_or_roi, n_perm, out_thresh, consensus_thresh, most_resec_thresh,...
        det_meths, save_fig, save_folder, file_type, out_grps] =...
        argument_validation(final_output, parc, n_perm, out_thresh, consensus_thresh,...
        most_resec_thresh, det_meths, save_fig, save_folder, file_type); % validate all optional arguments
    
    %% Additional parameters used throughout analyses (not to be adjused)
    atl_inds = 2*(atlas.scale');
    if contains(chan_or_roi, "roi")
        atl = atlas(atl_inds == str2double(extractAfter(chan_or_roi, "_")),:);
        vols = atl.vol{:};
        names = atl.name{:};
        dists = atl.dists{:};
    else
        atl = NaN;
        vols = NaN;
        names = NaN;
        dists = NaN;
    end
    
    %% Add outcome category column to final_output table
    % Outcome category (binary based on selected outcome threshold)
    final_output.outcome_cat = categorical(final_output.outcome>out_thresh,[0,1], out_grps);
    %% 3.1 icEEG Seizure onset regions tend to be resected, but more complete 
    % resections are not associated with more favourable surgical
    %outcomes
    
    % We will compute the proportion of the seizure onset zone (CLO and IOLA)
    % that was subsequently resected then compare across outcome groups
    
    is_onset_resected
    
    %% 3.2 Larger onsets and resections are not associated with surgical outcomes
    
    % Here we will compute the size of seizure onsets (both CLO and
    % automatically captured) using both count of regions and volume (based on
    % controls) and compare across surgical outcome groups
    
    ons_resec_as_prop
        
    %% SUPPLEMENTARY
    %% S2 Subject Metadata
    % Here we will look at various metadata variables to report in the 
    % supplementary results. Further, we will investigate if any metadata
    % groups have significant differences in surgical outcomes
    % have confounding effect on surgical outcomes
    metadata
    
    % Look at tables for results:
    % meta_tab_cont
    % cat_tabs
    
    %% S4.2 ALO is not highly concordant with CLO
    % We will consider concordance between the consensus onset and the maximum
    % consensus between the IOLA onset for any one seizure and the clinically
    % labelled onset (CLO)
    
    concordance_with_clo
    % concord_tab: subject level concordance between consensus onset and CLO
    % and maximum concordance between CLO and any one IOLA onset
    
    % auc_clo_comp: table of AUCs and associated p-values for distinguishing
    % surgical outcome groups based on concordance 

    close all
end



%% Comparing parcellations
clear combined_tab
for chan_or_roi = ["chan", "roi_120", "roi_250"]
    for comp = ["clo_concord", "onset_resec"]
        fields = string(fieldnames(auc_tabs.(sprintf(chan_or_roi)).(sprintf(comp))));
        for field_ind = 1:length(fields)
            field = fields(field_ind);
            value_tab = auc_tabs.(sprintf(chan_or_roi)).(sprintf(comp)).(sprintf(field));
            field_tab = table(chan_or_roi, comp, field);
            if exist("combined_tab","var")
                combined_tab = [combined_tab; [repmat(field_tab,size(value_tab,1),1), value_tab]];
            else
                combined_tab = [repmat(field_tab,size(value_tab,1),1), value_tab];
            end
        end
    end
    for comp = "onset_size"
        comp_tab = auc_tabs.(sprintf(chan_or_roi)).(sprintf(comp));
        field_tab = table(extractBefore(comp_tab.Comparison, "_"), extractAfter(comp_tab.Comparison, "_"), 'VariableNames', ["field", "Comparison"]);
        field_tab = [repmat(table(chan_or_roi, comp),size(field_tab,1),1), field_tab];
        % Add to combined table
        combined_tab = [combined_tab; field_tab, comp_tab(:,2:3)];
        
    end
end

save_folder = "figures/paper_figures/new_results";

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare results across parcellation schemes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for comp = ["clo_concord", "onset_resec", "onset_size"]
    comp_tab = combined_tab(combined_tab.comp == comp,:);
end

%% Proportion of onset resected
parcs = ["chan", "roi_120", "roi_250"];
comp_parc_tab.onset_resec = array2table(nan(0,7), 'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "t", "DoF", "p"]);

% Is onset resected (Section 3.2)
for det = det_meths
    resec_tab= table(repelem(parcs(1:length(parcs)-1), 1, length(parcs)*2)', ...
    repelem(repmat(parcs, 1, (length(parcs)-1)),1,2)', repmat(["Perc", "Perc_vol"], 1, (length(parcs)-1)*length(parcs))',...
    repmat(det, 2*(length(parcs)-1)*length(parcs), 1 ),...
    nan(2*(length(parcs)-1)*length(parcs),1), nan(2*(length(parcs)-1)*length(parcs),1), ...
    nan(2*(length(parcs)-1)*length(parcs),1),...
    'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "t", "DoF", "p"]);
    for parc1 = 1:length(parcs)
        for parc2 = 1:length(parcs)
            if parc1 == parc2 | parc1 > parc2
                continue
            else
                tab1 = comp_with_resec.(sprintf(parcs(parc1))).(sprintf(det));
                tab2 = comp_with_resec.(sprintf(parcs(parc2))).(sprintf(det));

                % Join tables together to match patients 
                across_parc_tab = innerjoin(tab1, tab2, 'Keys', ...
                    ["Patient_id", "Outcome", "Sz_count"]);

                if contains(parcs(parc1), "roi") & contains(parcs(parc1), "roi")
                    comps = ["Perc", "Perc_vol"];
                else
                    comps = "Perc";
                end

                for comp = comps
                    diffs = across_parc_tab.(sprintf("%s_tab1", comp)) - across_parc_tab.(sprintf("%s_tab2", comp));
                    figure()
                    sgtitle(strrep(sprintf("Comparing %s %s across %s and %s", det, comp, parcs(parc1), parcs(parc2) ), "_", " "))
                    subplot(121)
                    scatter(across_parc_tab.(sprintf("%s_tab1", comp)),...
                        across_parc_tab.(sprintf("%s_tab2", comp)), 'filled',...
                        'XJitter', 'randn', 'XJitterWidth', 0.05,...
                        'YJitter', 'randn', 'YJitterWidth', 0.05)
                    lsline()
                    xlabel(parcs(parc1))
                    ylabel(parcs(parc2))
                    subplot(122)
                    [~,p,~, stats] = ttest(diffs);
                    histogram(diffs)
                    title(sprintf("Histogram of differences \n t(%d) = %.3f, p = %.3f", length(diffs)-1, stats.tstat, p))

                    % Store test statistic, degrees of freedom, and
                    % p-values in table
                    resec_tab(resec_tab.Scheme_1 == parcs(parc1) &...
                        resec_tab.Scheme_2 == parcs(parc2) & ...
                        resec_tab.comp == comp & resec_tab.det == det, 5:7) = ...
                        table(stats.tstat, length(diffs)-1, p);
                end
            end
        end
    end

    resec_tab = resec_tab(~isnan(resec_tab.t),:);
    comp_parc_tab.onset_resec = [comp_parc_tab.onset_resec; resec_tab];
end


%% Onset size
save_loc = sprintf("%s/onset_size", save_folder);
comp_parc_tab.onset_size = array2table(nan(0,7), 'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "rho", "DoF", "p"]);

% Is onset resected (Section 3.1)
for det = det_meths
    resec_tab = table(repelem(parcs(1:length(parcs)-1), 1, length(parcs)*3)', ...
        repelem(repmat(parcs, 1, (length(parcs)-1)),1,3)', ...
        repmat(["count", "prop", "vol"], 1, (length(parcs)-1)*length(parcs))',...
        repmat(det, 3*(length(parcs)-1)*length(parcs), 1 ),...
        nan(3*(length(parcs)-1)*length(parcs),1),...
        nan(3*(length(parcs)-1)*length(parcs),1), ...
        nan(3*(length(parcs)-1)*length(parcs),1),...
        'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "rho", "DoF", "p"]);
    for parc1 = 1:length(parcs)
        for parc2 = 1:length(parcs)
            if parc1 == parc2 | parc1 > parc2
                continue
            else
                if contains(parcs(parc1), "roi") & contains(parcs(parc1), "roi")
                    comps = ["count", "prop", "vol"];
                else
                    comps = ["count", "prop"];
                end

                tab1 = onset_parc_tab.(sprintf(parcs(parc1)));
                tab2 = onset_parc_tab.(sprintf(parcs(parc2)));

                % Join tables together to match patients 
                across_parc_tab = innerjoin(tab1, tab2, 'Keys', ...
                    ["Patient_id", "Outcome"]);

                for comp = comps
                    grp_1_vals = across_parc_tab.(sprintf("%s_%s_tab1", det, comp));
                    grp_2_vals = across_parc_tab.(sprintf("%s_%s_tab2", det, comp));

                    rm_obs = isnan(grp_1_vals)| isnan(grp_2_vals);
                    grp_1_vals = grp_1_vals(~rm_obs);
                    grp_2_vals = grp_2_vals(~rm_obs);

                    f = figure();
                    f.Position = [10,10,700,900];
                    sgtitle(strrep(sprintf("Comparing %s %s across %s and %s", det, comp, parcs(parc1), parcs(parc2) ), "_", " "))
%                     subplot(121)
                    scatter(grp_1_vals, grp_2_vals,'filled',...
                        'XJitter', 'randn', 'XJitterWidth', range(grp_1_vals)/20,...
                        'YJitter', 'randn', 'YJitterWidth', range(grp_2_vals)/20)
                    lsline()    
                    axis square
                    xlabel(parcs(parc1))
                    ylabel(parcs(parc2))

                    [r, p] = corr(grp_1_vals, grp_2_vals, "type","Spearman");

                    title(sprintf("Spearman's rank \n r(%d) = %.3f, p = %.3f", length(grp_1_vals)-2, r, p))
                    saveas(f, sprintf('%s/onset_size_%s_%s(%s_v_%s).svg', save_loc, det, comp, parcs(parc1), parcs(parc2)))

                    % Store test statistic, degrees of freedom, and
                    % p-values in table
                    resec_tab(resec_tab.Scheme_1 == parcs(parc1) &...
                        resec_tab.Scheme_2 == parcs(parc2) & ...
                        resec_tab.comp == comp & resec_tab.det == det, 5:7) = ...
                        table(r, length(grp_1_vals)-2, p);
                end
            end
        end
    end

    resec_tab = resec_tab(~isnan(resec_tab.rho),:);
    comp_parc_tab.onset_size = [comp_parc_tab.onset_size; resec_tab];
end

%% Compare AUCs across parcellation schemes
parcs = ["chan", "roi_120", "roi_250"];
det_meths = ["clo", "imprint", "EI"]; 
comps = ["onset_resec", "onset_size", "clo_concord"]; % We will look at "resection_size" separately for now 

for comp = comps
    auc_comp_mat = nan(length(det_meths), length(parcs));
    p_comp_mat = auc_comp_mat;

    for det_meth = det_meths
        for parc = parcs
            row = det_meths == det_meth;
            col = parcs == parc;
            
            if comp == "onset_resec" 
                auc_and_p = auc_tabs.(sprintf(parc)).(sprintf(comp)).(sprintf(det_meth));
                auc_and_p = auc_and_p(~contains(auc_and_p.Comparison, "_vol"),:);
                auc_comp_mat(row,col) = auc_and_p.AUC;
                p_comp_mat(row,col) = auc_and_p.p;

            elseif comp == "clo_concord"
                if det_meth == "clo"
                    continue
                end
                auc_and_p = auc_tabs.(sprintf(parc)).(sprintf(comp)).(sprintf(det_meth));
                auc_and_p_consensus = auc_and_p(contains(auc_and_p.Comparison, "consensus"),:);
                auc_comp_mat_cons(row,col) = auc_and_p_consensus.AUC;
                p_comp_mat_cons(row,col) = auc_and_p_consensus.p;

                auc_and_p_max = auc_and_p(contains(auc_and_p.Comparison, "max_clo_concord"),:);
                auc_comp_mat_max(row,col) = auc_and_p_max.AUC;
                p_comp_mat_max(row,col) = auc_and_p_max.p;

            elseif comp == "onset_size" 
                % We will split the results into onset size and onset as a 
                % proportion of the implantation 
                auc_and_p = auc_tabs.(sprintf(parc)).(sprintf(comp));
                auc_and_p_count = auc_and_p(contains(auc_and_p.Comparison, sprintf("%s_count", det_meth)),:);
                auc_comp_mat_count(row,col) = auc_and_p_count.AUC;
                p_comp_mat_count(row,col) = auc_and_p_count.p;

                auc_and_p_prop = auc_and_p(contains(auc_and_p.Comparison, sprintf("%s_prop", det_meth)),:);
                auc_comp_mat_prop(row,col) = auc_and_p_prop.AUC;
                p_comp_mat_prop(row,col) = auc_and_p_prop.p;
            end
        end
    end
    if comp == "onset_resec" 
        auc_comp_tab.(sprintf(comp)) = array2table(auc_comp_mat,'RowNames', det_meths , 'VariableNames', parcs);
        p_comp_tab.(sprintf(comp)) = array2table(p_comp_mat,'RowNames', det_meths , 'VariableNames', parcs);
    elseif comp == "onset_size"
        auc_comp_tab.(sprintf("%s_count", comp)) = array2table(auc_comp_mat_count,'RowNames', det_meths, 'VariableNames', parcs);
        p_comp_tab.(sprintf("%s_count", comp)) = array2table(p_comp_mat_count,'RowNames', det_meths , 'VariableNames', parcs);
        auc_comp_tab.(sprintf("%s_prop", comp)) = array2table(auc_comp_mat_prop,'RowNames', det_meths, 'VariableNames', parcs);
        p_comp_tab.(sprintf("%s_prop", comp)) = array2table(p_comp_mat_prop,'RowNames', det_meths, 'VariableNames', parcs);
    elseif comp == "clo_concord"
        auc_comp_tab.(sprintf("%s_consensus", comp)) = array2table(auc_comp_mat_cons(2:end,:),'RowNames', det_meths(2:end), 'VariableNames', parcs);
        p_comp_tab.(sprintf("%s_consensus", comp)) = array2table(p_comp_mat_cons(2:end,:),'RowNames', det_meths(2:end) , 'VariableNames', parcs);
        auc_comp_tab.(sprintf("%s_max", comp)) = array2table(auc_comp_mat_max(2:end,:),'RowNames', det_meths(2:end), 'VariableNames', parcs);
        p_comp_tab.(sprintf("%s_max", comp)) = array2table(p_comp_mat_max(2:end,:),'RowNames', det_meths(2:end), 'VariableNames', parcs);
    else
        auc_comp_tab.(sprintf(comp)) = array2table(auc_comp_mat(1,:), 'VariableNames', parcs);
        p_comp_tab.(sprintf(comp)) = array2table(p_comp_mat(1,:), 'VariableNames', parcs);
    end
end
%% Create heatmaps (for supplementary)
comps = ["onset_resec", "onset_size_count", "onset_size_prop", "clo_concord_consensus", "clo_concord_max"];

f = figure();
f.Position = [10,10,1000,1000];
cm = [[163 2 52]/255  % red
    [0 118 192]/255 % blue
    [163 2 52]/255];% red
cm = interp1(cm, 1:0.01:size(cm,1));

tiledlayout(5,1)

for comp = comps
    nexttile
    h_auc = heatmap(table2array(auc_comp_tab.(sprintf(comp))));
    clim([0,1])
    title(comp)
   % colorbar off

    if comp ~= "resection_size"
        h_auc.YDisplayLabels = auc_comp_tab.(sprintf(comp)).Properties.RowNames;
    end
        h_auc.XDisplayLabels = auc_comp_tab.(sprintf(comp)).Properties.VariableNames;
        colormap(flipud(cm))
        clim([0,1])
        %colorbar off
         h_auc.CellLabelFormat = '%.3f';
        

end
sgtitle("AUCs across parcellations and onset detection methods")
saveas(f, sprintf("%s/comparing_parcellations_auc.svg", save_folder))
%%
% Plot p-values
f = figure();
f.Position = [10,10,1000,1000];
cm = [1,0,0;1,1,1];
cm = interp1(cm, 1:0.01:size(cm,1));
tiledlayout(5,1)
for comp = comps
    nexttile
    h_p = heatmap(table2array(p_comp_tab.(sprintf(comp))));
    title(comp)

    if comp ~= "resection_size"
        h_p.YDisplayLabels = p_comp_tab.(sprintf(comp)).Properties.RowNames;
    end
        h_p.XDisplayLabels = p_comp_tab.(sprintf(comp)).Properties.VariableNames;
        colormap(cm)
        clim([0,0.05])
        %colorbar off
        h_p.CellLabelFormat = '%.3f';

end
sgtitle("p-values across parcellations and onset detection methods")
saveas(f, sprintf("%s/comparing_parcellations_p.svg", save_folder))

%% Scatterplots across parcellations (not saved)
% figure()
% tiledlayout(7,3)
% col_names = string(onset_parc_tab.roi_120.Properties.VariableNames(3:end));
% col_names = col_names(~contains(col_names, "vol")); % We have excluded volume from our analyses
% 
% parcs = ["chan", "roi_120", "roi_250"];
% 
% for comp = col_names
%     for parc_ind_1 = 1:length(parcs)
%         for parc_ind_2 = 1:length(parcs)
%             if parc_ind_1 >= parc_ind_2
%                 continue
%             end
%             parc1 = parcs(parc_ind_1);
%             parc2 = parcs(parc_ind_2);
%             nexttile
%             gscatter(onset_parc_tab.(sprintf(parc1)).(sprintf(comp)),...
%                 onset_parc_tab.(sprintf(parc2)).(sprintf(comp)), onset_parc_tab.(sprintf(parc1)).Outcome ,...
%                 'filled')
%             lsline()
%             title(sprintf("Comparing count of %s regions", extractBefore(comp,"_")))
%             xlabel(strrep(sprintf(parc1), "_", " "))
%             ylabel(strrep(sprintf(parc2), "_", " "))
%         end
%     end
% end

