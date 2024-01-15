%% Are diffuse onsets associated with less favorable surgical outcomes?
% Here we will label seizures with at least 75% of regions in onset as
% diffuse

% We will perform chi-square test of those with/without diffuse onsets
% across outcome groups 

%%
save_loc = [save_folder 'diffuse/'];
mkdir(save_loc)
save_loc = [save_loc  char(chan_or_roi) '_'];

for det = det_meths
    % Create a structure diffuse_tab to store binary of diffuse onset
    % per seizure
    if det ~= "clo"
        sz_count = cellfun(@length, final_output.Segment_ids);
        sz_count(~cellfun(@iscell,final_output.Segment_ids)) = 1;
        tab_fill = nan(sum(sz_count),1);
        diffuse_tab.(sprintf(det)) = table(repelem(final_output.Patient_id,sz_count,1),...
            repelem(final_output.outcome_cat,sz_count,1),tab_fill,...
            tab_fill,tab_fill,...
            'VariableNames', ["Subj_id", "Outcome", "sz", "region_count", "is_diffuse"]);
    else
        tab_fill = nan(size(final_output,1),1);
        diffuse_tab.(sprintf(det)) = table(final_output.Patient_id, final_output.outcome_cat,...
            repelem("",size(final_output,1),1), tab_fill, tab_fill,...
            'VariableNames', ["Subj_id", "Outcome", "sz", "region_count", "is_diffuse"]);
    end

    % Iterate through each subject and compute diffusivity measures for
    % each seizure 
    for subj = 1:size(final_output,1)
        subj_onset = final_output(subj,:);
        onset = subj_onset.(sprintf("%s_%s", det, chan_or_roi)){:}; 

        reg_count = nan(size(onset,2),1);
        is_diffuse = nan(size(onset,2),1);

        % Seizure specific process
        seizures = 1:size(onset,2);
        for sz = seizures % Assess one seizure at a time 
            sz_onset = onset(:,sz);
            if sum(sz_onset) == 0
                continue
            end
            
            reg_count(sz) = sum(sz_onset);
            is_diffuse(sz) = mean(sz_onset) >= diffuse_thresh;
        end
        if det == "clo"
            seizures = "clo";
        end
        subj_tab = table(seizures',reg_count, is_diffuse,...
            'VariableNames', ["sz", "region_count", "is_diffuse"]);
        diffuse_tab.(sprintf(det))(diffuse_tab.(sprintf(det)).Subj_id ==...
            string(subj_onset.Patient_id),3:end)= subj_tab;
    end
end

%% Compute subject-level summary of diffusivity across seizures
summ_measures = ["mean", "max"];
tab_fill = nan(length(det_meths)*length(summ_measures),1);
chi2diffuse = table(repelem(det_meths,1,2)',...
    repmat(["has_diffuse", "most_diffuse"],1,length(det_meths))',...
    tab_fill,tab_fill,repmat("",length(tab_fill),1), 'VariableNames',...
    ["det", "comp", "Test-stat", "p", "Test"]);

% Repeat analyses for each detection method 
for det = det_meths
    summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)) = groupsummary(diffuse_tab.(sprintf(det)),...
    ["Subj_id", "Outcome"],summ_measures, "is_diffuse");

    summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)).has_diffuse =...
        summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)).max_is_diffuse == 1;
    summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)).most_diffuse =...
        summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)).mean_is_diffuse >= 0.5;
    
    figure()
    tiledlayout(1,2)
    for comp = ["has_diffuse", "most_diffuse"]
        nexttile
        heatmap(summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)), comp, "Outcome" )
        [cont_table, chi2,p] = crosstab(summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)).(sprintf(comp)),...
            summarised_diffuse_tab.(sprintf(chan_or_roi)).(sprintf(det)).Outcome);
        if any(any(cont_table <5))
            test = "Fisher";
            [~, p] = fishertest(cont_table);
            chi2diffuse(chi2diffuse.det == det & chi2diffuse.comp == comp,3:end) =  table(NaN, p, test);
        else
            test = "Chi-sq";
             chi2diffuse(chi2diffuse.det == det & chi2diffuse.comp == comp,3:end) = table(chi2, p, test); 
        end

    end
    sgtitle(sprintf(det))
end

chi2_tab.(sprintf(chan_or_roi)).is_diffuse = chi2diffuse;


















% %% Are more diffuse onsets associated with less favorable surgical outcomes?
% % Here we will estimate the diffusivity of seizure onsets based on the
% % maximum distance between the centers of regions in onset and onset volume
% 
% %% NATHAN IS WRITING CODE TO GET REGION DISTANCES AND VOLUMES TO REPLACE CURRENT ATLAS DISTANCES 
% 
% %%
% save_loc = [save_folder, 'more_diffuse/'];
% mkdir(save_loc)
% 
% addpath(genpath("../Nathan Code/SimpleRouteExample/"))
% 
% for det = det_meths
%     % Create a structure diffuse_tab to store seizure-wise diffusivity
%     % values 
%     if det ~= "clo"
%         sz_count = cellfun(@length, final_output.Segment_ids);
%         sz_count(~cellfun(@iscell,final_output.Segment_ids)) = 1;
%         tab_fill = nan(sum(sz_count),1);
%         diffuse_tab.(sprintf(det)) = table(repelem(final_output.Patient_id,sz_count,1),...
%             repelem(final_output.outcome_cat,sz_count,1),tab_fill,...
%             tab_fill, tab_fill, tab_fill, tab_fill, tab_fill,...
%             'VariableNames', ["Subj_id", "Outcome", "sz", "region_count", "max_dist", "onset_vol", "prop_max_dist",...
%             "prop_onset_vol"]);
%     else
%         tab_fill = nan(size(final_output,1),1);
%         diffuse_tab.(sprintf(det)) = table(final_output.Patient_id, final_output.outcome_cat,...
%             repelem("",size(final_output,1),1), tab_fill, tab_fill, tab_fill,...
%             tab_fill, tab_fill, 'VariableNames', ["Subj_id", "Outcome",...
%             "sz", "region_count", "max_dist", "onset_vol", "prop_max_dist",...
%             "prop_onset_vol"]);
%     end
%     % Iterate through each subject and compute diffusivity measures for
%     % each seizure 
%     for subj = 1:size(final_output,1)
%         subj_onset = final_output(subj,:);
%         onset = subj_onset.(sprintf("%s_%s", det, chan_or_roi)){:}; 
%         regions = subj_onset.(sprintf("roi_names_%s", extractAfter(chan_or_roi, "roi_"))){:};
% 
%         region_index = find_preserved(names,regions);
%         onset = logical(onset);
% 
%         reg_count = nan(size(onset,2), 1);
%         max_dist = nan(size(onset,2), 1);
%         onset_vol = nan(size(onset,2), 1);
% 
%         subj_max_dist = max(max(dists(region_index,region_index)));
%         subj_max_vol = sum(vols(region_index));
% 
%         % Seizure specific process
%         seizures = 1:size(onset,2);
%         for sz = seizures% Assess one seizure at a time 
%             % Extend onset to include all regions in relevant hemisphere
%             onset_full = zeros(length(names),1);
%             onset_full(region_index(onset(:,sz))) = 1;
% 
%             if sum(onset_full) == 0
%                 continue
%             end
%         
%             reg_count(sz) = sum(onset_full);
%             max_dist(sz) = max(max(dists(onset_full==1,onset_full==1)));
%             onset_vol(sz) = sum(vols(onset_full==1));
%         end
% 
%         if det == "clo"
%             seizures = "clo";
%         end
% 
%         % Replace distances of zero with NaN as they are seizure in just one
%         % region, therefore a distance of zero is mirepresentative
%         max_dist(max_dist == 0) = NaN;
% 
%         subj_tab = table(seizures',reg_count, max_dist, onset_vol, ...
%             max_dist/subj_max_dist, onset_vol/subj_max_vol,...
%             'VariableNames', ["sz", "region_count", "max_dist", "onset_vol",...
%             "prop_max_dist", "prop_onset_vol"]);
% 
%         diffuse_tab.(sprintf(det))(diffuse_tab.(sprintf(det)).Subj_id == string(subj_onset.Patient_id),3:end)=...
%             subj_tab;
%     end
% end
% 
% %% Compute subject-level summary of diffusivity across seizures
% % We will enforce inclusion criteria of at least 5 seizures when computing
% % median
% 
% % Repeat analyses for each detection method 
% for det = det_meths
%     if det == "clo" % we can only compute the maximum (i.e. raw) for CLOs 
%                     % as there is only one onset per subject
%         summ_measures = "max";
%     else 
%         summ_measures = ["median", "max"];
%     end
% 
%     summarised_diffuse_tab.(sprintf(det)) = groupsummary(diffuse_tab.(sprintf(det)),...
%     ["Subj_id", "Outcome"],summ_measures, ["max_dist", "onset_vol", "prop_max_dist", "prop_onset_vol"]);
% 
% %    % Replace values with NaN where inclusion criteria have not been met
% %     for col = find(contains(summarised_diffuse_tab.(sprintf(det)).Properties.VariableNames, "median"))
% %         vals = table2array(summarised_diffuse_tab.(sprintf(det))(:,col));
% %         vals(summarised_diffuse_tab.(sprintf(det)).GroupCount <5) = NaN;
% %         summarised_diffuse_tab.(sprintf(det))(:,col) = table(vals);
% %     end
%     
% end
% %% VISUALISATIONS
% 
% % Create look-up table for comparisons and associated y limits
% y_lim_lookup = table(string(summarised_diffuse_tab.(sprintf(det)).Properties.VariableNames(4:end))',...
%     [repmat([-20,160],2,1); repmat([0,4.5*10^5],2,1); repmat([-0.2,1.2],4,1)],...
%     'VariableNames', ["comp", "y_lim"]);
% 
% % Create violin plots comparing outcome groups
% for det = det_meths
%     comps = (string(summarised_diffuse_tab.(sprintf(det)).Properties.VariableNames(4:end)));
% 
%     for comp_var = ["dist", "vol"]
%         comp_ind = contains(comps, comp_var);
%         % Maximum distance between centres of regions
%         half_violin_upd(summarised_diffuse_tab.(sprintf(det)), comps(comp_ind),...
%             summarised_diffuse_tab.(sprintf(det)).Outcome == out_grps(2),...
%             "y_lim",y_lim_lookup(contains(y_lim_lookup.comp,comps(comp_ind)),:).y_lim,...
%             "save_fig", save_fig, "file_type","svg",...
%             "save_loc",sprintf("%s%s_%s_comp", save_loc, det, comp_var), "grp_names", out_grps)
%         sgtitle(sprintf("%s diffusivity using %s", det, comp_var))
%     end
% 
% end
% 
% %% AUC for separating outcome groups based on diffusivity markers
% for det = det_meths
%     comps = (string(summarised_diffuse_tab.(sprintf(det)).Properties.VariableNames(4:end)));
%     auc_tabs.diffuse.(sprintf(det)) = outcome_auc(summarised_diffuse_tab.(sprintf(det)), comps,...
%         summarised_diffuse_tab.(sprintf(det)).Outcome == out_grps(2), "n_perm",n_perm,"plot_fig",1,"save_fig",save_fig,...
%         "save_loc",save_loc, "file_type", "svg", "main_title", sprintf("AUCs for comparing outcome groups based on diffusivity of %s onsets", det));
% end
