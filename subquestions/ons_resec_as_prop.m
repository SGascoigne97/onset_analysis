%% 3.3 Larger onsets are not associated with surgical outcomes 

% Here we will compute the size of seizure onsets (both CLO and
% automatically captured) using both count of regions and volume (based on
% controls) and compare across surgical outcome groups

save_loc = [save_folder 'onset_size/'];

if ~exist(save_loc, "dir")
    mkdir(save_loc)
end

save_loc = [save_loc  char(chan_or_roi) '_'];

tab_fill = nan(size(final_output,1),1);
% Create a table to store output (just counts for channels, counts and
% volumes for rois)
if chan_or_roi == "chan"
    onset_size_tab = table(final_output.Patient_id, final_output.outcome_cat,...
        tab_fill, tab_fill, tab_fill,tab_fill,tab_fill,tab_fill, 'VariableNames',...
        ["Patient_id", "Outcome", det_meths + "_count", det_meths + "_prop"]);
else
    onset_size_tab = table(final_output.Patient_id, final_output.outcome_cat,...
        tab_fill, tab_fill, tab_fill, tab_fill, tab_fill, tab_fill,tab_fill,tab_fill,tab_fill,...
        'VariableNames', ["Patient_id", "Outcome", ...
        repelem(det_meths,3)+ repmat(["_count", "_vol", "_prop"],1,3)]);
end

impl_size = cell2mat(cellfun(@size, final_output.(sprintf("imprint_%s", parc)), 'UniformOutput',false));
onset_size_tab.impl_size = impl_size(:,1);
     

for subj = 1:size(final_output)
    subj_onset = final_output(subj,:);
   
    for meth = det_meths
        onset_mat = subj_onset.(sprintf("%s_%s", meth, chan_or_roi)){:};
%         onset_cons = mean(onset_mat,2) >= consensus_thresh;
        onset_cons = mean(onset_mat(:,sum(onset_mat,1)>0),2)>=consensus_thresh;
        if sum(onset_cons) >0
            onset_size_tab(subj,:).(sprintf("%s_count", meth)) = sum(onset_cons);
            onset_size_tab(subj,:).(sprintf("%s_prop", meth)) = sum(onset_cons)/onset_size_tab(subj,:).impl_size;
        end

        if contains(chan_or_roi, "roi")
            regions = subj_onset.(sprintf("roi_names_%s", extractAfter(chan_or_roi, "roi_"))){:};
            onset_regions = regions(onset_cons);
            onset_regions_in_atlas = contains(names, onset_regions);
            if sum(vols(onset_regions_in_atlas))>0
                onset_size_tab(subj,:).(sprintf("%s_vol", meth))  = sum(vols(onset_regions_in_atlas));
            end   
        end
    end
end

onset_parc_tab.(sprintf(chan_or_roi)) = onset_size_tab;

%% VISUALISATIONS
% Creating half-violin plots comparing onset sizes between outcome groups
comp_vars = string(onset_size_tab.Properties.VariableNames);
comp_vars = comp_vars(3:end);
if chan_or_roi == "chan"
    comp_vars = [comp_vars(contains(comp_vars, "count")), comp_vars(contains(comp_vars, "prop"))];
     half_violin(onset_size_tab, comp_vars(1:3),...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [1,1,1],"save_fig",save_fig,...
        "y_lim", repelem([0,20],3,1) ,...
        "file_type","svg", "save_loc", sprintf("%sviolins_count", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
    half_violin(onset_size_tab, comp_vars(4:6),...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [0.1,0.1,0.1],"save_fig",save_fig,...
        "y_lim", repelem([0,1],3,1) ,...
        "file_type","svg", "save_loc", sprintf("%sviolins_prop", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
else
    comp_vars = [comp_vars(contains(comp_vars, "count")),...
        comp_vars(contains(comp_vars, "prop")),comp_vars(contains(comp_vars, "vol"))];
     half_violin(onset_size_tab, comp_vars(1:3),...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [1,1,1],"save_fig",save_fig,...
        "y_lim", repelem([0,20],3,1),...
        "file_type","svg", "save_loc", sprintf("%sviolins_count", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
     half_violin(onset_size_tab, comp_vars(4:6),...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [0.1,0.1,0.1],"save_fig",save_fig,...
        "y_lim", repelem([0,1],3,1) ,...
        "file_type","svg", "save_loc", sprintf("%sviolins_prop", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
    half_violin(onset_size_tab, comp_vars(7:9),...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [NaN, NaN, NaN],"save_fig",save_fig,...
        "y_lim", repelem([0,4*10^5],3,1),...
        "file_type","svg", "save_loc", sprintf("%sviolins_vol", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
end


%% AUC for separating outcome groups based on onset size
auc_tabs.(sprintf(chan_or_roi)).onset_size = outcome_auc(onset_size_tab, comp_vars,...
    onset_size_tab.Outcome == out_grps(2), "n_perm",n_perm,"plot_fig",1,"save_fig",save_fig,...
    "save_loc",save_loc, "file_type", "svg", "main_title", "AUCs comparing outcome groups based on onset size");

%% Group medians to report in paper
comps =  string(onset_size_tab.Properties.VariableNames);
med_tab = table(comps(3:end)', nan(length(comps)-2,1), nan(length(comps)-2,1),...
    nan(length(comps)-2,1), 'VariableNames', ["Comparison", "median" + ["", out_grps]]);
for comp = comps(3:end)
    med = nan(1,3);
    med(1) = median(onset_size_tab.(sprintf(comp)), 'omitnan');
    for grp = 1:2
        med(grp+1) = median(onset_size_tab(onset_size_tab.Outcome == out_grps(grp),:).(sprintf(comp)), 'omitnan');
    end
    med_tab(med_tab.Comparison == comp,2:end) = array2table(med);
end

median_tabs.(sprintf(chan_or_roi)).onset_size = med_tab;