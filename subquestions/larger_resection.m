%% Does a larger resection tend to result in more favourable outcomes?
% Is a larger resection (i.e., more regions resected) associated with better
% post-surgical outcomes?

save_loc = [save_folder 'resection_size/'];
mkdir(save_loc)
save_loc = [save_loc  char(chan_or_roi) '_'];


tab_fill = nan(size(final_output,1),1);
% Create a table to store output (just counts for channels, counts and
% volumes for rois)
if chan_or_roi == "chan"
    resec_size_tab = table(final_output.Patient_id, final_output.outcome_cat,...
        tab_fill, 'VariableNames', ["Patient_id", "Outcome", "resected_count"]);
else
    resec_size_tab = table(final_output.Patient_id, final_output.outcome_cat,...
        tab_fill, tab_fill, 'VariableNames', ["Patient_id", "Outcome", ...
        repelem("resected",2)+ ["_count", "_vol"]]);
end

for subj = 1:size(final_output)
    subj_onset = final_output(subj,:);
  
    resec = subj_onset.(sprintf("resected_%s", chan_or_roi)){:};
    if sum(resec) > 0
        resec_size_tab(subj,:).resected_count = sum(resec);
    end

    if contains(chan_or_roi, "roi")
        regions = subj_onset.(sprintf("roi_names_%s", extractAfter(chan_or_roi, "roi_"))){:};
        resec_regions = regions(logical(resec));
        resec_regions_in_atlas = contains(names, resec_regions);
        if sum(vols(resec_regions_in_atlas))>0
            resec_size_tab(subj,:).resected_vol  = sum(vols(resec_regions_in_atlas));
        end   
    end
end

resec_size_tab = resec_size_tab(~isnan(resec_size_tab.resected_count),:);
resec_parc_tab.(sprintf(chan_or_roi)) = resec_size_tab;

%% VISUALISATIONS
% Creating half-violin plots comparing resection sizes between outcome groups
comp_vars = string(resec_size_tab.Properties.VariableNames);

if chan_or_roi == "chan"
    comp_vars = comp_vars(contains(comp_vars, "count"));
    half_violin(resec_size_tab, comp_vars,...
        double(resec_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", 1,"save_fig",save_fig, "y_lim", [0,14],...
        "file_type","svg", "save_loc", sprintf("%sviolins", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
else
    comp_vars = [comp_vars(contains(comp_vars, "count")), comp_vars(contains(comp_vars, "vol"))];
     half_violin(resec_size_tab, comp_vars,...
        double(resec_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [1,NaN],"save_fig",save_fig,...
        "y_lim", [0,14; 0, 3*10^5],...
        "file_type","svg", "save_loc", sprintf("%sviolins", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
end

%% T-tests to compare mean onset size across outcome groups
t_resec_comp = table(comp_vars', nan(length(comp_vars),1), nan(length(comp_vars),1),...
    nan(length(comp_vars),1),'VariableNames',["Comparison", "t", "df", "p"]);
for comp = comp_vars
    fig = figure("Position",[10,10,500,900]);
    subplot(3,2,1:2)
    boxchart(resec_size_tab.Outcome, resec_size_tab.(sprintf(comp)), 'MarkerStyle','none')
    hold on
    swarmchart(resec_size_tab.Outcome, resec_size_tab.(sprintf(comp)), 'filled')
    hold off
    title(sprintf("%s across outcome groups", comp))
    subplot(3,2,3:4)
    boxchart(resec_size_tab.Outcome, log(resec_size_tab.(sprintf(comp))+1), 'MarkerStyle','none')
    hold on
    swarmchart(resec_size_tab.Outcome, log(resec_size_tab.(sprintf(comp))+1), 'filled')
    hold off
    title(sprintf("log(%s+1) across outcome groups", comp))
    subplot(3,2,5)
    histogram(log(resec_size_tab(resec_size_tab.Outcome == out_grps(1),:).(sprintf(comp))+1))
    title(sprintf("Log(resection size) (%s)",out_grps(1)))
    ylim([0 14])
    subplot(3,2,6)
    histogram(log(resec_size_tab(resec_size_tab.Outcome == out_grps(2),:).(sprintf(comp))+1))
    title(sprintf("Log(resection size) (%s)",out_grps(2)))
    ylim([0 14])
    [~,p, ~,st] = ttest2(log(resec_size_tab(resec_size_tab.Outcome == out_grps(1),:).(sprintf(comp))+1),...
        log(resec_size_tab(resec_size_tab.Outcome == out_grps(2),:).(sprintf(comp))+1));
    sgtitle(sprintf("t(%d)=%.3f, p=%.3f",st.df ,st.tstat,p))
    t_resec_comp(t_resec_comp.Comparison == comp,2:end) = table(st.tstat,st.df,p);
    saveas(fig, sprintf("%s_%s_t_tests.%s", save_loc, comp, file_type))
end

%% AUC for separating outcome groups based on resection size
auc_tabs.(sprintf(chan_or_roi)).resection_size  = outcome_auc(resec_size_tab,comp_vars,...
    resec_size_tab.Outcome == out_grps(2), "n_perm",n_perm,"plot_fig",1,"save_fig",save_fig,...
    "save_loc",save_loc, "file_type", "svg", "main_title", "AUCs comparing outcome groups based on resection size");

%% Logistic regression model
% OUtcome_group ~ proportion resected + size of resection


%% Group medians to report in paper
comps =  string(resec_size_tab.Properties.VariableNames);
med_tab = table(comps(3:end)', nan(length(comps)-2,1), nan(length(comps)-2,1),...
    nan(length(comps)-2,1), 'VariableNames', ["Comparison", "median" + ["", out_grps]]);
for comp = comps(3:end)
    med = nan(1,3);
    med(1) = median(resec_size_tab.(sprintf(comp)), 'omitnan');
    for grp = 1:2
        med(grp+1) = median(resec_size_tab(resec_size_tab.Outcome == out_grps(grp),:).(sprintf(comp)), 'omitnan');
    end
    med_tab(med_tab.Comparison == comp,2:end) = array2table(med);
end

median_tabs.(sprintf(chan_or_roi)).resection_size = med_tab;