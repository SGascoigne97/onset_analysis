%% 3.3 Larger onsets are not associated with surgical outcomes 

% Here we will compute the size of seizure onsets (both CLO and
% automatically captured) using both count of regions and volume (based on
% controls) and compare across surgical outcome groups

save_loc = [save_folder 'onset_size/'];
mkdir(save_loc)
save_loc = [save_loc  char(chan_or_roi) '_'];

tab_fill = nan(size(final_output,1),1);
% Create a table to store output (just counts for channels, counts and
% volumes for rois)
if chan_or_roi == "chan"
    onset_size_tab = table(final_output.Patient_id, final_output.outcome_cat,...
        tab_fill, tab_fill, tab_fill, 'VariableNames',...
        ["Patient_id", "Outcome", det_meths + "_count"]);
else
    onset_size_tab = table(final_output.Patient_id, final_output.outcome_cat,...
        tab_fill, tab_fill, tab_fill, tab_fill, tab_fill, tab_fill, ...
        'VariableNames', ["Patient_id", "Outcome", ...
        repelem(det_meths,2)+ repmat(["_count", "_vol"],1,3)]);
end

for subj = 1:size(final_output)
    subj_onset = final_output(subj,:);
   
    for meth = det_meths
        onset_mat = subj_onset.(sprintf("%s_%s", meth, chan_or_roi)){:};
%         onset_cons = mean(onset_mat,2) >= consensus_thresh;
        onset_cons = mean(onset_mat(:,sum(onset_mat,1)>0),2)>=consensus_thresh;
        if sum(onset_cons) >0
            onset_size_tab(subj,:).(sprintf("%s_count", meth)) = sum(onset_cons);
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

onset_size_tab.impl_size =...
     cellfun(@length, final_output.(sprintf("imprint_%s", parc)));


onset_parc_tab.(sprintf(chan_or_roi)) = onset_size_tab;

%% VISUALISATIONS
% Creating half-violin plots comparing onset sizes between outcome groups
comp_vars = string(onset_size_tab.Properties.VariableNames);
comp_vars = comp_vars(3:end);
if chan_or_roi == "chan"
    comp_vars = comp_vars(contains(comp_vars, "count"));
    half_violin(onset_size_tab, comp_vars,...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [1,1,1],"save_fig",save_fig,...
        "y_lim", repelem([NaN, NaN],3,1),...
        "file_type","svg", "save_loc", sprintf("%sviolins_count", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
else
    comp_vars = [comp_vars(contains(comp_vars, "count")), comp_vars(contains(comp_vars, "vol"))];
     half_violin(onset_size_tab, comp_vars(1:3),...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [1,1,1],"save_fig",save_fig,...
        "y_lim", repelem([0,20],3,1),...
        "file_type","svg", "save_loc", sprintf("%sviolins_count", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
    half_violin(onset_size_tab, comp_vars(4:6),...
        double(onset_size_tab.Outcome == out_grps(2)), ...
        "bin_wds", [NaN, NaN, NaN],"save_fig",save_fig,...
        "y_lim", repelem([0,4*10^5],3,1),...
        "file_type","svg", "save_loc", sprintf("%sviolins_vol", save_loc), ...
        "grp_names", out_grps, "legend_pos", "northeast")
end

%% T-tests to compare mean onset size across outcome groups
t_onset_comp = table(comp_vars', nan(length(comp_vars),1), nan(length(comp_vars),1),...
    nan(length(comp_vars),1),'VariableNames',["Comparison", "t", "df", "p"]);
for comp = comp_vars
    data_tab = table(onset_size_tab.Patient_id, onset_size_tab.Outcome,...
        onset_size_tab.(sprintf(comp)), 'VariableNames',["Patient_id", ...
        "Outcome", "Comp"]);
    data_tab = data_tab(~isnan(data_tab.Comp),:);

    if comp == "imprint_count"
        fprintf("T-test for imprint count is not possible due to non-normality \n")
        continue
    end

    fig = figure("Position",[10,10,500,900]);
    subplot(3,2,1:2)
    boxchart(data_tab.Outcome, data_tab.Comp, 'MarkerStyle','none')
    hold on
    swarmchart(data_tab.Outcome, data_tab.Comp, 'filled')
    hold off
    title(strrep(sprintf("%s", comp), "_", " "))
    subplot(3,2,3:4)
    boxchart(data_tab.Outcome, log(data_tab.Comp+1), 'MarkerStyle','none')
    hold on
    swarmchart(data_tab.Outcome, log(data_tab.Comp+1), 'filled')
    hold off
    title(strrep(sprintf("log(%s+1)", comp), "_", " "))

    subplot(3,2,5)
    histogram(log((data_tab(data_tab.Outcome == out_grps(1),:).Comp)+1))%, 'BinWidth',bin_wd)
    title(strrep(sprintf("log(%s+1) (%s)", comp, out_grps(1)), "_", " "))
    ylim([0 15])
    
    subplot(3,2,6)
    histogram(log((data_tab(data_tab.Outcome == out_grps(2),:).Comp)+1))%, 'BinWidth',bin_wd)
    title(strrep(sprintf("log(%s+1) (%s)", comp, out_grps(2)), "_", " "))
    ylim([0 15])
    [~,p, ~,st] = ttest2(log(data_tab(data_tab.Outcome == out_grps(1),:).Comp+1),...
        log(data_tab(data_tab.Outcome == out_grps(2),:).Comp+1));
    sgtitle(sprintf("t(%d)=%.3f, p=%.3f",st.df ,st.tstat,p))
    t_onset_comp(t_onset_comp.Comparison == comp,2:end) = table(st.tstat,st.df,p);
    saveas(fig, sprintf("%s_%s_t_tests.%s", save_loc, comp, file_type))
end

t_tabs.(sprintf(chan_or_roi)).onset_size = t_onset_comp;

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

median_tabs.(sprintf(chan_or_roi)).resection_size = med_tab;

%% SUPPLEMENTARY
% Compare volumes between CLO and automatic consensus onsets (paired t)
if chan_or_roi == "chan"
    comps = "count";
else
    comps = ["count", "vol"];
end
for comp = comps
    clo_size = log(onset_size_tab.(sprintf("clo_%s", comp))+1);
    for det = det_meths(2:end)
        auto_size = log(onset_size_tab.(sprintf("%s_%s", det, comp))+1);

        fig = figure('Position',[10,10,1000,500]);
        sgtitle(sprintf("Comparing CLO and %s consensus onset volumes", det))
        subplot(121)
        diffs = auto_size-clo_size;
        diffs_clean = diffs(~isnan(diffs));
        histogram(diffs_clean)
        title("Histogram of differences")
        subplot(122)
        boxplot([clo_size, auto_size])
        hold on
        for i = 1:length(auto_size)
            plot([1,2],[clo_size(~isnan(diffs)),auto_size(~isnan(diffs))],'-x')
        end
        hold off
        set(gca, "XTickLabel", ["CLO", "Consensus"])
        [~,p, ~,st] = ttest(diffs_clean);
        title(sprintf("%s t(%d) = %.2f, p = %.2f", comp, length(diffs_clean)-1, st.tstat, p ))
        saveas(fig, sprintf("%sclo_vs_%s_%s.%s", save_loc, det, comp, file_type))
    end
end
