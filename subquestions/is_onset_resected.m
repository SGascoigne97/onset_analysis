%% 3.1 Seizure onset tends to be resected, however resecting a larger 
% proportion of the SOZ is not associated with more favourable surgical 
% outcomes
save_loc = [save_folder 'onset_resected/'];
mkdir(save_loc)
save_loc = [save_loc  char(chan_or_roi) '_'];

% Create structure with all comparisons between onset (CLO and automatic)
% and resection

if contains(chan_or_roi, "roi")
    atl = atlas(atl_inds == str2double(extractAfter(chan_or_roi, "_")),:);
    vols = atl.vol{:};
    comp_measures = ["Perc", "Perc_vol"];
else % If analysis is channel-wise, we cannot compute volume
    atl = NaN;
    vols = NaN;
    comp_measures = "Perc";
end

for det = det_meths
    fprintf("%s \n", det)
    comp_meth_tab = compute_resec_comp(final_output,atl,vols,...
        "chan_or_roi", chan_or_roi, "det_meth", det,...
        "consensus_thresh",consensus_thresh,"onset_across",1);
    comp_with_resec.(sprintf("%s", chan_or_roi)).(sprintf("%s", det)) = comp_meth_tab;
    clear comp_meth_tab
end

%% Chi-square test between outcome groups and whether subject had 'most' or all
% of onset resected
tab_fill = nan(length(comp_measures)*length(det_meths),1);
onset_resec_chi2 = table(repelem(det_meths,length(comp_measures))',...
    repmat(comp_measures,1,length(det_meths))',tab_fill,tab_fill,...
    tab_fill,tab_fill, 'VariableNames', ["det_meth", "comp", "chi2", "df", "n", "p"]);

for amount_resec = ["most", "all"]
    for det = det_meths    
        comp_tab = comp_with_resec.(sprintf("%s", chan_or_roi)).(sprintf("%s", det));
        out = comp_tab.Outcome;
        for comp_meas = comp_measures
            data = comp_tab.(sprintf(comp_meas));
            if amount_resec == "most"
                comp_tab.most_resec = categorical(data>=most_resec_thresh, [0,1], ["No", "Yes"]);
            else
                comp_tab.most_resec = categorical(data==1, [0,1], ["No", "Yes"]);
            end
    
            [cont_tab,chi2,p,cats] = crosstab(comp_tab.most_resec, comp_tab.Outcome);
            df = (size(cats,1)-1)*(size(cats,2)-1);
            fig = figure();
            heatmap(comp_tab, "most_resec", "Outcome")
            colorbar off
            
            % We will use a chi^2 test if all values are greater than 5,
            % otherwise we will use Fisher's exact test
            if any(any(cont_tab <5))
                [~, p, ~] = fishertest(cont_tab);
                title(strrep(sprintf("%s %s (%s) \n Fisher's exact test p = %.2f",...
                    det, comp_meas,amount_resec,p), "_", " "))
                % Save in outcome table
                onset_resec_chi2(onset_resec_chi2.comp == comp_meas &...
                            onset_resec_chi2.det_meth == det,3:end) = ...
                            table(NaN, NaN, size(comp_tab,1),p);
    
            else 
                title(strrep(sprintf("%s %s (%s)\nchi^2(%d,%d) = %.2f, p = %.2f",...
                    det, comp_meas, amount_resec, df, size(comp_tab,1), chi2,p), "_", " "))
                % Save in outcome table
                onset_resec_chi2(onset_resec_chi2.comp == comp_meas &...
                            onset_resec_chi2.det_meth == det,3:end) = ...
                            table(chi2, df, size(comp_tab,1),p);
            end
           
            if save_fig == 1
                saveas(fig,sprintf("%s%s_%s_chi2_%s.%s", save_loc, det, comp_meas, amount_resec, file_type))
            end
            clear fig
        end
    end
    
    chi2_tab.(sprintf(chan_or_roi)).(sprintf("%s_resec", amount_resec)) = onset_resec_chi2;
end
%% VISUALISATIONS
% Violin plots comparing outcome groups
for comp_ind = 1:length(comp_measures)
    fig = figure();
    fig.Position = [500,500,400*length(det_meths),450];
    for det_ind = 1:length(det_meths)
        comp_tab = comp_with_resec.(sprintf(chan_or_roi)).(sprintf(det_meths(det_ind)));
        subplot(1,3,det_ind)
        half_violin(comp_tab, comp_measures(comp_ind), comp_tab.Outcome,...
            "new_fig",0, "grp_names", out_grps, "y_lim", [0,1.2], "save_fig", 0)
        title(strrep(sprintf("%s %s (n=%d)", det_meths(det_ind), ...
            comp_measures(comp_ind), size(comp_tab,1)), "_"," "))
        clear comp_tab
        if save_fig == 1
            saveas(fig, sprintf("%s%s_violin.svg", save_loc, comp_measures(comp_ind)))
        end
    end
end

%% Compute AUC and p-values comparing outcome groups
for det = ["clo", "imprint", "EI"]
    if chan_or_roi == "chan"
        comp_measures = "Perc";
    else
        comp_measures = ["Perc", "Perc_vol"];
    end

    data_tab = comp_with_resec.(sprintf(chan_or_roi)).(sprintf(det));
    auc_tabs.(sprintf(chan_or_roi)).onset_resec.(sprintf(det)) = ...
        outcome_auc(data_tab, comp_measures, data_tab.Outcome == out_grps(2),...
            "n_perm", n_perm, "plot_fig",1, "save_fig",save_fig,...
            "file_type", file_type, "save_loc",save_loc+det+"_", "main_title", ...
            sprintf("AUCs for comparing proportion of onset (%s) resected", det));
end

%% SUPPLEMENTARY: Look at association between proportion resected based on count of regions and volumes 
if contains(chan_or_roi, "roi")
    for det = det_meths
        comp_meth_tab = comp_with_resec.roi_120.(sprintf(det));
    
        fig = figure("Position", [100,1000,800,800]);
        sgtitle(det)
        subplot(3,3,[2,3,5,6])
        hold on
        scatter(comp_meth_tab(comp_meth_tab.Outcome == out_grps(1),:).Perc,...
            comp_meth_tab(comp_meth_tab.Outcome == out_grps(1),:).Perc_vol,...
            'b',"filled", "XJitter","randn", "YJitter","randn",...
            "XJitterWidth",0.1,"YJitterWidth",0.1)
        scatter(comp_meth_tab(comp_meth_tab.Outcome == out_grps(2),:).Perc,...
            comp_meth_tab(comp_meth_tab.Outcome == out_grps(2),:).Perc_vol, ...
            'r',"filled", "XJitter","randn", "YJitter","randn",...
            "XJitterWidth",0.1,"YJitterWidth",0.1)
        plot([0,1], [0,1], 'k--')   
        hold off
        xlabel("Perc (binary)")
        ylabel("Perc (vol)")
        xlim([0,1.1])
        ylim([0,1.1])
    
        subplot(3,3,[1,4])
        hold on
        histogram(comp_meth_tab(comp_meth_tab.Outcome == out_grps(1),:).Perc_vol, "FaceColor", 'b', "BinWidth",0.05)
        histogram(comp_meth_tab(comp_meth_tab.Outcome == out_grps(2),:).Perc_vol, "FaceColor", 'r', "BinWidth",0.05)
        hold off
        xlim([0,1.1])
        set(gca, 'view',[90, -90],'YDir','reverse')
        axis off
      
        subplot(3,3,[8,9])
        hold on
        histogram(comp_meth_tab(comp_meth_tab.Outcome == out_grps(1),:).Perc, "FaceColor", 'b', "BinWidth",0.05)
        histogram(comp_meth_tab(comp_meth_tab.Outcome == out_grps(2),:).Perc, "FaceColor", 'r', "BinWidth",0.05)
        hold off
        xlim([0,1.1])
        set(gca, 'YDir','reverse')
        axis off
    
        if save_fig == 1
            saveas(fig, sprintf("%sscatter_count_vol_%s.%s", save_loc, det, file_type))
        end
        clear fig
    end
end
