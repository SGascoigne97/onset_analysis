%% Is resecting a larger proportion of the clinically labelled onset 
% associated with more favourable outcomes?
% tables were created in is_onset_resected.m (see comp_meth_with_resec
% struct)

save_loc = "../figures/subquestions/larger_prop_resec/";

atl_inds = [70, 120, 250, 500];
atl = atlas(atl_inds == str2double(extractAfter(chan_or_roi, "_")),:);
vols = atl.vol{:};

% % Create half-violin plots
% for parc = ["chan", "roi_120"]
%     if parc == "chan"
%         comp_measures = "Perc";
%     else 
%         comp_measures = ["Perc", "Perc_vol"];
%     end
%     half_violin(comp_meth_with_resec, det_meths, comp_measures, "parc", parc,...
%         "onset_across", onset_across, "vis_plot", "on",...
%         "save_fig", 1, "save_loc", save_loc, "file_type", "svg")
% end

det_meths = ["clo", "imprint", "EI"];
f = figure();
f.Position = [500,500,2400,450];
for det_ind = 1:3
    subplot(1,3,det_ind)
    half_violin_upd(comp_meth_with_resec.across_sz.chan.(sprintf(det_meths(det_ind))),...
        "Perc", comp_meth_with_resec.across_sz.chan.(sprintf(det_meths(det_ind))).outcome_cat,...
        "new_fig",0, "grp_names", ["ILAE 1-2", "ILAE 3+"], "y_lim",[0,1.2])
    xlim([-10,10])
    title(det_meths(det_ind))
end
saveas(f, sprintf("%sviolin_clo.svg", save_loc))
%%
comp_vars = ["Perc", "Perc_vol"];
f2 = figure();
f2.Position = [500,500,800,800];
for det_ind = 1:3
    for comp_ind = 1:2
        subplot(2,3,((comp_ind-1)*3)+det_ind)
        half_violin_upd(comp_meth_with_resec.across_sz.roi_120.(sprintf(det_meths(det_ind))),...
            comp_vars(comp_ind), comp_meth_with_resec.across_sz.roi_120.(sprintf(det_meths(det_ind))).outcome_cat,...
            "new_fig",0, "grp_names", ["ILAE 1-2", "ILAE 3+"], "y_lim",[0,1.2])
        title(strrep(sprintf("%s %s", det_meths(det_ind), comp_vars(comp_ind)), "_"," "))
        xlim([-15,15])
    end
end

saveas(f2, sprintf("%sviolin_roi_120.svg", save_loc))




%% Compute AUC and p-values comparing outcome groups
n_perm = 5000;
resec_thresh = 0.5;

for chan_or_roi = ["chan", "roi_120"]
    for det_meth = ["clo", "imprint"]
        if chan_or_roi == "chan"
            comp_measures = "Perc";
        else
            comp_measures = ["Perc", "Perc_vol"];
        end

        for comp_measure = comp_measures
            f = figure();
            f.Position = [10,10,1000,500];
            tiledlayout(1,3)
            data_tab = comp_meth_with_resec.across_sz.(sprintf(chan_or_roi)).(sprintf(det_meth));
            expl = data_tab.(sprintf(comp_measure));
            out = data_tab.Outcome>2;
            mod = fitglm(expl, out,'Distribution','binomial','Link','logit');

            probs = mod.Fitted.Probability;
            [X,Y,T,AUC] = perfcurve(out,probs,1);
        
            perm_auc = nan(n_perm,1);
            for perm = 1:n_perm
                rng(perm)
                expl_perm = expl(randperm(length(expl)));
                mod_perm = fitglm(expl_perm, out,...
                    'Distribution','binomial','Link','logit');
                probs_perm = mod_perm.Fitted.Probability;
                [~,~,~,AUC_perm] = perfcurve(out,probs_perm,1);
                perm_auc(perm) = AUC_perm;
            end
        
            nexttile
            plot(X,Y)
            hold on
            plot([0,1], [0,1])
            hold off
            xlabel('False positive rate') 
            ylabel('True positive rate')
            title(sprintf('AUC = %.3f', AUC))
        
            nexttile
            histogram(perm_auc, 'BinWidth', 0.025)
            hold on
            xline(AUC)
            hold off
            title(sprintf("p = %.3f", mean(AUC<perm_auc)))

            nexttile
            scatter(log(probs./(1-probs)), expl)
            lsline()
            title("Linearity assumption")

            sgtitle(strrep(sprintf("%s %s %s", chan_or_roi, det_meth,comp_measure), "_", " "))
            saveas(f, sprintf("%s%s_%s_%s_AUC.svg", save_loc,chan_or_roi, det_meth,comp_measure ))

            %% Perform chi-square test to determine if there is a signifcant
            % association between 'most' being resected and outcome
            f2 = figure();
            f2.Position = [10,10,600,500];
            out_str = string();
            out_str(out==1) = "Bad";
            out_str(out==0) = "Good";
            resec_str = string();
            resec_str(expl>=resec_thresh) = "Yes";
            resec_str(expl<resec_thresh) = "No";
            [contingency_tab, chi2, p, lab] = crosstab(resec_str, out_str);
            out_tab = table(resec_str', out_str');
            out_tab.Properties.VariableNames = ["resec_thresh", "Outcome"];
            % Create heatmap
            heatmap(out_tab, "resec_thresh", "Outcome")
            colorbar off

            % Perform test 
            if any(any(contingency_tab <5))
                [~, p, ~] = fishertest(contingency_tab);
                test_name = "Fishers exact";
                title(sprintf("%s (thresh = %d%%) \n %s (%s p=%.3f)", strrep(comp_measure, "_", " "),resec_thresh*100, det_meth, test_name, p))
          
            else 
                test_name = "Chi-square";
                title(sprintf("%s (thresh = %d%%) \n %s (%s X^2=%.2f, p=%.3f)", strrep(comp_measure, "_", " "),resec_thresh*100, det_meth, test_name, chi2, p))
          
            end
            saveas(f2, sprintf("%scontingency_%s_%s_%s_%d.%s", save_loc, chan_or_roi, det_meth,comp_measure,resec_thresh*100,file_type))
  
        end
    end
end

%% Look at association between proportion resected based on count of regions and volumes 
for det_meth = det_meths
    comp_meth_tab = comp_meth_with_resec.across_sz.roi_120.(sprintf(det_meth));

    f3 = figure("Position", [100,1000,800,800]);
    sgtitle(det_meth)
    subplot(3,3,[2,3,5,6])
    hold on
    scatter(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc,...
        comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc_vol,...
        'b',"filled", "XJitter","randn", "YJitter","randn",...
        "XJitterWidth",0.1,"YJitterWidth",0.1)
    scatter(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc,...
        comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc_vol, ...
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
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc_vol, "FaceColor", 'b', "BinWidth",0.05)
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc_vol, "FaceColor", 'r', "BinWidth",0.05)
    hold off
    xlim([0,1.1])
    set(gca, 'view',[90, -90],'YDir','reverse')
    axis off
  

    subplot(3,3,[8,9])
    hold on
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc, "FaceColor", 'b', "BinWidth",0.05)
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc, "FaceColor", 'r', "BinWidth",0.05)
    hold off
    xlim([0,1.1])
    set(gca, 'YDir','reverse')
    axis off

    saveas(f3, sprintf("%scomp_count_vol_%s.svg", save_loc, det_meth))

end
