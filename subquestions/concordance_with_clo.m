%% 3.1 IOLA onsets have moderate agreement with CLO in most subjects
% We will consider concordance between the consensus onset and the maximum
% consensus between the IOLA onset for any one seizure and the clinically
% labelled onset (CLO)

save_loc = [save_folder 'clo_concord/'];
mkdir(save_loc)
save_loc = [save_loc  char(chan_or_roi) '_'];

for det_meth = det_meths
    if det_meth == "clo"
        continue
    else
        % Create a table to store concordance values
        concord_tab.(sprintf(det_meth)) = table(final_output.Patient_id, nan(size(final_output,1), 1),...
            nan(size(final_output,1), 1), final_output.outcome_cat,...
            'VariableNames', ["Patient_id", "consensus", "max_clo_concord", "outcome"]);
        
        % Iterate through subjects
        for subj = 1:size(final_output)
            % Extract the CLO, and automatically detected onset
            clo = final_output(subj,:).(sprintf("clo_%s", chan_or_roi)){:};
            onset_mat = final_output(subj,:).(sprintf("%s_%s", det_meth, chan_or_roi)){:};
    
            % Compute consensus onset
            onset_cons = mean(onset_mat(:,sum(onset_mat,1)>0),2)>=consensus_thresh;
        
            % Skip subjects with no consensus onset or CLO
             if sum(clo) == 0
                 continue
             end
        
             % First we will compute Cohen's Kappa as measure of concordance
             % between the CLO and each seizure from which we will compute the
             % maximum
             per_sz_concord = nan(size(onset_mat,1),1);
             for sz = 1:size(onset_mat,2)
                per_sz_concord(sz) = cohensKappa(logical(onset_mat(:,sz)), logical(clo));
             end
             concord_tab.(sprintf(det_meth)).max_clo_concord(subj) = max(per_sz_concord);
        
             % Next, we will compute Cohen's Kappa between CLO and the consensus
             % onset
             % NOTE: we must exclude subjects with no consensus onset 
              if sum(onset_cons) == 0
                continue
             end
             concord_tab.(sprintf(det_meth)).consensus(subj) = cohensKappa(logical(clo), onset_cons);
        end
        
        % Compute AUC and p-values comparing outcome groups
        auc_tabs.(sprintf(chan_or_roi)).clo_concord.(sprintf(det_meth)) = outcome_auc(concord_tab.(sprintf(det_meth)), ["consensus", "max_clo_concord"],...
            concord_tab.(sprintf(det_meth)).outcome == out_grps(2), "n_perm",n_perm,"plot_fig",1,"save_fig",save_fig,...
            "save_loc",save_loc, "file_type", "svg", "main_title", sprintf("AUCs comparing outcome groups based on %s concordance with CLO %s", det_meth, chan_or_roi));
%   
        %% VISUALISATIONS 
        % Half violins plot distributions of markers for each outcome group 
        % (ILAE 1-2 compared against ILAE 3+)
        half_violin(concord_tab.(sprintf(det_meth)),  ["consensus", "max_clo_concord"],...
            concord_tab.(sprintf(det_meth)).outcome==out_grps(2), "grp_names", out_grps,...
            "y_lim", [-0.55, 1.05;-0.55, 1.05],"save_fig", save_fig, "file_type", "svg",...
            "save_loc", sprintf("%shalf_violins_%s",save_loc, det_meth),...
            "legend_pos", "southeast")
        
        % Chi-sq test comparing counts of subjects with moderate concordance
        % across outcome groups (not included in paper)
        f = figure("Position",[10,10,1000,500]);
        tiledlayout(1,2)
        for concord_summ = ["consensus", "max_clo_concord"]
            nexttile
            non_nan = ~isnan(concord_tab.(sprintf(det_meth)).(sprintf(concord_summ)));
            concord_var = categorical(concord_tab.(sprintf(det_meth)).(sprintf(concord_summ))>0.4,[0,1],["Below 0.4", "Above 0.4"]);
            outcome_var = categorical(concord_tab.(sprintf(det_meth)).outcome == out_grps(2),[0,1],out_grps);
            var_tab = table(concord_var(non_nan), outcome_var(non_nan), 'VariableNames', ["Concordance", "Outcome"]);
            [~,chi2,p,~] = crosstab(concord_var(non_nan), outcome_var(non_nan));
            heatmap(var_tab,"Concordance", "Outcome")
            colorbar off
            title(sprintf("%s %s (X^2 =%.2f , p = %.3f)", det_meth, concord_summ, chi2,p))
        end
        saveas(f, sprintf("%s%s_contingency_moderate_concord.%s", save_loc, det_meth, file_type))
    end
end

