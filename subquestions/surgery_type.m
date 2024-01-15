%% Chi square test to determine association between surgery type and outcome group
final_output.outcome_cat = categorical(final_output.outcome>2,[0,1], ["ILAE 1-2", "ILAE 3+"]);
ep_type_tab = final_output(string(final_output.("Op type")) == "F Lx" | string(final_output.("Op type")) == "T Lx",:);

[~,chi2,p,~] = crosstab(ep_type_tab.outcome_cat, ep_type_tab.("Op type"));

figure()
heatmap(ep_type_tab,'outcome_cat','Op type')
colorbar off
title(sprintf("chi^2(1,55) = %.3f, p = %.3f", chi2, p))

%% Repeating previous analyses including surgery type as a covariate
% Comparing proportion of onset resected
comp_meth_tab_op = innerjoin(comp_meth_tab, final_output(:,[1,27]), 'Keys','Patient_id');
comp_meth_tab_op = comp_meth_tab_op(string(comp_meth_tab_op.("Op type")) == "F Lx" |...
    string(comp_meth_tab_op.("Op type")) == "T Lx",:);

figure()
sgtitle("Comparing proportion resected across outcome groups and surgery types")
tiledlayout(2,2)
for op = ["T Lx", "F Lx"]
    for comp = ["Perc", "Perc_vol"]
        nexttile
        hold on
        boxchart(categorical(comp_meth_tab_op(comp_meth_tab_op.("Op type") == op,:).outcome_cat),...
            comp_meth_tab_op(comp_meth_tab_op.("Op type") == op,:).(sprintf(comp)), "MarkerStyle","none")
        swarmchart(categorical(comp_meth_tab_op(comp_meth_tab_op.("Op type") == op,:).outcome_cat),...
            comp_meth_tab_op(comp_meth_tab_op.("Op type") == op,:).(sprintf(comp)), "filled")
        hold off
        title(sprintf("%s for %s", comp, op))
    end
end


%%
figure()
sgtitle("Comparing proportion resected across outcome groups and surgery types")
tiledlayout(1,2)
for comp = ["Perc", "Perc_vol"]
    nexttile
    hold on
    boxchart(comp_meth_tab_op.outcome_cat, comp_meth_tab_op.(sprintf(comp)), "MarkerStyle","none")
    for op = ["T Lx", "F Lx"]
        swarmchart(categorical(comp_meth_tab_op(comp_meth_tab_op.("Op type") == op,:).outcome_cat),...
            comp_meth_tab_op(comp_meth_tab_op.("Op type") == op,:).(sprintf(comp)), "filled")
    end
        hold off
        title(sprintf("%s", comp))
end
legend([NaN, "TLE","FLE"], "Location","northeastoutside")
%%
comp_meth_tab_op.("Op_type") = categorical(comp_meth_tab_op.("Op type"));
comp_meth_tab_op.out_cat_bin = comp_meth_tab_op.Outcome>2;

mod_t = table(nan(2,1), nan(2,1), 'VariableNames',["Perc", "Perc_vol"], 'RowNames', ["variable", "op_type"]);
mod_p = mod_t;

n_perm = 2000;
for comp_measure = ["Perc", "Perc_vol"]
    f = figure();
    f.Position = [10,10,1000,500];
    tiledlayout(1,2)

    modelspec = sprintf('out_cat_bin ~ %s + Op_type',comp_measure);
    mod = fitglm(comp_meth_tab_op, modelspec,'Distribution','binomial','Link','logit', 'CategoricalVars', "Op_type");
    probs = mod.Fitted.Probability;
    [X,Y,T,AUC] = perfcurve(comp_meth_tab_op.out_cat_bin,probs,1);

    % Extract t-statistics and p-values to determine if there is evidence
    % of significant associations
    mod_t.(sprintf(comp_measure)) = mod.Coefficients.tStat(2:3);
    mod_p.(sprintf(comp_measure)) = mod.Coefficients.pValue(2:3);

    perm_auc = nan(n_perm,1);
    for perm = 1:n_perm
        comp_meth_tab_op_perm = comp_meth_tab_op;
        
        rng(perm)
        expl = comp_meth_tab_op_perm.(sprintf(comp_measure));
        comp_meth_tab_op_perm.(sprintf(comp_measure)) = expl(randperm(length(expl)));
        mod_perm = fitglm(comp_meth_tab_op_perm, modelspec,'Distribution','binomial','Link','logit', 'CategoricalVars', "Op_type");
        probs_perm = mod_perm.Fitted.Probability;
        [~,~,~,AUC_perm] = perfcurve(comp_meth_tab_op_perm.out_cat_bin,probs_perm,1);
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

    sgtitle(strrep(sprintf("%s",comp_measure), "_", " "))
    %saveas(f, sprintf("%s%s_AUC.svg", save_loc,comp_measure ))

end

%% Section 3.3: Comparing resection size
resec_size_tab = resec_tab;
resec_size_tab.Properties.VariableNames{1} = 'Patient_id';
resec_size_tab = innerjoin(resec_size_tab, final_output(:,[1,27,37]), 'Keys','Patient_id');
resec_size_tab = resec_size_tab(string(resec_size_tab.("Op type")) == "F Lx" |...
    string(resec_size_tab.("Op type")) == "T Lx",:);

comp_op_type_box(resec_size_tab, ["resec_count", "resec_vol"])
% Checking AUC using model that takes epilepsy type into account
% AUC for separating outcome groups based on onset size
resec_size_tab.out_cat_bin = resec_size_tab.outcome>2;
resec_size_tab = resec_size_tab(~isnan(resec_size_tab.resec_count),:);
mod_summaries.resec = comp_log_reg_auc(resec_size_tab,"out_cat_bin", ["resec_count", "resec_vol"], "n_perm",1000,"plot_fig",1);


%% Section 3.4: Comparing onset size (based on clo)
onset_size_tab = clo_tab;
onset_size_tab.Properties.VariableNames{1} = 'Patient_id';
onset_size_tab = innerjoin(onset_size_tab, final_output(:,[1,27,37]), 'Keys','Patient_id');
onset_size_tab = onset_size_tab(string(onset_size_tab.("Op type")) == "F Lx" |...
    string(onset_size_tab.("Op type")) == "T Lx",:);

comp_op_type_box(onset_size_tab, ["clo_count", "clo_vol"])

% Checking AUC using model that takes epilepsy type into account
% AUC for separating outcome groups based on onset size
onset_size_tab.out_cat_bin = onset_size_tab.outcome>2;
onset_size_tab = onset_size_tab(~isnan(onset_size_tab.clo_count),:);
mod_summaries.clo_onset_size = comp_log_reg_auc(onset_size_tab,"out_cat_bin", ["clo_count", "clo_vol"], "n_perm",1000,"plot_fig",1);

%% Comparing onset size (based on imprint)
onset_size_tab = impr_tab;
onset_size_tab.Properties.VariableNames{1} = 'Patient_id';
onset_size_tab = innerjoin(onset_size_tab, final_output(:,[1,27,37]), 'Keys','Patient_id');
onset_size_tab = onset_size_tab(string(onset_size_tab.("Op type")) == "F Lx" |...
    string(onset_size_tab.("Op type")) == "T Lx",:);

comp_op_type_box(onset_size_tab, ["impr_count", "impr_vol"])

% Checking AUC using model that takes epilepsy type into account
% AUC for separating outcome groups based on onset size
onset_size_tab.out_cat_bin = onset_size_tab.outcome>2;
onset_size_tab = onset_size_tab(~isnan(onset_size_tab.impr_count),:);
mod_summaries.imprint_onset_size = comp_log_reg_auc(onset_size_tab,"out_cat_bin", ["impr_count", "impr_vol"], "n_perm",1000,"plot_fig",1);


%% Section 3.5: Diffuse onsets
diffuse_tab = summary_tab;
diffuse_tab.Properties.VariableNames{1} = 'Patient_id';
diffuse_tab = innerjoin(diffuse_tab, final_output(:,[1,27,37]), 'Keys','Patient_id');
diffuse_tab = diffuse_tab(string(diffuse_tab.("Op type")) == "F Lx" |...
    string(diffuse_tab.("Op type")) == "T Lx",:);

% Median comparisons
med_comps = ["med_max_dist","med_onset_vol",...
    "subj_med_prop_dist","subj_med_prop_vol"];
comp_op_type_box(diffuse_tab, med_comps)
% Checking AUC using model that takes epilepsy type into account
% AUC for separating outcome groups based on onset size
diffuse_tab.out_cat_bin = diffuse_tab.outcome>2;
diffuse_tab = diffuse_tab(~isnan(diffuse_tab.med_max_dist),:);
mod_summaries.resec = comp_log_reg_auc(diffuse_tab,"out_cat_bin", med_comps, "n_perm",1000,"plot_fig",1);

% Maximum comparisons
max_comps = ["max_max_dist","max_onset_vol",...
    "subj_max_prop_dist","subj_max_prop_vol"];
comp_op_type_box(diffuse_tab, max_comps)
% Checking AUC using model that takes epilepsy type into account
% AUC for separating outcome groups based on onset size
diffuse_tab.out_cat_bin = diffuse_tab.outcome>2;
diffuse_tab = diffuse_tab(~isnan(diffuse_tab.max_max_dist),:);
mod_summaries.resec = comp_log_reg_auc(diffuse_tab,"out_cat_bin", max_comps, "n_perm",1000,"plot_fig",1);
%%

