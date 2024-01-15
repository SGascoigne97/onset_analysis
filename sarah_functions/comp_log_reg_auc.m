function [mod_details] = comp_log_reg_auc(data_tbl, response_var, comp_measures, opts)
% Compute AUC and associated p-value for logistic regression models
% predicting surgical outcome

% Sarah Jane Gascoigne
% 24/08/2023

% input:
%   - data_tbl: full data table

% output
%   - mod_details: struct with relevant information for logistic regression
%   models

    arguments
        data_tbl table
        response_var (1,1) string
        comp_measures (1,:) string
        opts.plot_fig {mustBeMember(opts.plot_fig, [0,1])} = 0
        opts.save_fig (1,1) double {mustBeNumeric} = 0
        opts.save_loc (1,1) string = "figures/"
        opts.n_perm (1,1) double = 1000
    end
    
    %fill in optional arguments
    plot_fig = opts.plot_fig;
    save_fig = opts.save_fig;
    save_loc = opts.save_loc;
    n_perm = opts.n_perm;

    data_tbl.("Op_type") = categorical(data_tbl.("Op type"));

for comp_measure = comp_measures
    
    modelspec = sprintf('%s ~ %s + Op_type', response_var, comp_measure);
    mod = fitglm(data_tbl, modelspec,'Distribution','binomial','Link','logit', 'CategoricalVars', "Op_type");
    probs = mod.Fitted.Probability;
    [X,Y,~,AUC] = perfcurve(data_tbl.(sprintf(response_var)),probs,1);

    % Extract model details to determine if there is evidence
    % of significant associations with surgical outcomes
    mod_details.(sprintf(comp_measure)).model = mod;
    mod_chisq = devianceTest(mod);
    mod_details.(sprintf(comp_measure)).chi2 = mod_chisq.chi2Stat(2);
    mod_details.(sprintf(comp_measure)).chi2_p = mod_chisq.pValue(2);

    perm_auc = nan(n_perm,1);
    for perm = 1:n_perm
        data_tbl_perm = data_tbl;
        rng(perm)
        expl = data_tbl_perm.(sprintf(comp_measure));
        data_tbl_perm.(sprintf(comp_measure)) = expl(randperm(length(expl)));
        mod_perm = fitglm(data_tbl_perm, modelspec,'Distribution','binomial','Link','logit', 'CategoricalVars', "Op_type");
        probs_perm = mod_perm.Fitted.Probability;
        [~,~,~,AUC_perm] = perfcurve(data_tbl_perm.(sprintf(response_var)),probs_perm,1);
        perm_auc(perm) = AUC_perm;
    end

      mod_details.(sprintf(comp_measure)).AUC = AUC;
      mod_details.(sprintf(comp_measure)).AUC_p = mean(perm_auc>AUC);

    if plot_fig == 1
        f = figure();
        f.Position = [10,10,1000,500];
        tiledlayout(1,2)
    
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
        title(sprintf("p = %.3f", mean(perm_auc>AUC)))
    
        sgtitle(strrep(sprintf("%s",comp_measure), "_", " "))
        if save_fig == 1
            saveas(f, sprintf("%s%s_AUC.svg", save_loc,comp_measure ))
        end
    end

end