function [auc_tab] = outcome_auc(data_tbl, comp_measures, outcome, opts)
    % Compute AUC and associated p-value for logistic regression models
    % predicting surgical outcome
    
    % Sarah Jane Gascoigne
    % 29/08/2023
    
    % input:
    %   - data_tbl: full data table
    
    % output
    %   - 
    
    arguments
        data_tbl table
        comp_measures (1,:) string 
        outcome (:,1) double % Binary separating subjects into outcome groups 
        opts.n_perm (1,1) double {mustBePositive} = 5000
        opts.plot_fig (1,1) double {mustBeMember(opts.plot_fig, [0,1])} = 1
        opts.save_fig (1,1) double {mustBeMember(opts.save_fig, [0,1])}  = 0
        opts.save_loc (1,1) string = "figures/"
        opts.file_type (1,1) string {mustBeMember(opts.file_type, ["png", "svg"])} = "svg"
        opts.main_title (1,1) string = "AUCs for comparing outcome groups"
    end

    if ~all(contains(comp_measures, data_tbl.Properties.VariableNames))
        fprintf("comp_measures must be variable names within data_tbl \n")
        auc_tab = NaN;
        return
    end
    
    %fill in optional arguments
    n_perm = opts.n_perm;
    plot_fig = opts.plot_fig;
    save_fig = opts.save_fig;
    save_loc = opts.save_loc;
    file_type = opts.file_type;
    main_title = opts.main_title;

    % Create a table to store outputs
    auc_tab = table(comp_measures', nan(length(comp_measures),1),...
        nan(length(comp_measures),1), 'VariableNames', ["Comparison", "AUC", "p"]);

    if plot_fig == 1
        f = figure();
        f.Position = [10,10,1000,1500];
        tiledlayout(length(comp_measures),2)
    end

    for comp = comp_measures
        expl = data_tbl.(sprintf(comp));
        outcome_clean = outcome(~isnan(expl));
        expl_clean = expl(~isnan(expl));
        mod = fitglm(expl_clean, outcome_clean,...
        'Distribution','binomial','Link','logit');
    
        probs = mod.Fitted.Probability;
        [X,Y,~,AUC] = perfcurve(outcome_clean,probs,1);

        if AUC <0.5
            AUC = 1-AUC;
        end
    
        perm_auc = nan(n_perm,1);
        for perm = 1:n_perm
            rng(perm) % set seed to ensure reproducible permutations and results
            expl_perm = expl_clean(randperm(length(expl_clean)));
            mod_perm = fitglm(expl_perm, outcome_clean,...
                'Distribution','binomial','Link','logit');
            probs_perm = mod_perm.Fitted.Probability;
            [~,~,~,AUC_perm] = perfcurve(outcome_clean,probs_perm,1);
            if AUC_perm <0.5
                AUC_perm = 1-AUC_perm;
            end
            perm_auc(perm) = AUC_perm;
        end

        % p-value is the proportion of permuted AUC values with higher
        % value than observed AUC

        if AUC>=0.5
            p_val = mean(AUC<perm_auc);
        else
            p_val = mean(AUC>perm_auc);
        end

        auc_tab(auc_tab.Comparison == comp,:).AUC = AUC;
        auc_tab(auc_tab.Comparison == comp,:).p = p_val;
    
        if plot_fig == 1
            nexttile
            plot(X,Y)
            hold on
            plot([0,1], [0,1])
            hold off
            xlabel('False positive rate') 
            ylabel('True positive rate')
            title(sprintf('%s \n AUC = %.3f', strrep(comp, "_", " "), AUC))
        
            nexttile
            histogram(perm_auc, 'BinWidth', 0.025)
            hold on
            xline(AUC)
            hold off
            title(sprintf("p = %.3f", p_val))

%             nexttile
%             P = mod.Fitted.Probability;
%             scatter(expl, log(P./(1-P)), "filled", "XJitter","rand", "YJitter", "rand", "XJitterWidth", range(expl)/25, "YJitterWidth", range(log(P./(1-P)))/25)
%             lsline()
%             xlabel(comp)
%             ylabel("Log(P/1-P)")
%             title("Linearity assumption")
        end
    end

    if plot_fig == 1
        sgtitle(main_title)
        if save_fig == 1
            save_loc = char(save_loc);
            saveas(f, sprintf("%sAUC.%s", save_loc, file_type))
        end
    end
end