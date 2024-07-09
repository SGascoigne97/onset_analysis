clear combined_tab
for chan_or_roi = ["chan", "roi_120", "roi_250"]
    for comp = ["clo_concord", "onset_resec"]
        fields = string(fieldnames(auc_tabs.(sprintf(chan_or_roi)).(sprintf(comp))));
        for field_ind = 1:length(fields)
            field = fields(field_ind);
            value_tab = auc_tabs.(sprintf(chan_or_roi)).(sprintf(comp)).(sprintf(field));
            field_tab = table(chan_or_roi, comp, field);
            if exist("combined_tab","var")
                combined_tab = [combined_tab; [repmat(field_tab,size(value_tab,1),1), value_tab]];
            else
                combined_tab = [repmat(field_tab,size(value_tab,1),1), value_tab];
            end
        end
    end

    for comp = ["onset_size", "resection_size"]
        comp_tab = auc_tabs.(sprintf(chan_or_roi)).(sprintf(comp));
        field_tab = table(extractBefore(comp_tab.Comparison, "_"), extractAfter(comp_tab.Comparison, "_"), 'VariableNames', ["field", "Comparison"]);
        field_tab = [repmat(table(chan_or_roi, comp),size(field_tab,1),1), field_tab];
        % Add to combined table
        combined_tab = [combined_tab; field_tab, comp_tab(:,2:3)];
        
    end
end

save_loc = "figures/paper_figures/supplementary";

%% Compare results across parcellation schemes

for comp = ["clo_concord", "onset_resec", "onset_size", "resection_size"]
    comp_tab = combined_tab(combined_tab.comp == comp,:);
end

%%
figure()
tiledlayout(3,3)
col_names = string(onset_parc_tab.roi_120.Properties.VariableNames(3:end));
col_names = col_names(~contains(col_names, "vol")); % We have excluded volume from our analyses

parcs = ["chan", "roi_120", "roi_250"];

for comp = col_names
    for parc_ind_1 = 1:length(parcs)
        for parc_ind_2 = 1:length(parcs)
            if parc_ind_1 >= parc_ind_2
                continue
            end
            parc1 = parcs(parc_ind_1);
            parc2 = parcs(parc_ind_2);
            nexttile
            gscatter(onset_parc_tab.(sprintf(parc1)).(sprintf(comp)),...
                onset_parc_tab.(sprintf(parc2)).(sprintf(comp)), onset_parc_tab.(sprintf(parc1)).Outcome ,...
                'filled')
            lsline()
            title(sprintf("Comparing count of %s regions", extractBefore(comp,"_")))
            xlabel(strrep(sprintf(parc1), "_", " "))
            ylabel(strrep(sprintf(parc2), "_", " "))
        end
    end
end

%%

% Paired sample t-test for proportion of onset resected between parcellation schemes 
% As distributions are guaranteerd to have different means across
% parcellations schemes, we will use a correlation analysis to determine if
% there are notable deviations from expected relationships between
% parcellation schemes
parcs = ["chan", "roi_120", "roi_250"];

%% Proportion of onset resected
comp_parc_tab.onset_resec = array2table(nan(0,7), 'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "t", "DoF", "p"]);

% Is onset resected (Section 3.2)
for det = det_meths
    resec_tab= table(repelem(parcs(1:length(parcs)-1), 1, length(parcs)*2)', ...
    repelem(repmat(parcs, 1, (length(parcs)-1)),1,2)', repmat(["Perc", "Perc_vol"], 1, (length(parcs)-1)*length(parcs))',...
    repmat(det, 2*(length(parcs)-1)*length(parcs), 1 ),...
    nan(2*(length(parcs)-1)*length(parcs),1), nan(2*(length(parcs)-1)*length(parcs),1), ...
    nan(2*(length(parcs)-1)*length(parcs),1),...
    'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "t", "DoF", "p"]);
    for parc1 = 1:length(parcs)
        for parc2 = 1:length(parcs)
            if parc1 == parc2 | parc1 > parc2
                continue
            else
                tab1 = comp_with_resec.(sprintf(parcs(parc1))).(sprintf(det));
                tab2 = comp_with_resec.(sprintf(parcs(parc2))).(sprintf(det));

                % Join tables together to match patients 
                across_parc_tab = innerjoin(tab1, tab2, 'Keys', ...
                    ["Patient_id", "Outcome", "Sz_count"]);

                if contains(parcs(parc1), "roi") & contains(parcs(parc1), "roi")
                    comps = ["Perc", "Perc_vol"];
                else
                    comps = "Perc";
                end

                for comp = comps
                    diffs = across_parc_tab.(sprintf("%s_tab1", comp)) - across_parc_tab.(sprintf("%s_tab2", comp));
                    figure()
                    sgtitle(strrep(sprintf("Comparing %s %s across %s and %s", det, comp, parcs(parc1), parcs(parc2) ), "_", " "))
                    subplot(121)
                    scatter(across_parc_tab.(sprintf("%s_tab1", comp)),...
                        across_parc_tab.(sprintf("%s_tab2", comp)), 'filled',...
                        'XJitter', 'randn', 'XJitterWidth', 0.05,...
                        'YJitter', 'randn', 'YJitterWidth', 0.05)
                    lsline()
                    xlabel(parcs(parc1))
                    ylabel(parcs(parc2))
                    subplot(122)
                    [~,p,~, stats] = ttest(diffs);
                    histogram(diffs)
                    title(sprintf("Histogram of differences \n t(%d) = %.3f, p = %.3f", length(diffs)-1, stats.tstat, p))

                    % Store test statistic, degrees of freedom, and
                    % p-values in table
                    resec_tab(resec_tab.Scheme_1 == parcs(parc1) &...
                        resec_tab.Scheme_2 == parcs(parc2) & ...
                        resec_tab.comp == comp & resec_tab.det == det, 5:7) = ...
                        table(stats.tstat, length(diffs)-1, p);
                end
            end
        end
    end

    resec_tab = resec_tab(~isnan(resec_tab.t),:);
    comp_parc_tab.onset_resec = [comp_parc_tab.onset_resec; resec_tab];
end


%% Onset size
comp_parc_tab.onset_size = array2table(nan(0,7), 'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "rho", "DoF", "p"]);

% Is onset resected (Section 3.2)
for det = det_meths
    resec_tab = table(repelem(parcs(1:length(parcs)-1), 1, length(parcs)*2)', ...
    repelem(repmat(parcs, 1, (length(parcs)-1)),1,2)', repmat(["count", "vol"], 1, (length(parcs)-1)*length(parcs))',...
    repmat(det, 2*(length(parcs)-1)*length(parcs), 1 ),...
    nan(2*(length(parcs)-1)*length(parcs),1), nan(2*(length(parcs)-1)*length(parcs),1), ...
    nan(2*(length(parcs)-1)*length(parcs),1),...
    'VariableNames', ["Scheme_1", "Scheme_2", "comp", "det", "rho", "DoF", "p"]);
    for parc1 = 1:length(parcs)
        for parc2 = 1:length(parcs)
            if parc1 == parc2 | parc1 > parc2
                continue
            else
                if contains(parcs(parc1), "roi") & contains(parcs(parc1), "roi")
                    comps = ["count", "vol"];
                else
                    comps = "count";
                end

                tab1 = onset_parc_tab.(sprintf(parcs(parc1)));
                tab2 = onset_parc_tab.(sprintf(parcs(parc2)));

                % Join tables together to match patients 
                across_parc_tab = innerjoin(tab1, tab2, 'Keys', ...
                    ["Patient_id", "Outcome"]);

                for comp = comps
                    grp_1_vals = across_parc_tab.(sprintf("%s_%s_tab1", det, comp));
                    grp_2_vals = across_parc_tab.(sprintf("%s_%s_tab2", det, comp));

                    rm_obs = isnan(grp_1_vals)| isnan(grp_2_vals);
                    grp_1_vals = grp_1_vals(~rm_obs);
                    grp_2_vals = grp_2_vals(~rm_obs);

                    f = figure();
                    f.Position = [10,10,700,900];
                    sgtitle(strrep(sprintf("Comparing %s %s across %s and %s", det, comp, parcs(parc1), parcs(parc2) ), "_", " "))
%                     subplot(121)
                    scatter(grp_1_vals, grp_2_vals,'filled',...
                        'XJitter', 'randn', 'XJitterWidth', range(grp_1_vals)/20,...
                        'YJitter', 'randn', 'YJitterWidth', range(grp_2_vals)/20)
                    lsline()    
                    axis square
                    xlabel(parcs(parc1))
                    ylabel(parcs(parc2))

                    [r, p] = corr(grp_1_vals, grp_2_vals, "type","Spearman");

                    title(sprintf("Spearman's rank \n r(%d) = %.3f, p = %.3f", length(grp_1_vals)-2, r, p))
                    saveas(f, sprintf('%s/onset_size_%s_%s(%s_v_%s).svg', save_loc, det, comp, parcs(parc1), parcs(parc2)))

                    % Store test statistic, degrees of freedom, and
                    % p-values in table
                    resec_tab(resec_tab.Scheme_1 == parcs(parc1) &...
                        resec_tab.Scheme_2 == parcs(parc2) & ...
                        resec_tab.comp == comp & resec_tab.det == det, 5:7) = ...
                        table(r, length(grp_1_vals)-2, p);
                end
            end
        end
    end

    resec_tab = resec_tab(~isnan(resec_tab.rho),:);
    comp_parc_tab.onset_size = [comp_parc_tab.onset_size; resec_tab];
end

%%
%% Resection size
comp_parc_tab.resection_size = table(repelem(parcs(1:length(parcs)-1), 1, length(parcs)*2)', ...
repelem(repmat(parcs, 1, (length(parcs)-1)),1,2)', repmat(["count", "vol"], 1, (length(parcs)-1)*length(parcs))',...
nan(2*(length(parcs)-1)*length(parcs),1), nan(2*(length(parcs)-1)*length(parcs),1), ...
nan(2*(length(parcs)-1)*length(parcs),1),...
'VariableNames', ["Scheme_1", "Scheme_2", "comp", "rho", "DoF", "p"]);
for parc1 = 1:length(parcs)
    for parc2 = 1:length(parcs)
        if parc1 == parc2 | parc1 > parc2
            continue
        else
            if contains(parcs(parc1), "roi") & contains(parcs(parc1), "roi")
                comps = ["count", "vol"];
            else
                comps = "count";
            end

            tab1 = onset_parc_tab.(sprintf(parcs(parc1)));
            tab2 = onset_parc_tab.(sprintf(parcs(parc2)));

            % Join tables together to match patients 
            across_parc_tab = innerjoin(tab1, tab2, 'Keys', ...
                ["Patient_id", "Outcome"]);

            for comp = comps
                grp_1_vals = across_parc_tab.(sprintf("%s_%s_tab1", det, comp));
                grp_2_vals = across_parc_tab.(sprintf("%s_%s_tab2", det, comp));

                rm_obs = isnan(grp_1_vals)| isnan(grp_2_vals);
                grp_1_vals = grp_1_vals(~rm_obs);
                grp_2_vals = grp_2_vals(~rm_obs);

                f = figure();
                f.Position = [10,10,700,900];
                sgtitle(strrep(sprintf("Comparing resection %s across %s and %s", comp, parcs(parc1), parcs(parc2) ), "_", " "))
%                     subplot(121)
                scatter(grp_1_vals, grp_2_vals, 'filled')
                axis square
                lsline()          
                xlabel(parcs(parc1))
                ylabel(parcs(parc2))

                [r, p] = corr(grp_1_vals, grp_2_vals, "type","Spearman");

                title(sprintf("Spearman's rank \n r(%d) = %.3f, p = %.3f", length(grp_1_vals)-2, r, p))
                saveas(f, sprintf('%s/resection_%s(%s_v_%s).svg', save_loc, comp, parcs(parc1), parcs(parc2)))

                % Store test statistic, degrees of freedom, and
                % p-values in table
                comp_parc_tab.resection_size(comp_parc_tab.resection_size.Scheme_1 == parcs(parc1) &...
                    comp_parc_tab.resection_size.Scheme_2 == parcs(parc2) & ...
                    comp_parc_tab.resection_size.comp == comp, 4:6) = ...
                    table(r, length(grp_1_vals)-2, p);
            end
        end
    end
end

comp_parc_tab.resection_size = comp_parc_tab.resection_size(~isnan(comp_parc_tab.resection_size.rho),:);

%% Compare AUCs across parcellation schemes
parcs = ["chan", "roi_120", "roi_250"];
det_meths = ["clo", "imprint", "EI"]; 
comps = ["onset_resec", "onset_size", "resection_size", "clo_concord"]; % We will look at "resection_size" separately for now 

for comp = comps
    auc_comp_mat = nan(length(det_meths), length(parcs));
    p_comp_mat = auc_comp_mat;

    for det_meth = det_meths
        for parc = parcs
            row = det_meths == det_meth;
            col = parcs == parc;
            
            if comp == "onset_resec" 
                auc_and_p = auc_tabs.(sprintf(parc)).(sprintf(comp)).(sprintf(det_meth));
                auc_and_p = auc_and_p(~contains(auc_and_p.Comparison, "_vol"),:);
                auc_comp_mat(row,col) = auc_and_p.AUC;
                p_comp_mat(row,col) = auc_and_p.p;

            elseif comp == "clo_concord"
                if det_meth == "clo"
                    continue
                end
                auc_and_p = auc_tabs.(sprintf(parc)).(sprintf(comp)).(sprintf(det_meth));
                auc_and_p_consensus = auc_and_p(contains(auc_and_p.Comparison, "consensus"),:);
                auc_comp_mat(row,col) = auc_and_p_consensus.AUC;
                p_comp_mat(row,col) = auc_and_p_consensus.p;

                auc_and_p_max = auc_and_p(contains(auc_and_p.Comparison, "max_clo_concord"),:);
                auc_comp_mat_max(row,col) = auc_and_p_max.AUC;
                p_comp_mat_max(row,col) = auc_and_p_max.p;

            elseif comp == "onset_size" 
                auc_and_p = auc_tabs.(sprintf(parc)).(sprintf(comp));
                auc_and_p = auc_and_p(~contains(auc_and_p.Comparison, "_vol"),:);
                auc_and_p = auc_and_p(contains(auc_and_p.Comparison, det_meth),:);
                auc_comp_mat(row,col) = auc_and_p.AUC;
                p_comp_mat(row,col) = auc_and_p.p;

            elseif comp == "resection_size"
                auc_and_p = auc_tabs.(sprintf(parc)).(sprintf(comp));
                auc_and_p = auc_and_p(~contains(auc_and_p.Comparison, "_vol"),:);
                auc_comp_mat(1,col) = auc_and_p.AUC;
                p_comp_mat(1,col) = auc_and_p.p;

            end
        end
    end
    if comp == "onset_resec" | comp == "onset_size" 
        auc_comp_tab.(sprintf(comp)) = array2table(auc_comp_mat,'RowNames', det_meths , 'VariableNames', parcs);
        p_comp_tab.(sprintf(comp)) = array2table(p_comp_mat,'RowNames', det_meths , 'VariableNames', parcs);
    elseif comp == "clo_concord"
        auc_comp_tab.(sprintf("%s_consensus", comp)) = array2table(auc_comp_mat(2:end,:),'RowNames', det_meths(2:end), 'VariableNames', parcs);
        p_comp_tab.(sprintf("%s_consensus", comp)) = array2table(p_comp_mat(2:end,:),'RowNames', det_meths(2:end) , 'VariableNames', parcs);
        auc_comp_tab.(sprintf("%s_max", comp)) = array2table(auc_comp_mat_max(2:end,:),'RowNames', det_meths(2:end), 'VariableNames', parcs);
        p_comp_tab.(sprintf("%s_max", comp)) = array2table(p_comp_mat_max(2:end,:),'RowNames', det_meths(2:end), 'VariableNames', parcs);
    else
        auc_comp_tab.(sprintf(comp)) = array2table(auc_comp_mat(1,:), 'VariableNames', parcs);
        p_comp_tab.(sprintf(comp)) = array2table(p_comp_mat(1,:), 'VariableNames', parcs);
  
    end
end
%% Create heatmaps (for supplementary)
comps = ["onset_resec", "onset_size", "resection_size", "clo_concord_consensus", "clo_concord_max"];

f = figure();
f.Position = [10,10,1000,1000];
cm = [[163 2 52]/255  % red
    [0 118 192]/255 % blue
    [163 2 52]/255];% red
cm = interp1(cm, 1:0.01:size(cm,1));
tiledlayout(5,1)

for comp = comps
    nexttile
    h_auc = heatmap(table2array(auc_comp_tab.(sprintf(comp))));
    clim([0,1])
    title(comp)
   % colorbar off

    if comp ~= "resection_size"
        h_auc.YDisplayLabels = auc_comp_tab.(sprintf(comp)).Properties.RowNames;
    end
        h_auc.XDisplayLabels = auc_comp_tab.(sprintf(comp)).Properties.VariableNames;
        colormap(flipud(cm))
        clim([0,1])
        %colorbar off
         h_auc.CellLabelFormat = '%.3f';
        

end
sgtitle("AUCs across parcellations and onset detection methods")
saveas(f, sprintf("%s/comparing_parcellations_auc.svg", save_loc))

% Plot p-values
f = figure();
f.Position = [10,10,1000,1000];
cm = [repmat([1,1,1],95,1); repmat([1,0,0],4,1)];
tiledlayout(5,1)
for comp = comps
    nexttile
    h_p = heatmap(table2array(p_comp_tab.(sprintf(comp))));
    title(comp)

    if comp ~= "resection_size"
        h_p.YDisplayLabels = p_comp_tab.(sprintf(comp)).Properties.RowNames;
    end
        h_p.XDisplayLabels = p_comp_tab.(sprintf(comp)).Properties.VariableNames;
        colormap(flipud(cm))
        clim([0,1])
        colorbar off
        h_p.CellLabelFormat = '%.3f';

end
sgtitle("p-values across parcellations and onset detection methods")
saveas(f, sprintf("%s/comparing_parcellations_p.svg", save_loc))
