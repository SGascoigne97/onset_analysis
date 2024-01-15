function [comp_tab] = compute_resec_comp_measures(onset_tab, opts)
% input:
%   - onset_tab
%   - Optional arguments
%       - det_meth
%       - chan_or_roi
%       - onset_across
%       - onset_acr_thresh
%       - n_perm

% output
%   - comp_with_clo_tab

% Sarah Jane Gascoigne
% 20/07/2023

arguments
    onset_tab
    opts.det_meth = "imprint"
    opts.chan_or_roi = "roi_120"
    opts.onset_across = 1
    opts.onset_acr_thresh = 0.5
    opts.n_perm = 1000
end

% Set optional arguments
det_meth = opts.det_meth;
chan_or_roi = opts.chan_or_roi;
onset_across = opts.onset_across;
onset_acr_thresh = opts.onset_acr_thresh;
n_perm = opts.n_perm;

if det_meth == "clo" & onset_across == 0
    comp_tab = [];
    return
end

% Create tables to store results
col_names = {'Patient_id', 'Sz_count', 'Outcome','Jacc', 'Jacc_z', 'Perc', 'Perc_z', ...
    'Coh', 'Adj_Coh_1', 'Adj_Coh_2'};
if onset_across == 1
     comp_tab = nan(size(onset_tab,1), length(col_names)); % Store double per comparison per subject
else
    comp_tab = cell(size(onset_tab,1), length(col_names)); % Store cell of values per comparison per subject
end
comp_tab = array2table(comp_tab, 'VariableNames',col_names);
        
comp_tab.Patient_id = onset_tab.Patient_id;
comp_tab.Sz_count = cellfun(@length, onset_tab.Segment_ids);
comp_tab.Sz_count(~cellfun(@iscell, onset_tab.Segment_ids)) = 1;
comp_tab.Outcome = onset_tab.outcome;

for pat = 1:size(onset_tab,1)
    clear jacc jacc_z perc perc_z coh coh_z
    % Extract automatically detected onset and CLO
    pat_onset = onset_tab(pat,:);
    onset = pat_onset.(sprintf("%s_%s",det_meth, chan_or_roi)){:};
    resec = pat_onset.(sprintf("resected_%s", chan_or_roi)){:};
    if sum(resec) == 0
        fprintf("%s No resection included \n", pat_onset.Patient_id{:})
        continue
    end
        
    if onset_across == 1 % compute one onset across all seizures
        onset_acr = sum(onset,2)/size(onset,2) >= onset_acr_thresh;

        if sum(onset_across) == 0
            fprintf("%s No onset across found \n", pat_onset.Patient_id{:})
            continue
        end
    
        jacc = jaccard(logical(resec), logical(onset_acr));
        perc = sum(resec + onset_acr == 2)/sum(onset_acr);
        coh = cohensKappa(logical(resec),logical(onset_acr));

        onset_acr_adj = double(onset_acr);
        resec_adj = double(resec);
        onset_acr_adj((2*resec + onset_acr) ==2) = NaN;
        resec_adj((2*resec + onset_acr) ==2) = NaN;
        adj_coh_1 = cohensKappa(resec_adj,onset_acr_adj);
    
        resec_adj = double(resec);
        resec_adj((2*resec + onset_acr) ==2) = 0;
        adj_coh_2 = cohensKappa(logical(resec_adj),logical(onset_acr)); 

        % Compute comparison measures
        jacc_perm = zeros(1,n_perm);
        perc_perm = jacc_perm; 

        for perm = 1:n_perm
            rng(perm)
            perm_onset = onset_acr(randperm(length(onset_acr)));
            rng(perm+1)
            perm_resec = resec(randperm(length(resec)));
            jacc_perm(perm) = jaccard(double(perm_onset),double(perm_resec));
            perc_perm(perm) = sum(perm_resec + perm_onset == 2)/sum(perm_onset);
        end

        jacc_z = (jacc - mean(jacc_perm))/std(jacc_perm);
        perc_z = (perc - mean(perc_perm))/std(perc_perm);

        comp_tab(pat, 4:10) = table(jacc, jacc_z, perc,...
        perc_z, coh, adj_coh_1, adj_coh_2);
    else % if onset_across == 0 (do not compute onset across seizures)
        jacc = nan(size(onset,2),1);
        perc = jacc;
        coh = jacc; 
        adj_coh = jacc;
        jacc_z = jacc;
        perc_z = jacc;
        coh_z = jacc; 
        adj_coh_1 = jacc;
        adj_coh_2 = jacc;

        for sz = 1:size(onset,2)
            sz_onset = onset(:,sz);
            if sum(sz_onset) == 0
                fprintf("%s seizure %d No seizure onset found \n", pat_onset.Patient_id{:}, sz)
                continue
            end
            jacc(sz) = jaccard(logical(resec), logical(sz_onset));
            perc(sz) = sum(resec + sz_onset == 2)/sum(sz_onset);
            coh(sz) = cohensKappa(logical(resec),logical(sz_onset));
            sz_onset_adj = double(sz_onset);
            resec_adj = double(resec);
            sz_onset_adj((2*resec + sz_onset) ==2) = NaN;
            resec_adj((2*resec + sz_onset) ==2) = NaN;
            adj_coh_1(sz) = cohensKappa(resec_adj,sz_onset_adj);
            resec_adj = double(resec);
            resec_adj((2*resec + sz_onset) ==2) = 0;
            adj_coh_2(sz) = cohensKappa(resec_adj,sz_onset); 

            % Compute comparison measures
            jacc_perm = zeros(1,n_perm);
            perc_perm = jacc_perm; 
            for perm = 1:n_perm
                rng(perm)
                perm_onset = sz_onset(randperm(length(sz_onset)));
                rng(perm+1)
                perm_resec = resec(randperm(length(resec)));
                jacc_perm(perm) = jaccard(double(perm_onset),double(perm_resec));
                perc_perm(perm) = sum(perm_resec + perm_onset == 2)/sum(perm_onset);
            end

            jacc_z(sz) = (jacc(sz) - mean(jacc_perm))/std(jacc_perm);
            perc_z(sz) = (perc(sz) - mean(perc_perm))/std(perc_perm);
        end
        comp_tab(pat, 4:10) = [{jacc}, {jacc_z}, {perc}, {perc_z}, {coh}, {adj_coh_1}, {adj_coh_2}];
    end
    
end
if iscell(comp_tab(1,:).Jacc_z)
    comp_tab = comp_tab(~cellfun(@isempty,comp_tab.Jacc_z),:);
else
    comp_tab = comp_tab(~isnan(comp_tab.Jacc_z),:);
end