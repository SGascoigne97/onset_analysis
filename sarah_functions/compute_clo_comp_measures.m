function [comp_with_clo_tab] = compute_clo_comp_measures(onset_tab, opts)
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
% 19/07/2023

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

% Create tables to store results
col_names = {'Patient_id', 'Sz_count', 'Outcome','Jacc', 'Jacc_z', 'Perc', 'Perc_z', ...
    'Coh', 'Coh_z'};
if onset_across == 1
     comp_with_clo_tab = zeros(size(onset_tab,1), length(col_names)); % Store double per comparison per patient
else
    comp_with_clo_tab = cell(size(onset_tab,1), length(col_names)); % Store cell of values per comparison per patient
end
comp_with_clo_tab = array2table(comp_with_clo_tab, 'VariableNames',col_names);
        
comp_with_clo_tab.Patient_id = onset_tab.Patient_id;
comp_with_clo_tab.Sz_count = cellfun(@length, onset_tab.Segment_ids);
comp_with_clo_tab.Sz_count(~cellfun(@iscell, onset_tab.Segment_ids)) = 1;
comp_with_clo_tab.Outcome = onset_tab.outcome;

for pat = 1:size(onset_tab,1)
    clear jacc jacc_z perc perc_z coh coh_z
    % Extract automatically detected onset and CLO
    pat_onset = onset_tab(pat,:);
    onset = pat_onset.(sprintf("%s_%s",det_meth, chan_or_roi)){:};
    clo = pat_onset.(sprintf("clo_%s", chan_or_roi)){:};
    if sum(clo) == 0
        fprintf("%s No CLO included \n", pat_onset.Patient_id{:})
        continue
    end
        
    if onset_across == 1 % compute one onset across all seizures
        onset_acr = sum(onset,2)/size(onset,2) >= onset_acr_thresh;

        if sum(onset_across) == 0
            fprintf("%s No onset across found \n", pat_onset.Patient_id{:})
            continue
        end
    
        jacc = jaccard(logical(clo), logical(onset_acr));
        perc = sum(clo + onset_acr == 2)/sum(clo);
        coh = cohensKappa(logical(clo),logical(onset_acr));

        % Compute comparison measures
        jacc_perm = zeros(1,n_perm);
        perc_perm = jacc_perm; 
        coh_perm = jacc_perm;
        for perm = 1:n_perm
            rng(perm)
            perm_onset = onset_acr(randperm(length(onset_acr)));
            rng(perm+1)
            perm_clo = clo(randperm(length(clo)));
            jacc_perm(perm) = jaccard(double(perm_onset),double(perm_clo));
            perc_perm(perm) = sum(perm_clo + perm_onset == 2)/sum(perm_clo);
            coh_perm(perm) = cohensKappa(logical(perm_onset),logical(perm_clo));
        end

        jacc_z = (jacc - mean(jacc_perm))/std(jacc_perm);
        perc_z = (perc - mean(perc_perm))/std(perc_perm);
        coh_z = (coh- mean(coh_perm))/std(coh_perm);

        comp_with_clo_tab(pat, 4:9) = table(jacc, jacc_z, perc,...
        perc_z, coh, coh_z);
    else % if onset_across == 0 (do not compute onset across seizures)
        jacc = nan(size(onset,2),1);
        perc = jacc;
        coh = jacc; 
        jacc_z = jacc;
        perc_z = jacc;
        coh_z = jacc; 
        for sz = 1:size(onset,2)
            sz_onset = onset(:,sz);
            if sum(sz_onset) == 0
                fprintf("%s seizure %d No seizure onset found \n", pat_onset.Patient_id{:}, sz)
                continue
            end
            jacc(sz) = jaccard(logical(clo), logical(sz_onset));
            perc(sz) = sum(clo + sz_onset == 2)/sum(clo);
            coh(sz) = cohensKappa(logical(clo),logical(sz_onset));

            % Compute comparison measures
            jacc_perm = zeros(1,n_perm);
            perc_perm = jacc_perm; 
            coh_perm = jacc_perm;
            for perm = 1:n_perm
                rng(perm)
                perm_onset = sz_onset(randperm(length(sz_onset)));
                rng(perm+1)
                perm_clo = clo(randperm(length(clo)));
                jacc_perm(perm) = jaccard(double(perm_onset),double(perm_clo));
                perc_perm(perm) = sum(perm_clo + perm_onset == 2)/sum(perm_clo);
                coh_perm(perm) = cohensKappa(logical(perm_onset),logical(perm_clo));
            end

            jacc_z(sz) = (jacc(sz) - mean(jacc_perm))/std(jacc_perm);
            perc_z(sz) = (perc(sz) - mean(perc_perm))/std(perc_perm);
            coh_z(sz) = (coh(sz)- mean(coh_perm))/std(coh_perm);
        end
        comp_with_clo_tab(pat, 4:9) = [{jacc}, {jacc_z}, {perc},...
        {perc_z}, {coh}, {coh_z}];

    end
end