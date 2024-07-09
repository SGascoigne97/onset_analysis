function [comp_tab] = compute_resec_comp(onset_tab,atlas,vols, opts)
% input:
%   - onset_tab
%   - Optional arguments
%       - det_meth
%       - chan_or_roi
%       - onset_across
%       - onset_acr_thresh
%       - n_perm

% output
%   - comp_tab

% Sarah Jane Gascoigne
% 10/08/2023

% Computes percentage resected for binary arrays and volumes 

arguments
    onset_tab
    atlas
    vols
    opts.det_meth = "imprint"
    opts.chan_or_roi = "roi_120" % comparing volumes can only be done on ROI-wise onset/resection
    opts.onset_across = 1
    opts.consensus_thresh = 0.5
end

% Set optional arguments
det_meth = opts.det_meth;
chan_or_roi = opts.chan_or_roi;
onset_across = opts.onset_across;
consensus_thresh = opts.consensus_thresh;

if det_meth == "clo" & onset_across == 0
    comp_tab = [];
    return
end

% Create tables to store results
if chan_or_roi == "chan"
    col_names = {'Patient_id', 'Sz_count', 'Outcome', 'Perc'}; 
else
    col_names = {'Patient_id', 'Sz_count', 'Outcome', 'Perc', 'Perc_vol'}; 
    names = atlas.name{:};
end

if onset_across == 1
     comp_tab = nan(size(onset_tab,1), length(col_names)); % Store double per comparison per subject
else
    comp_tab = cell(size(onset_tab,1), length(col_names)); % Store cell of values per comparison per subject
end
comp_tab = array2table(comp_tab, 'VariableNames',col_names);
        
comp_tab.Patient_id = onset_tab.Patient_id;
chan_and_sz_count = cell2mat(cellfun(@size, onset_tab.imprint_chan, 'UniformOutput', false));
comp_tab.Sz_count = chan_and_sz_count(:,2);
comp_tab.Outcome = onset_tab.outcome_cat;

for pat = 1:size(onset_tab,1)
    clear perc perc_vol names_pat 

    % Extract automatically detected onset and CLO
    pat_onset = onset_tab(pat,:);
    onset = pat_onset.(sprintf("%s_%s",det_meth, chan_or_roi)){:};
    resec = pat_onset.(sprintf("resected_%s", chan_or_roi)){:};

    if sum(resec) == 0
        fprintf("%s No resection included \n", pat_onset.Patient_id{:})
        continue
    end

    if det_meth == "clo" & sum(onset,1) == 0
        fprintf("%s No CLO included \n", pat_onset.Patient_id{:})
        continue
    end
    
    if chan_or_roi~="chan"
        regions = pat_onset.(sprintf("roi_names_%s", extractAfter(chan_or_roi, "roi_"))){:};
        % Find index of regions within atlas
        reg_ind = nan(length(regions),1);
        for reg = 1:length(regions)
            reg_ind(reg) = find(contains(names, regions(reg)));
        end
    end

    % Compute consensus onset (if onset_across == 1)
    if onset_across == 1 % compute one onset across all seizures
        if det_meth ~= "clo"
            onset_acr = mean(onset(:,sum(onset,1)>0),2)>=consensus_thresh;
            %onset_acr = mean(onset,2) >= consensus_thresh;
        else
            onset_acr = onset;
        end

        if sum(onset_across) == 0
            fprintf("%s No onset across found \n", pat_onset.Patient_id{:})
            continue
        end

        perc = sum((onset_acr+resec) == 2)/sum(onset_acr);

        if chan_or_roi ~= "chan"
            vol_ons = vols(reg_ind(logical(onset_acr)));
            vol_ons_resec = vols(reg_ind((onset_acr+resec)==2));
            perc_vol = sum(vol_ons_resec)/sum(vol_ons); 
    
            comp_tab(pat, 4:5) = table(perc, perc_vol); % Removed z-score as distributions were skewed by large spike at 0
        else
             comp_tab(pat, 4) = table(perc);
        end

    else % if onset_across == 0 (do not compute onset across seizures)
        perc = nan(size(onset,2),1);
        if chan_or_roi ~= "chan"
            perc_vol = perc;
        end

        for sz = 1:size(onset,2)
            sz_onset = onset(:,sz);
            if sum(sz_onset) == 0
                fprintf("%s seizure %d No seizure onset found \n", pat_onset.Patient_id{:}, sz)
                continue
            end
            perc(sz) = sum((sz_onset+resec)==2)/sum(sz_onset); 

            if chan_or_roi ~= "chan"
                vol_ons = vols(reg_ind(logical(onset_acr)));
                vol_ons_resec = vols(reg_ind((onset_acr+resec)==2));
                perc_vol(sz) = sum(vol_ons_resec)/sum(vol_ons); 
            end
        
        % Removed z-score as distributions were skewed by large spike at 0
        end
        if chan_or_roi ~= "chan"
            comp_tab(pat, 4:5) = [{perc}, {perc_vol}];   
        else
            comp_tab(pat, 4) = {perc}; 
        end
    end
end

if iscell(comp_tab(1,:).Perc)
    comp_tab = comp_tab(~cellfun(@isempty,comp_tab.Perc),:);
else
    comp_tab = comp_tab(~isnan(comp_tab.Perc),:);
end