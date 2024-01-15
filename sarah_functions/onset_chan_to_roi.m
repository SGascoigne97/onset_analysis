function [onset_output] = onset_chan_to_roi(pat_onset, json_data, onset_output, opts)
% Compute average EEG across regions of interest for each patient

% input:
%   - data_tbl: data table for ONE patient
%   - json_data: json data for same patient
%   - Onset_output: output table with onset detected on a channel level for
%   this patient
%   - optional inputs
%       - method: method of automatic detection used ("imprint", "EI", or
%       "PLHG")

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        pat_onset
        json_data
        onset_output
        opts.method (1,1) string {mustBeMember(opts.method, ["imprint", "EI", "PLHG"])} = "imprint" % state method of automatic onset detection
        opts.thresh (1,1) double = nan % define as value within 0 and 1 (i.e. 0.25 = 25% of channels) or NaN if only one channel is required
        opts.atlas (1,1) double {mustBeMember(opts.atlas, [36, 60, 125, 250])} = 60

    end
    
    %fill in optional arguments
    method = opts.method;
    thresh = opts.thresh;
    atlas = opts.atlas;

    if method == "imprint"
        onset = onset_output.imprint_chan{:};
    elseif method == "EI"
        onset = onset_output.EI_chan{:};
    elseif method == "PLHG"
        onset = onset_output.PLHG_chan{:};
    end
  
    % Add empty columns to store output
    onset_output = [onset_output cell(size(onset_output,1),1) cell(size(onset_output,1),1)];
    onset_output.Properties.VariableNames([-1,0]+size(onset_output,2)) = [sprintf("%s_roi", method), sprintf("%s_roi_names", method)];
    %%
    channel_details = json_data(1).channel_details;
    recorded_channels = channel_details(...
        ismember(strrep(channel_details.chan_name,' ', ''),...
        strrep(onset_output.channel_names{:},' ','')),:);
    recorded_roi_all_atlas = cat(2,recorded_channels.ROIname{:});
    
     atl_id =  [36; 60; 125; 250];
     col_id = [1;2;3;4];
     atl_tab = table(atl_id, col_id);
     % Find the row to pull atlas ID from    
     atl_row = table2array(atl_tab(atl_tab.atl_id == atlas,"col_id"));
     
     % Extract the name of included regions
     incl_roi = recorded_roi_all_atlas(atl_row,:)';
     unq_roi = unique(incl_roi, 'stable');
    
    onset_roi = zeros(length(unq_roi), size(pat_onset,1));
    nam_mat = strings(length(unq_roi), size(pat_onset,1));
    for sz = 1:size(onset,2)
        [onset_roi(:,sz), names] = chan_to_roi_crit(onset(:,sz), incl_roi, unq_roi, "threshold", thresh);
        nam_mat(1:length(names),sz) = names;

    end

    onset_output(:,size(onset_output,2)-1) =  mat2cell(onset_roi,size(onset_roi,1),size(onset_roi,2));
    onset_output(:,size(onset_output,2)) =  mat2cell(nam_mat,size(onset_roi,1),size(onset_roi,2));

%     recorded_channels = channel_details(ismember(strrep(channel_details.chan_name, ' ', ''),...
%         strrep(pat_onset.segment_channel_labels{1},' ', '')),:);
%     recorded_roi_all_atlas = cat(2,recorded_channels.ROIname{:});
%     incl_roi = recorded_roi_all_atlas(3,:)';
%     
%     pat_onset = onset{1,1};
%     unq_roi = onset_output.ROI_ids{1,1};
%    
%     onset_roi = zeros(length(unq_roi), size(pat_onset,1));
% 
%     for grp = 1:length(unq_roi)
%         onset_roi(grp,:) = sum(pat_onset(string(incl_roi) == string(unq_roi(grp)),:),1);
%     end
   
end