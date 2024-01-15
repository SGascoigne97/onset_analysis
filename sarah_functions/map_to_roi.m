function [mapping_ch2roi,ch_names,roi_names, channel_details_used]=map_to_roi(json_data, pat_data, atl)
%% channel details
channel_details=json_data(1,:).channel_details;

%% filter by the actual channels used
%channels_used=json_data(1,:).pre_eeg_channels;
channels_used = pat_data(1,:).segment_channel_labels{:};

tbl_ids=zeros(length(channels_used),1);
for i=1:length(channels_used)
    id = find(strrep(string(channel_details.chan_name),' ', '')...
        ==strrep(channels_used(i),' ',''));
    if ~isempty(id)
        tbl_ids(i)=id;
    else
        tbl_ids(i)=NaN;
    end
end


channel_details_used=channel_details( tbl_ids(~isnan(tbl_ids)),:);%this is now in the right order as the EEG data in terms of channels
ch_names=channel_details_used.chan_name;
%% create a matrix of dimensions n_chan x nROI, and a associated channel names list & ROI names list

allROIs=cat(2,channel_details_used.ROIname{:});
roi_names=unique(allROIs(atl,:))';

mapping_ch2roi=zeros(length(ch_names),length(roi_names));
for i=1:length(roi_names)
    chids=find(string(allROIs(atl,:))==roi_names{i});
    mapping_ch2roi(chids,i)=1;
end

mapping_ch2roi=mapping_ch2roi';
