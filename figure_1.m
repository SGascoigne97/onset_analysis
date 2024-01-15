%% Code for creating figure one
% Direct to location of saved subject data
data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];
basefolder = 'example_onset_detection';

addpath('lib_biomarkers/')
addpath('lib_dataflow/')
addpath(genpath('lib/'))

%% Load data and Lausanne atlases
%final_output_struct = load('tables/final_output.mat'); % Here any subjects without follow-up have been removed
%final_output = final_output_struct.final_output;
% clear final_output_struct
addpath(genpath('sarah_functions'))

load('roi_info/ATLAS.mat')
%%
% Select a subject to create figure for
pat = 24;
pat_onset = final_output(pat,:);
pat_id = string(pat_onset.Patient_id);

atlas = 2; %1=scale36,2=60,3=125,4=250
% Choose atlas
% Choose threshold for allowing propagated activity (seconds)
min_sz = 1;
min_sz_dur = 9; % minimum seizure duration for inclusion
chan_to_roi_thresh_type = "count";
chan_to_roi_thresh = 1; % One channel in region is sufficient to include region
window_size = 1; % size of sliding window (in seconds)
window_overlap = 7/8; % Overlap between imprint windows
det = 7; % count of windows beyond first detected onset which will also be considered as onset
rec_type = "sec"; % type of threshold used to validate activity is ictal activity (minimum duration) - can use "sec" or "prop"
rec_thresh = 9; % number of seconds for which activity must persist to be labelled as ictal
mad_thresh = 3;
ict_buffer = 10;

onset_calc_loc = "imprint_ons"; % Specify folder to store imprint values in
% Need a new folder if using a different subset of the data/different data
% as it will load previous save if folder is not empty

% We will show a subsection of EEG channels (25) for clearer visualisations
% We must include any onset channels so imprint onset across seizures has
% no empty columns
incl_channels = sum(final_output(pat,:).imprint_chan{:}, 2) > 0;
remaining_channels = find(~incl_channels);

rng(1)
random_channels = remaining_channels(randsample(length(remaining_channels), 25 - sum(incl_channels)));
incl_channels(random_channels) = 1;
incl_chan_id = find(incl_channels);

% Organise channels into order of regions
roi_names = pat_onset.roi_names_120{:}; 
% List channel ids in region order
chan_ind = [];
roi_ind = [];
for roi = 1:length(roi_names)
    chan_ind = [chan_ind, find(final_output.chan_2_roi_matrix_120{pat,1}(roi,:))];
    roi_ind = [roi_ind, repmat(roi, 1, length(find(final_output.chan_2_roi_matrix_120{pat,1}(roi,:))))];
end
roi_tbl = table(chan_ind', roi_ind', 'VariableNames', ["channel_id", "roi_id"]);

% roi_tbl is my channel order reference table
roi_tbl_included = roi_tbl(ismember(roi_tbl.channel_id, incl_chan_id),:);
unq_roi = unique(roi_tbl_included.roi_id);

incl_chan_ordered = flipud(roi_tbl_included.channel_id);

% Compute imprint
pat_data = load(sprintf('%s/%s.mat', data_location, pat_id));
pat_data = pat_data.data_export;

pat_onset = final_output(final_output.Patient_id == pat_id,:);
pat_data = pat_data(ismember(pat_data.segment_id, pat_onset.Segment_ids{1,1}),:);
pat_meta = pat_data(:, 1:end);

[pat_data, cell_imprint,  sz_count_pat] = calc_imprint_mahal(pat_data,...
    "window_size", window_size, "min_sz_count", min_sz,...
    "folder", onset_calc_loc, "mad_thresh", mad_thresh,...
    "rec_thresh", rec_thresh, "window_overlap", window_overlap);

cell_imprint = cell_imprint(find(ismember(cell_imprint.segment_id, pat_onset.Segment_ids{1,1})),:);

%%
% xlims = nan(size(pat_data,1), 2);
% xlims([1,9],:) = [720, 1200; 720, 1200];

for sz = [4, 19,  33]
    f = figure('Visible', 'on');
    f.Position = [100,100, 1600, 1000];
    subplot(1,9,1:8)
    eeg_dat = pat_data.segment_data{sz,1};
    offset = round((max(max(eeg_dat))-min(min(eeg_dat)))/3);
    imprint = cell_imprint.cell_imprint{sz,1};
    imprint_bin = [zeros(size(imprint,1),110*8), imprint];
    imprint_bin = repelem(imprint_bin(incl_chan_ordered,:),offset,1);
    for chan = 1:size(eeg_dat,1)
        top = 900 + (chan-1)*offset;
        bottom = chan*offset;
        imprint_bin(top:bottom,:) = 0;
    end
    imprint_bin = [imprint_bin; zeros(round(offset/2),size(imprint_bin,2))];
    imprint_bin = imprint_bin(round(offset/2):end,:);

    % Extract segments of EEG which are included in imprint
    % Add preictal and ictal segments to imprint
    imprint_full = [zeros(size(imprint,1),110*8), imprint, zeros(size(imprint,1),120*8)];
    imprint_eeg = repelem(imprint_full, 1, round(size(eeg_dat,2)/size(imprint_full,2)));
    % Add in zeros at each side so the sizes are equivalent 
    dif = size(eeg_dat,2) - size(imprint_eeg,2);
    if dif > 0
        imprint_eeg = [zeros(size(imprint_eeg,1), round(dif/2)), imprint_eeg, zeros(size(imprint_eeg,1), round(dif/2)-1)];
    elseif dif < 0
        dif = abs(dif);
        imprint_eeg = imprint_eeg(:, round(dif/2):end-round(dif/2));
    end

    imprint_eeg(imprint_eeg == 0) = NaN;
    imprint_eeg_dat = eeg_dat.*imprint_eeg;

    if ~isempty(find(sum(imprint_bin,1)>0, 1, 'first'))
        onset_heatmap = zeros(size(imprint_bin,1), size(imprint_bin,2));
        onset_heatmap(:,find(sum(imprint_bin,1)>0, 1, 'first') +[0:7]) = ...
            imprint_bin(:,find(sum(imprint_bin,1)>0, 1, 'first') +[0:7]);
    else 
%         continue
    end

    hold on
    imagesc(imprint_bin, 'AlphaData', 0.5)
    imagesc(2*onset_heatmap, "AlphaData", 0.7)
    vis_eeg(eeg_dat(incl_chan_ordered,:),...
        512/8, "ChannelNames", pat_data.segment_channel_labels{1,1}(incl_chan_ordered),...
        "PlotNewFig", false, "Color", [0.2,0.2,0.2], "Offset", offset);
    vis_eeg(imprint_eeg_dat(incl_chan_ordered,:),...
        512/8, "ChannelNames", pat_data.segment_channel_labels{1,1}(incl_chan_ordered),...
        "PlotNewFig", false, "Color", [0,0,0.7], "Offset", offset);
    hold off
    colormap([1,1,1;0,0.2,1;0,0.5,0])
    clim([0,2])
    xline(120*8, LineWidth=2, Color="red")
    xline(find(sum(imprint_bin),1, 'first'), LineWidth=2, Color="blue")
    xline((120+pat_data.duration(sz))*8, LineWidth=2, Color="red")
    set(gca, "XTick", 0:10*8:size(eeg_dat,2)/8/8, "XTickLabel", -120+(0:10:size(eeg_dat,2)/8^3))
    xlim([720, 1200])
    ylim([-offset, offset*35.5])
    
    subplot(1,9,9)
    imagesc(flipud(2*final_output(pat,:).imprint_chan{:}(incl_chan_ordered,sz)), 'AlphaData', 0.5)
    set(gca, "YTick", 1:length(incl_chan_ordered), "YTickLabel", flipud(string(final_output(pat,:).channel_names{:}(incl_chan_ordered))))
    clim([0,2])
    

    sgtitle(sprintf("%s seizure %d (%s)", string(pat_data(sz,:).patient_id), sz, pat_data.ilae_sz_type(sz)))
    saveas(f,sprintf("figures/paper_figures/Figure 1/%s_sz%d.png",  string(pat_data(1,:).patient_id), sz))

end
 

%% 
%% Panel C: Channel-wise imprint onset 
onsets = final_output(pat,:).imprint_chan{:}(incl_chan_ordered,:);
f = figure();
heatmap(flipud(onsets), 'CellLabelColor','none')
colormap([0.8,0.8,0.8;1,0.5,0.5])
colorbar off
saveas(f,sprintf("figures/paper_figures/Figure 1/%s_all_chan_onsets.svg", string(pat_onset.Patient_id)))

%% Panel D: ROI-120 imprint onset (regions included in channels from Panels A-C)
f = figure();
onsets = final_output(pat,:).imprint_roi_120{:}(unq_roi,:);
subplot(1,5,1:4)
heatmap(onsets, 'CellLabelColor','none')
colorbar off
subplot(1,5,5)
heatmap(double(mean(final_output(pat,:).imprint_roi_120{:}(unq_roi,:),2)>=0.5), 'CellLabelColor','none')
colormap([0.8,0.8,0.8;1,0.5,0.5])
colorbar off
saveas(f,sprintf("figures/paper_figures/Figure 1/%s_all_roi_onsets.svg", string(pat_onset.Patient_id)))


%% Panel E: Onset regions for two example seizures highlighted on brain
roi_names = pat_onset.roi_names_120{:};
roi_names = strrep(roi_names, 'r.', 'ctx-rh-');
roi_names = strrep(roi_names, 'l.', 'ctx-lh-');
cm = [0.8,0.8,0.8;1,0.5,0.5];

plotBrain(roi_names, double(mean(final_output(pat,:).imprint_roi_120{:},2)>=0.5),...
    cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1], 'savePath', 'figures/paper_figures/Figure 1/consensus')

cm = [0.8,0.8,0.8;0,0,0];
plotBrain(roi_names, final_output(pat,:).resected_roi_120{:},...
    cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1], 'savePath', 'figures/paper_figures/Figure 1/resected')

cm = [0.8,0.8,0.8;0.8549,0.4392,.8392];
plotBrain(roi_names, final_output(pat,:).clo_roi_120{:},...
    cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1], 'savePath', 'figures/paper_figures/Figure 1/clo')

% plotBrain(roi_names, pat_onset.imprint_roi_120{:}(:,3),cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1])%, 'savePath', 'figures/paper_figures/sz3')
% plotBrain(roi_names, pat_onset.imprint_roi_120{:}(:,7),cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1], 'savePath', 'figures/paper_figures/sz7')
