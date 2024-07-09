% This code must be run on an individual patient basis
%%
load('roi_info/ATLAS_.mat')
load('roi_info/retains.mat')
load('roi_info/atlasinfo.mat')
load('roi_info/ATLAS.mat')

%% Select patient to compute values for
patient = "U10";
patient_old_lab = string(pat_id(pat_id(:,2) == patient,1));
channel_folder = '../../channels_ROI/';
channel_info = load(sprintf('%s%s/channels.mat',channel_folder,patient_old_lab));
pat_data = data(string(data.patient_id) == patient,:);
incl_chan = ismember(string(channel_info.channels.name), string(pat_data.segment_channel_labels{1,1}));
ROI_incl_channels = channel_info.channels.ROIids(incl_chan,:);

unq_roi = unique(ROI_incl_channels(:,3));

% Extract onset information
pat_onset = onset_output(onset_output.Patient_id ==patient,:);
labelled_onset = unq_roi.*cell2mat(pat_onset.Labelled_onset);
labelled_onset = labelled_onset(labelled_onset~=0);

%% Plot ROIs on crystal brain and highlight labelled onset channels 
load('roi_info/lhSurf.mat') 
load('roi_info/rhSurf.mat')

h = figure(1);
%for i = 1
    %h(i)=subplot(1,1,i);
    
    Hl = patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    camlight('headlight','infinite')
    hold on
    Hv = patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    camlight('headlight','infinite')
    axis equal

%end

plot3(xyz125(unq_roi,1),xyz125(unq_roi,2),xyz125(unq_roi,3), 'o', ...
    'MarkerSize', 12, 'MarkerEdgeColor', 'k')
plot3(xyz125(labelled_onset,1),xyz125(labelled_onset,2),...
    xyz125(labelled_onset,3), 'o', ...
    'MarkerSize', 12, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r')
title(sprintf('Patient %s ROIs with labelled onset (based on reports)',patient))
hold off


%% Plot ROIs on crystal brain andhighlight onset channels based on each onset detection method
sz = 3;
h = figure(sz);
sgtitle(sprintf('Patient %s seizure %s onset channels',patient, string(sz)))
for i = 1:4
    h(i)=subplot(2,2,i);
    Hl = patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    camlight('headlight','infinite')
    hold on
    Hv = patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    camlight('headlight','infinite')
    axis equal

    plot3(xyz125(unq_roi,1),xyz125(unq_roi,2),xyz125(unq_roi,3), 'o', ...
            'MarkerSize', 6, 'MarkerEdgeColor', 'k')

    if i == 1
        plot3(xyz125(labelled_onset,1),xyz125(labelled_onset,2),...
            xyz125(labelled_onset,3), 'o', ...
            'MarkerSize', 6, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r')
        title("Labelled onset")
        hold off
    else
        onset_by_sz = onset_output(onset_output.Patient_id == patient,3+i);
        onset_by_sz = cell2mat(onset_by_sz{1,1});
        onset = onset_by_sz(:,sz);
        onset_roi = unq_roi.*onset;
        onset_roi = onset_roi(onset_roi~=0);
        plot3(xyz125(onset_roi,1),xyz125(onset_roi,2),...
            xyz125(onset_roi,3), 'o', ...
            'MarkerSize', 6, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'r')
        title(sprintf("%s", onset_output.Properties.VariableNames{3+i}) )
        hold off
    end

end

%% Peter's code
MasterChannelTable = channel_info.channels; 
%MasterChannelTable.ROIID_inAtlas=double(MasterChannelTable.ROIID_inAtlas);
MasterChannelTable.ROI_inAtlas = double(MasterChannelTable.ROIids(:,3));
%MasterChannelTable.ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas==0)=nan;
count=1;
ROIID_inAtlas=nan(size(MasterChannelTable.ROIID_inAtlas,1),4);
for i=1:size(scale36retain)
    if scale36retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,1)==i,1)=count;
        count=count+1;
    end
end
count=1;
for i=1:size(scale60retain)
    if scale60retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,2)==i,2)=count;
        count=count+1;
    end
end
count=1;
for i=1:size(scale125retain)
    if scale125retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,3)==i,3)=count;
        count=count+1;
    end
end
count=1;
for i=1:size(scale250retain)
    if scale250retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,4)==i,4)=count;
        count=count+1;
    end
end
MasterChannelTable.ROIID_inAtlas=ROIID_inAtlas;
clear i ROIID_inAtlas count

%%

ATLAS.name36=ATLAS.name36(scale36retain==1);
ATLAS.dists36=dists36(scale36retain==1,scale36retain==1);
ATLAS.xyz36=xyz36(scale36retain==1,:);
ATLAS.vol36=vol36(scale36retain==1);
clear scale36retain dists36 xyz36 vol36

ATLAS.name60=ATLAS.name60(scale60retain==1);
ATLAS.dists60=dists60(scale60retain==1,scale60retain==1);
ATLAS.xyz60=xyz60(scale60retain==1,:);
ATLAS.vol60=vol60(scale60retain==1);
clear scale60retain dists60 xyz60 vol60

ATLAS.name125=ATLAS.name125(scale125retain==1);
ATLAS.dists125=dists125(scale125retain==1,scale125retain==1);
ATLAS.xyz125=xyz125(scale125retain==1,:);
ATLAS.vol125=vol125(scale125retain==1);
clear scale125retain dists125 xyz125 vol125

ATLAS.name250=ATLAS.name250(scale250retain==1);
ATLAS.dists250=dists250(scale250retain==1,scale250retain==1);
ATLAS.xyz250=xyz250(scale250retain==1,:);
ATLAS.vol250=vol250(scale250retain==1);
clear scale250retain dists250 xyz250 vol250
%%
scale=[36;60;125;250];
name={ATLAS.name36;ATLAS.name60;ATLAS.name125;ATLAS.name250};
dists={ATLAS.dists36;ATLAS.dists60;ATLAS.dists125;ATLAS.dists250;};
vol={ATLAS.vol36;ATLAS.vol60;ATLAS.vol125;ATLAS.vol250;};
xyz={ATLAS.xyz36;ATLAS.xyz60;ATLAS.xyz125;ATLAS.xyz250;};
T=table(scale,name,dists,vol,xyz);
clear scale name dists vol xyz

%%

nbroi=size(T.vol{parc},1);

