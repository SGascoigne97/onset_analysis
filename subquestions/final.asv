%% Code for results in Incomplete resection of the icEEG seizure onset zone
% is not associated with post-surgical outcomes"
% (S J Gascoigne et al., 2023)

%% Load data and Lausanne atlases
% Set working directory to subquestions folder
cd('subquestions')
addpath(genpath('../sarah_functions'))

load("final_output.mat")

load('../roi_info/ATLAS.mat')
final_output = final_output(final_output.outcome ~= 8,:); % ILAE 8 is used 
% to encode unknown outcome - remove such subjects BEFORE completing downstram analkysis

%% Set parameters for analyses
parc = "roi_120"; % All analyses are based on Lausanne 120 atlas
n_perm = 1000; % All AUC p-values are computed based on permutation test 
               % with 1000 permutations
out_thresh = 2;%2; % ILAE 1-2 is considered as a 'favourable' outcome whilst 
                % ILAE 3+ is considered as an 'unfavourable' outcome
consensus_thresh = 0.5; % Set threshold for inclusion of regions in 
                        % consensus onset, here we will include regions 
                        % present in at least half of the subject's seizures
most_resec_thresh = 0.5; % Set threshold above which subjects are considered 
                         % as having 'most' of their onset resected

det_meths = ["clo", "imprint", "EI"]; % List onset detection methods you 
                                      % are interested in analysing
                                      

% Set parameters for saving figures
save_fig = 1; % Each figure that is created will be saved 
save_folder = '../figures/paper_figures/Figure 3/'; % Location where subfolders will
                                          % be created to store figures and tables

file_type = "svg"; % Figures will be saved as svgs so they can be pulled 
                   % into illustrator for paper figures
    
% for parc = ["chan", "roi_120", "roi_250"]
    [chan_or_roi, n_perm, out_thresh, consensus_thresh, most_resec_thresh,...
        det_meths, save_fig, save_folder, file_type, out_grps] =...
        argument_validation(final_output, parc, n_perm, out_thresh, consensus_thresh,...
        most_resec_thresh, det_meths, save_fig, save_folder, file_type);
    
    %% Additional parameters used throughout analyses (not to be adjused)
    atl_inds = 2*(atlas.scale');
    if contains(chan_or_roi, "roi")
        atl = atlas(atl_inds == str2double(extractAfter(chan_or_roi, "_")),:);
        vols = atl.vol{:};
        names = atl.name{:};
        dists = atl.dists{:};
    else
        atl = NaN;
        vols = NaN;
        names = NaN;
        dists = NaN;
    end
    
    %% Add outcome category column to final_output table
    % Outcome category (binary based on selected outcome threshold)
    final_output.outcome_cat = categorical(final_output.outcome>out_thresh,[0,1], out_grps);
    
    %% 3.1 icEEG Seizure onset regions tend to be resected, but more complete 
    % resections are not associated with more favourable surgical
    %outcomes
    
    % We will compute the proportion of the seizure onset zone (CLO and IOLA)
    % that was subsequently resected then compare across outcome groups
    
    is_onset_resected
    
    %% 3.2 Larger onsets and resections are not associated with surgical outcomes
    
    % Here we will compute the size of seizure onsets (both CLO and
    % automatically captured) using both count of regions and volume (based on
    % controls) and compare across surgical outcome groups
    
    larger_onset
    
    % Here we will compute the size of resections using both count of regions
    % and volume (if using ROIS, based on controls) and compare across surgical
    % outcome groups
    
    larger_resection
% end
    
%% SUPPLEMENTARY
%% S2 Subject Metadata
% Here we will look at various metadata variables to report in the 
% supplementary results. Further, we will investigate if any metadata
% groups have significant differences in surgical outcomes
% have confounding effect on surgical outcomes
metadata

%% S4.2 ALO is not highly concordant with CLO
% We will consider concordance between the consensus onset and the maximum
% consensus between the IOLA onset for any one seizure and the clinically
% labelled onset (CLO)

concordance_with_clo
% concord_tab: subject level concordance between consensus onset and CLO
% and maximum concordance between CLO and any one IOLA onset

% auc_clo_comp: table of AUCs and associated p-values for distinguishing
% surgical outcome groups based on concordance 
    
