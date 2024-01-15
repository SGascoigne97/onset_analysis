function [val_tbl, t] = comp_features(data_tbl, opts)

% input:
%   - data_tbl: full data table
%   - metadata_tbl: full metadata table
%   - optional inputs
%       - window_size: window size for which imprint will be computed
%       - min_sz_count: minimum number of seizures to have been recorded per patient
%       - folder: folder to store markers in

% output
%   - data_tbl: full data table (seizures with no imprint have been
%   removed)
%   - metadata_tbl:
%   - cell_imprint
%   - sz_count_pat: count of focal seizures for all patients meeting
%   inclusion criteria for the minimum number of seizures recorded

    arguments
        data_tbl
        opts.folder = 'onset_calcs'; % folder to store markers in
        opts.window_size (1,1) double {mustBeNumeric} = 1; % Decide window size (in seconds) for which markers are computed
        opts.window_overlap (1,1) double {mustBeNumeric} = 0;
        opts.bands = [1 4; 4 8; 8 13; 13 30; 30 60; 60 100];
    end
    
    %fill in optional arguments
    folder = opts.folder;
    window_size = opts.window_size;
    window_overlap = opts.window_overlap;
    bands = opts.bands;
    
    % Set basefolder to store markers
    subject = data_tbl.patient_id{1};
    basefolder = sprintf('%s/%s', folder, string(subject));
    sampling_rate=data_tbl.segment_fs(1) ;%just using the first, as all sampling was the same in our data after preproc.
    
    linelength_db = LL_db([basefolder '/LL_db/']); %setup folder for all Line Length measures
    linelength_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    linelength_db.paramset_tbl                       % display all currently tracked paramsets
    linelength_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    energy_db = Energy_db([basefolder '/Energy_db']);
    energy_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    energy_db.paramset_tbl                       % display all currently tracked paramsets
    energy_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    %bands: [1 4; 4 8; 8 13; 13 30; 30 60; 60 100]; 
    bandpower_db = BP_db([basefolder '/BP_db']);
    for band = 1:size(bands,1)
        bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',bands(band,:));
    end
    bandpower_db.paramset_tbl                       % display all currently tracked paramsets
    bandpower_db.calc(data_tbl,[]);   
    
    
    %% Import markers (line length, energy, bandpowers)
    calcs_ll = linelength_db.get(1);    % get all calculation outputs in variable. 
    calcs_energy = energy_db.get(1);
    calcs_bp = cell(size(bands, 1),1);
    for band = 1:size(bands,1)
        calcs_bp{band} = bandpower_db.get(band);
    end
    % Extract time windows
    t = calcs_ll.t_wndw;
    
    %% Create val_tbl then log transform all markers
    val_tbl=[calcs_ll.LL_ms calcs_energy.energy_ms];
    for band = 1:size(bands,1)
        val_tbl = [val_tbl, calcs_bp{band}.bp];
    end

    for sz = 1:size(val_tbl,1)
        for feat = 1:size(val_tbl,2)
            val_mat = val_tbl(sz,feat);
            val_tbl(sz,feat) = {log(val_mat{:}+1)};
        end
    end
    val_tbl

end
