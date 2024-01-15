function [onset_output] = compute_onset(pat_data, pat_imprints, opts)

% input:
%   - pat_data: % table with all seizures for subject
%   - pat_imprints: table of imprints for all seizures included in data table
%   - Optional arguments
%       - opts.wdw_size: Integer denoting window size for considering activity 
                       % as included in the onset (default 2 to reduce risk 
                       % of inclusion of spurious activity)
%       - opts.det: Integer denoting the number of seconds following onset 
                        % detection where activity will be included

% output
%   - onset_output: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
% pat_data

% Sarah Jane Gascoigne & Yujiang Wang
% 25/05/2023

    arguments
        pat_data 
        pat_imprints 
        opts.wdw_sz (1,1) double = 2 
        opts.det (1,1) double {mustBeNumeric} = 0
    end

    % Set optional arguments
    wdw_sz = opts.wdw_sz;
    det = opts.det;
    
    % We could include an optional argument to determine which onset detection
    % methods we will use
    % Additional arguents to change parameters in detection methods (e.g.,
    % lambda and v in EI)

    colNames = {'imprint_chan','EI_chan', 'PLHG_chan', 'when_onset'};
    onset_output = cell2table(cell(1, length(colNames)), 'VariableNames', colNames);
    
    %% Iterate through seizures and compute onset for each method
    sz_count = size(pat_imprints,1); % Number of recorded seizures
    n_chan =  size(pat_data.segment_data{1,1},1); % Number of channels recorded
        
    % Calculate onset based on EI 
    [EI_tbl,Nd_tbl,~,~] = ms_epi_ind(pat_data);
    tbl_plhg = table();

    % Create matrices to store onset channels for each method
    onset_imprint = zeros(n_chan, sz_count);
    onset_ei = onset_imprint;
    onset_plhg = onset_imprint;
    onset_time_imprint = zeros(1,sz_count);

    when_onset_ei = zeros(1, sz_count);
       
    % Iterate over all patient seizures
    for sz = 1:sz_count
        % Onset based on seizure imprint
        sz_imprint = pat_imprints.cell_imprint{sz,1};
        imprint_wind = zeros(size(sz_imprint,1),size(sz_imprint,2)-(wdw_sz-1));
        for epoch = 1:(size(sz_imprint,2)-(wdw_sz-1))
            imprint_wind(:,epoch) = sum(sz_imprint(:,epoch:epoch+(wdw_sz-1)),2);
        end
        imprint_wind = imprint_wind==wdw_sz;
        chan_sum = sum(imprint_wind);
        % Identify when first activity is found
        onset_time = find(chan_sum,1,'first');
       
        if isempty(onset_time)
            fprintf("Patient %s, seizure %s: No imprint onset detected \n", pat_data.patient_id{1}, string(pat_imprints.segment_id(sz)))
            onset_time_imprint(sz) = NaN;
        else
            onset_imprint(:,sz) = sum(imprint_wind(:,(onset_time:onset_time+det)),2) > 0;
            onset_time_imprint(sz) = onset_time; %Imprint
        end
        
        % Onset based on Epileptigenicity Index
        %[EI,Nd,time,ERmaster,Na] = epileptogenicityIndex(ict');
        onset_ei(EI_tbl{sz}>0.3,sz) = 1;
    
        % Onset based on Phase-locked high-gamma
        [tbl_plhg(sz,:)] = ms_PLHG(pat_data(sz,:));
%             [~,plhg_id] = mink(tbl_plhg(sz,:).when_plhg,4);
%             onset_plhg(plhg_id,sz) = 1;
        when_plhg = tbl_plhg(sz,:).when_plhg;
        when_plhg_na_omit = when_plhg(~isnan(when_plhg));
        % Onset within 0.5 of earliest PLHG included in seizure onset
        earliest = min(when_plhg_na_omit);
        latest = earliest + 0.5;
        onset_plhg(:,sz) = when_plhg >= earliest & when_plhg<=latest;

        % determine time of onset using EI
        when_onset_ei(sz) = Nd_tbl{sz}(EI_tbl{sz}==1); 

    end
    when_onset = [onset_time_imprint' when_onset_ei'  min(tbl_plhg.when_plhg,[],2,'omitnan')];

    % Populate table 
    onset_output.imprint_chan = mat2cell(onset_imprint,n_chan,sz_count);
    onset_output.EI_chan = mat2cell(onset_ei,n_chan,sz_count);
    onset_output.PLHG_chan = mat2cell(onset_plhg,n_chan,sz_count);
    onset_output.when_onset =  mat2cell(when_onset,sz_count,3);
end