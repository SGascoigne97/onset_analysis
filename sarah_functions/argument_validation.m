function [chan_or_roi, n_perm, out_thresh, consensus_thresh, most_resec_thresh,...
    det_meths, save_fig, save_folder, file_type, out_grps] = ...
    argument_validation(data_tbl, chan_or_roi, n_perm, out_thresh, consensus_thresh, ...
    most_resec_thresh, det_meths, save_fig, save_folder, file_type)

    % Validate that all parameters used throughout our analyses are in the
    % correct format for use in this code

    % Sarah Jane Gascoigne
    % 29/08/2023
    
    % input:
    %   - data_tbl
    %   - chan_or_roi 
    %   - n_perm
    %   - out_thresh
    %   - consensus_thresh
    %   - most_resec_thresh
    %   - det_meths
    %   - save_fig
    %   - save_folder
    %   - file_type
    
    % output
    %   - chan_or_roi 
    %   - n_perm
    %   - out_thresh
    %   - consensus_thresh
    %   - most_resec_thresh
    %   - det_meths
    %   - save_fig
    %   - save_folder
    %   - out_grps
    
    arguments
        data_tbl table
        chan_or_roi 
        n_perm 
        out_thresh
        consensus_thresh
        most_resec_thresh
        det_meths
        save_fig
        save_folder
        file_type (1,1) string {ismember(file_type, ["png", "svg"])}
    end

    %% Validating chan_or_roi
    if isnumeric(chan_or_roi)
        if any(chan_or_roi == [1,2,3,4])
            lausanne_atlases = ["roi_72", "roi_120", "roi_250", "roi_500"];
            chan_or_roi = lausanne_atlases(chan_or_roi);
        elseif any(chan_or_roi == [72,120,250,500])
            chan_or_roi = sprintf("roi_%d",chan_or_roi);
        end
    elseif isstring(chan_or_roi)
        if any(chan_or_roi == ["72","120","250","500"])
            chan_or_roi = sprintf("roi_%s",chan_or_roi);
        elseif contains(chan_or_roi, "chan")
            chan_or_roi = "chan";
        end
    end

    if ~ismember(string(chan_or_roi), ["chan", "roi_72", "roi_120", "roi_250", "roi_500"])
        fprintf("chan_or_roi = %s is not a valid argument, it should be chan, roi_72, roi_120, roi_250, or roi_500 \n", string(chan_or_roi))
        chan_or_roi = [];
    end

    %% Validating n_perm
    if isnumeric(n_perm)
        if n_perm ~= round(n_perm) % This checks if the value is an integer as non integers will have a ~= round(a)
            fprintf("n_perm = %s is not a valid argument, it must be a positive integer greater than 100 \n", string(n_perm))
            n_perm = [];
        elseif n_perm<100
            fprintf("n_perm = %s is not a valid argument, it must be a positive integer greater than 100 \n", string(n_perm))
            n_perm = [];
        end
    else
        fprintf("n_perm = %s is not a valid argument, it must be a positive integer greater than 100 \n", string(n_perm))
        n_perm = [];
    end

    %% Validating out_thresh
    if isnumeric(out_thresh)
        if ~ismember(out_thresh, 1:6)
            fprintf("out_thresh = %d is not a valid argument, it must be between one and six \n", out_thresh)
            out_thresh = [];
        end
    else
        fprintf("out_thresh = %d is not a valid argument, it must be a NUMERIC between one and six \n", out_thresh)
        out_thresh = [];
    end

    % Create a variable out_grps for tables and figures
    if out_thresh == 1
        out_grps = ["ILAE 1", "ILAE 2+"];
    else 
        out_grps = [sprintf("ILAE 1-%d", out_thresh), sprintf("ILAE %d+", out_thresh+1)];
    end

    %% Validating consensus_thresh
    if isnumeric(consensus_thresh)
        if consensus_thresh < 0 
            fprintf("consensus_thresh = %.2f is not a valid argument, it must be greater than zero and less than or equal to one \n", consensus_thresh)
            consensus_thresh = [];
        elseif consensus_thresh > 1
            if consensus_thresh <= 100
                consensus_thresh = consensus_thresh/100;
            else
                fprintf("consensus_thresh = %.2f is not a valid argument, it must be greater than zero and less than or equal to one \n", consensus_thresh)
                consensus_thresh = [];
            end
        end
    else
        fprintf("consensus_thresh = %.2f is not a valid argument, it must be greater than zero and less than or equal to one \n", consensus_thresh)
        consensus_thresh = [];
    end

     %% Validating most_resec_thresh
    if isnumeric(most_resec_thresh)
        if most_resec_thresh < 0 
            fprintf("most_resec_thresh = %.2f is not a valid argument, it must be greater than zero and less than or equal to one \n", most_resec_thresh)
            most_resec_thresh = [];
        elseif most_resec_thresh > 1
            if most_resec_thresh <= 100
                most_resec_thresh = most_resec_thresh/100;
            else
                fprintf("most_resec_thresh = %.2f is not a valid argument, it must be greater than zero and less than or equal to one \n", most_resec_thresh)
                most_resec_thresh = [];
            end
        end
    else
        fprintf("most_resec_thresh = %.2f is not a valid argument, it must be greater than zero and less than or equal to one \n", most_resec_thresh)
        most_resec_thresh = [];
    end

    %% Validating det_meths
    % Here we will validate that the required onsets are available in the
    % data table
    tbl_vars = data_tbl.Properties.VariableNames;
    valid_meth = nan(1,length(det_meths));
    for meth_ind = 1:length(det_meths)
        valid_meth(meth_ind) = any(contains(tbl_vars, sprintf("%s_%s", det_meths(meth_ind), chan_or_roi)));
        if valid_meth(meth_ind) == 0
            fprintf("%s is not a valid onset detection method as %s_%s is not included in data table \n", det_meths(meth_ind),det_meths(meth_ind), chan_or_roi )
        end
    end

    if any(valid_meth ~= 1)
        det_meths = [];
    end

    %% Validating save_fig
    if ~ismember(save_fig, 0:1)
        fprintf("save_fig = %d is not a valid argument, it must be 1 (save figures) or 0 (do not save figures) \n", save_fig)
        save_fig = [];
    elseif (string(save_fig) == "true" | string(save_fig) == "yes")
        save_fig = 1;
    elseif string(save_fig) == "false" | string(save_fig) == "no"
        save_fig = 0;
    end

    %% Validating save_folder
    if isstring(save_folder)
        save_folder = char(save_folder);
    end

    if ischar(save_folder)
        if ~(save_folder(end)=='/')
            save_folder(end+1) = '/';
        end
        if ~exist(save_folder, "dir")
            mkdir(save_folder)
        end
    else
        fprintf("save_folder is not in the correct format, please provide string or char \n")
        save_folder = [];
    end

    if any(cellfun(@isempty,{chan_or_roi, n_perm, out_thresh, consensus_thresh, ...
            most_resec_thresh, det_meths, save_fig, save_folder}))
        fprintf("Plese review arguments and ammend where necessary \n")
        chan_or_roi = [];
        n_perm = [];
        out_thresh = [];
        consensus_thresh = [];
        most_resec_thresh = [];
        det_meths = [];
        save_fig = [];
        save_folder = [];
    else
        % Print chan_or_roi
        fprintf("\nAll arguments are valid, please verify they match your requirements \n")
        if contains(string(chan_or_roi), "roi")
            fprintf("Analyses will be using %s atlas \n", chan_or_roi)
        elseif string(chan_or_roi) == "chan"
            fprintf("Analyses will be channel-wise \n")
        end
        % Print n_perm
        fprintf("P-values for AUCs will be computing using %d permutations \n", n_perm)
        % Print out_thresh
        fprintf("%s will be considered favourable outcomes, %s will be considered unfavourable \n", out_grps(1), out_grps(2))
        % Print consensus_thresh
        fprintf("Consensus onsets will include regions/channels included in at least %.d%% of seizures \n", floor(consensus_thresh*100))
        % Print most_resec_thresh
        fprintf("Subjects with at least %.d%% of onsets resected will be considered to have 'most' resected \n", floor(most_resec_thresh*100))
        % Print det_meths
        fprintf("Onsets will be detected using the following detection methods: \n ")
        for meth = det_meths
            fprintf("\t - %s \n", meth)
        end
        % Print save_fig
        if save_fig == 1
            fprintf("All plots will be saved as %ss in the following location: \n %s \n", file_type, string(save_folder))
        else
            fprintf("Plots will not be saved \n")
        end  
    end
end