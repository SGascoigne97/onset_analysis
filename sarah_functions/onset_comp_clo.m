function [all_pat_table] = onset_comp_clo(final_output, opts)
% Compute seizure imprint based on EEG recordings

% input:
%   - 

% output
%   - 

    arguments
        final_output
        opts.atl (1,1) double {mustBeMember(opts.atl, [72, 120, 250, 500])}= 120; % Lausanne atlas to convert channels to ROIs
        opts.per_or_all (1,1) {mustBeMember(opts.per_or_all, ["per_sz", "all_sz"])} = "per_sz" % Compute value for each seizure or for onset across seizures
        opts.det (1,1) {mustBeMember(opts.det, ["imprint", "EI", "PLHG"])} = "imprint" % Automatic seizure onset localisation algorithm
        opts.sz_prop_thresh (1,1) double = 0.25 % If computing imprint summarised across seizures, proportion of seizures in which region is required to be included in summary onset 
        opts.sz_types = "focal" % State the seizure types to include
    end
    
    %fill in optional arguments
    atl = opts.atl;
    per_or_all = opts.per_or_all;
    det = opts.det;
    sz_prop_thresh = opts.sz_prop_thresh;
    sz_types = opts.sz_types;

    % Print parameters to the console
    fprintf("Lausanne %d results for %s using %s onset detection \n", atl, per_or_all, det)

    % Itearate across patients
    for pat = 1:size(final_output,1)
        patient = final_output.Patient_id{pat};
        pat_onset = final_output(string(final_output.Patient_id) == patient,:);
        if all(isnan(pat_onset.clo_chan{:})) | sum(pat_onset.clo_chan{:}) ==0
            fprintf("%s does not have clinically labelled onset \n", patient)
             continue
        end
        
        pat_clo = pat_onset.(sprintf("clo_roi_%d", atl)){:};
        pat_auto = pat_onset.(sprintf("%s_roi_%d", det, atl)){:};
        
        % Remove any seizures with onset in >50% of regions or no onset
        % detected
        rm_sz =  sum(pat_auto,1) == 0 |...
            sum(pat_auto,1) >= size(pat_auto,1)/2;
        pat_auto = pat_auto(:,~rm_sz);
        if size(pat_auto,2) == 0
            fprintf("%s all seizures removed \n", patient)
            continue
        end
        pat_sz_ids = string(pat_onset.Segment_ids{:});
        pat_sz_id = pat_sz_ids(~rm_sz);
        pat_sz_types = pat_onset.sz_types{:}(~rm_sz);
        
        % Isolate the seizure types you're interesterested in
        if sz_types ~= "all"
            keep_sz = ismember(pat_sz_types, sz_types);
            pat_auto = pat_auto(:,keep_sz);
            pat_sz_id = pat_onset.Segment_ids{:}(keep_sz);
            pat_sz_types = pat_onset.sz_types{:}(keep_sz);
        end
        
        % Compare each imprint onset against CLO (with permutation test)
        pat_comp_tab = comp_auto_clo(pat_clo,...
            pat_auto, "per_sz_or_all_sz", per_or_all,...
            "sz_prop_thresh", sz_prop_thresh);

        if isempty(pat_comp_tab)
            continue
        end
        
        pat_surg_out = pat_onset.Surgery_outcome{:};
        pat_surg_year = pat_onset.("Surgery year");
        pat_out_year = pat_onset.("Outcome year"){:};
        year_1 = pat_out_year - pat_surg_year == 1;

        pat_comp_tab.patient_id = repmat(string(patient),size(pat_comp_tab,1),1);
        pat_comp_tab.outcome = repmat(pat_surg_out(year_1),size(pat_comp_tab,1),1);
        pat_comp_tab.op_type = repmat(pat_onset.("Op type"),size(pat_comp_tab,1),1);
        
        % Identify whether results are per seizure or summarised across seizures
        if per_or_all == "per_sz"
            pat_comp_tab.sz_id = pat_sz_id;
            pat_comp_tab.sz_type = pat_sz_types;
        else
            pat_comp_tab.sz_id = {pat_sz_id};
            pat_comp_tab.sz_type = {pat_sz_types};
        end

        % Add patient to table
        if exist('all_pat_table', 'var')
            all_pat_table = [all_pat_table; pat_comp_tab];
        else
            all_pat_table = pat_comp_tab;
        end
    end
end
