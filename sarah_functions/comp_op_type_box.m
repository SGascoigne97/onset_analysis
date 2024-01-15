function [] = comp_op_type_box(data_tbl, comp_measures, opts)
    % Compute AUC and associated p-value for logistic regression models
    % predicting surgical outcome
    
    % Sarah Jane Gascoigne
    % 24/08/2023
    
    % input:
    %   - data_tbl: full data table
    
    % output
    %   - mod_details: struct with relevant information for logistic regression
    %   models
    
        arguments
            data_tbl table
            comp_measures (1,:) string
            opts.save_fig (1,1) double {mustBeNumeric} = 0
            opts.save_loc (1,1) string = "figures/"
            opts.comp_name (1,1) string = ""
        end
        
        %fill in optional arguments
        save_fig = opts.save_fig;
        save_loc = opts.save_loc;
        comp_name = opts.comp_name;
    
    figure()
    tiledlayout(2,length(comp_measures))
    for op = ["T Lx", "F Lx"]
        for comp = comp_measures
            nexttile
            hold on
            boxchart(categorical(data_tbl(data_tbl.("Op type") == op,:).outcome_cat),...
                data_tbl(data_tbl.("Op type") == op,:).(sprintf(comp)), "MarkerStyle","none")
            swarmchart(categorical(data_tbl(data_tbl.("Op type") == op,:).outcome_cat),...
                data_tbl(data_tbl.("Op type") == op,:).(sprintf(comp)), "filled")
            hold off
            title(sprintf("%s for %s", comp, op))
        end
    end

    if save_fig == 1
        saveas(fig, sprintf("%s%s_boxplots.svg", save_loc, comp_name))
    end
end