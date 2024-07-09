function [] = half_violin_upd(data_tab, comp_var, grps, opts)
% Create half-violin plot to compare distributions between two groups

% input:
%   - data_tbl: full data table
%   - optional inputs
%       - sz_type: seizure type of interest
%       - min_sz_duration: minimum duration for each seizure 
%       - min_sz_count: minimum number of seizures to ahve been recorded per patient

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tab table
        comp_var string
        grps (1,:) double
        opts.bin_wds (:,1) double = NaN % One value per group 
        opts.y_lim (:,2) double = [NaN,NaN]
        opts.grp_names = ["Group One", "Group Two"]
        opts.legend_pos (1,1) string = "northeast"
        opts.save_fig = 0;
        opts.save_loc = "figures/"
        opts.file_type = "png"
        opts.new_fig (1,1) double = 1

    end
    
    %fill in optional arguments
    bin_wds = opts.bin_wds;
    y_lim = opts.y_lim;
    save_fig = opts.save_fig;
    save_loc = opts.save_loc;
    file_type = opts.file_type;
    grp_names = opts.grp_names;
    new_fig = opts.new_fig;
    legend_pos = opts.legend_pos;

    grp_vals = unique(grps);

    offset = (size(data_tab,1)+10)/2;
    max_x = offset + 10;
    x_lim = [-10, max_x];
    
    if new_fig == 1
        f = figure();
        f.Position = [500,500,400*length(comp_var),450];
        tiledlayout(1,length(comp_var))
    else 
        if length(comp_var)>1
            fprintf("Must be new figure if including multiple variables")
            return
        end
    end
    
    var_ind = 1;
    for var = comp_var
        data_tab_var = data_tab(~isnan(data_tab.(sprintf(var))),:);
        grps_var = grps(~isnan(data_tab.(sprintf(var))));
        if new_fig == 1
            nexttile
        end
        data = data_tab_var.(sprintf(var));

        hold on
        if length(bin_wds) == 1
            if isnan(bin_wds)
                bin_wd = 5*range(data)/length(data);
            else
                bin_wd = bin_wds;
            end
        else
            if isnan(bin_wds(var_ind))
                bin_wd = 5*range(data)/length(data);
            else
                bin_wd = bin_wds(var_ind);
            end
        end
        grp_a = data(grps_var == grp_vals(1));
        grp_b = data(grps_var == grp_vals(2));
                  
        min_val_a = round(min(grp_a)*(1/bin_wd))*bin_wd;
        min_val_b = round(min(grp_b)*(1/bin_wd))*bin_wd;
        max_val_a = round(max(grp_a)*(1/bin_wd))*bin_wd;
        max_val_b = round(max(grp_b)*(1/bin_wd))*bin_wd;
        
        if min_val_a > min(grp_a)
            min_val_a = min_val_a - bin_wd;
        end 
        if min_val_b > min(grp_b)
            min_val_b = min_val_b - bin_wd;
        end
        if max_val_a < max(grp_a)
            max_val_a = max_val_a + bin_wd;
        end
        if max_val_b < max(grp_b)
            max_val_b = max_val_b + bin_wd;
        end
        
        if min_val_a == max_val_a
            max_val_a = max_val_a+bin_wd;
        end
        if min_val_b == max_val_b
            max_val_b = max_val_b+bin_wd;
        end
        
        hist_vals_a = histcounts(grp_a, min_val_a:bin_wd:max_val_a);
        % Add scatter points
        non_zero_grps = hist_vals_a(hist_vals_a~=0);
        jitt = [];
        for grp_sz = non_zero_grps
            jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
        end
        scatter(-jitt, sort(grp_a), [], [0.3010 0.7450 0.9330], 'filled');
        
        hist_vals_b = histcounts(grp_b, min_val_b:bin_wd:max_val_b);
        % Add scatter points
        non_zero_grps = hist_vals_b(hist_vals_b~=0);
        jitt = [];
        for grp_sz = non_zero_grps
            jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
        end
        scatter(jitt, sort(grp_b), [],[0.8500 0.3250 0.0980], 'filled');
        
        if length(hist_vals_a) >1
            smooth_data_g = 1.05*(smoothdata(interp1(hist_vals_a, 1:0.01:length(hist_vals_a))));
            fill(-[0 smooth_data_g 0], ...
                rescale((1:length([0 smooth_data_g 0]))/length([0 smooth_data_g 0]),...
                min(grp_a) - (bin_wd/10), max(grp_a) + (bin_wd/10)),[0.3010 0.7450 0.9330], 'FaceAlpha',0.3,...
                'LineStyle','none')
        end 
        plot([-2,0], median(grp_a)*[1,1], "LineWidth", 2.5, "Color",  0.75*[0.3010 0.7450 0.9330] )
        
        if length(hist_vals_b) > 1
            smooth_data_b = 1.05*(smoothdata(interp1(hist_vals_b, 1:0.01:length(hist_vals_b))));
            fill([0 smooth_data_b 0],...
                rescale((1:length([0 smooth_data_b 0]))/length([0 smooth_data_b 0]),...
                min(grp_b) - (bin_wd/10), max(grp_b) + (bin_wd/10)), [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,...
                'LineStyle','none')
        end
        plot([2,0], median(grp_b)*[1,1], "LineWidth", 2.5, "Color", 0.75*[0.8500 0.3250 0.0980] )
        set(gca, "XTick", [])
        if any(~isnan(y_lim(:,1)))
            if size(y_lim,1) == 1
                ylim(y_lim);
            else
                if all(~isnan(y_lim(var_ind,:)))
                    ylim(y_lim(var_ind,:));
                end
            end
            
        end
        box off
        title(sprintf("%s (n=%d)", strrep(var, "_", " "), length(data)))

        var_ind = var_ind +1;
    end

    legend(grp_names, "Location", legend_pos, "box", "off")

    if save_fig == 1
        if new_fig == 0
            fprintf("Save figure outside of function if not creating a new figure \n")
        else
            saveas(f, sprintf("%s.%s", save_loc, file_type))
        end
    end
end