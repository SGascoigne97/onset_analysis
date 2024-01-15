% Create table summarising metadata
final_output.epilepsy_duration = final_output.age - final_output.onset_age;

% final_output_fav = final_output(final_output.outcome <3,:);
% final_output_unfav = final_output(final_output.outcome >2,:);

cont_vars = ["age", "epilepsy_duration"];
cat_vars = ["outcome", "sex", "Op type"];

tab_rows = ["outcome count", cont_vars]';

tab_fill = nan(length(cont_vars),1);

meta_tab_cont = table(cont_vars', repelem("", length(cont_vars),1),tab_fill,...
    tab_fill, tab_fill, tab_fill, tab_fill, tab_fill, tab_fill,...
    tab_fill, tab_fill,'VariableNames',...
    ["metadata", "var_type", "mean/count", "sd", sprintf("mean/count %s",out_grps(1)),...
    sprintf("sd %s",out_grps(1)), sprintf("mean/count %s",out_grps(2)),...
    sprintf("sd %s",out_grps(2)),"t","df","p"]);

for var = cont_vars
    data = final_output.(sprintf(var));
    data_fav = final_output(final_output.outcome_cat == out_grps(1),:).(sprintf(var));
    data_unfav = final_output(final_output.outcome_cat == out_grps(2),:).(sprintf(var));
    meta_tab_cont(meta_tab_cont.metadata == var,2:8) = table("continuous",...
        mean(data), std(data), mean(data_fav), std(data_fav),...
        mean(data_unfav), std(data_unfav));
    [~,p,~,st] = ttest2(data_fav,data_unfav);
    meta_tab_cont(meta_tab_cont.metadata == var,:).t = st.tstat;
    meta_tab_cont(meta_tab_cont.metadata == var,:).df = st.df;
    meta_tab_cont(meta_tab_cont.metadata == var,:).p = p;
end
%%
for var = cat_vars
    cate = categorical(final_output.(sprintf(var)));
    cat_f = categorical(final_output(final_output.outcome <3,:).(sprintf(var)));
    cat_u = categorical(final_output(final_output.outcome >2,:).(sprintf(var)));
    cats = categories(cate);
    cat_var_tab = table(cats, countcats(cate),...
        nan(length(cats),1), nan(length(cats),1),...
        'VariableNames', ["Group", "All",out_grps]);

    for cat_ind = 1:length(cats)
        cat_a_ind = cate(cate == cats{cat_ind});
        cat_f_ind = cat_f(cat_f == cats{cat_ind});
        cat_u_ind = cat_u(cat_u == cats{cat_ind});
        cat_var_tab(cat_ind,2:4) = table(length(cat_a_ind),...
            length(cat_f_ind),length(cat_u_ind));
    end
    cat_tabs.(strrep(sprintf(var)," ", "_")) = cat_var_tab;
end



