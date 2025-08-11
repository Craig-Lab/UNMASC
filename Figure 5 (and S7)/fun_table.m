function tb = fun_table(p,interactions,mono,sz,varTypes,varNames,rowNames,results)

    tb = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames,'RowNames',rowNames);

    for i = 1:sz(2)
        p.name = interactions(i);
        p.monotonic = mono.(p.days).(p.cell).(p.name);

        new_vec = reshape(results.(p.days).(p.cell).(p.name),[],1);

        if p.monotonic == 0 %none of the drugs are monotonic
            tb(:,i) = array2table(new_vec);
        elseif p.monotonic == 1 %first drug is monotonic
            tb(1:5,i) = array2table(new_vec(1:5));
            tb(6,i) = {NaN};
            tb(7:9,i) = array2table(new_vec(6:8));
            tb(10,i) = {NaN};
            tb(11:14,i) = array2table(new_vec(9:12));
        elseif p.monotonic ==2 %second drug is monotonic 
            tb(1:7,i) = array2table(new_vec(1:7));
            tb(8,i) = {NaN};
            tb(9:11,i) = array2table(new_vec(8:10));
            tb(12,i) = {NaN};
            tb(13:14,i) = array2table(new_vec(11:12));
        elseif p.monotonic == 3 %both drugs are monotonic
            tb(1:5,i) = array2table(new_vec(1:5));
            tb(6,i) = {NaN};
            tb(7,i) = array2table(new_vec(6));
            tb(8,i) = {NaN};
            tb(9,i) = array2table(new_vec(7));
            tb(10,i) = {NaN};
            tb(11,i) = array2table(new_vec(8));
            tb(12,i) = {NaN};
            tb(13:14,i) = array2table(new_vec(9:10));
        end 
    end 

end 