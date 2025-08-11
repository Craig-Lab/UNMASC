function tb = table_fun(p,cytokines,sz,varTypes,varNames,rowNames,results)

    tb = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames,'RowNames',rowNames);
    
    for i = 1:sz(2)
      p.name = num2str(cell2mat(cytokines(i)));
      new_vec = reshape(results.(p.days).(p.cell).(p.name),[],1);
      if length(new_vec) == 12
          new_vec = [new_vec;zeros(12,1)];
      end 
      tb(:,i) = array2table(new_vec);
    end 

end 