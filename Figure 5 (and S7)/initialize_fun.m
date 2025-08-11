function p = initialize_fun(p,cytokines,cytConc,coeff)
    %Find index of agent name
    p.idx1 =  find(contains(cytokines,p.drug1_name));
    p.idx2 =  find(contains(cytokines,p.drug2_name));
    
    %Find the c0 of the scaled concentration
    p.c01 = cytConc(3,p.idx1);
    p.c02 = cytConc(3,p.idx2);
    
    %Find the interaction name
    p.int_name = strcat(p.drug1_name,"_",p.drug2_name);
    
    %Find the concentration of each agent (ng/mL)
    p.drug1_conc = cytConc(1:5,p.idx1);
    p.drug2_conc = cytConc(1:5,p.idx2);
    
    %Get the polynomial coefficients for each agent and their interaction
    p.drug1_coeff = coeff.(p.days).(p.cell).(p.drug1_name);
    p.drug2_coeff = coeff.(p.days).(p.cell).(p.drug2_name);
    p.int_coeff = coeff.(p.days).(p.cell).(p.int_name);
end 