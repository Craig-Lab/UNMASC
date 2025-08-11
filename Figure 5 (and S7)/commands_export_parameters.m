clear variables

load("Results/fit_parameters.mat")

interactions = ["SCF_IL7", "IL3_IL7","IL7_TNFa","TNFa_CXCL12"];


%% Set monotonicity of interactions

%Days 7-14
mono.D7.PROT.SCF_IL7 = 0;
mono.D7.PROT.IL3_IL7 = 0;
mono.D7.PROT.IL7_TNFa = 0;
mono.D7.PROT.TNFa_CXCL12 = 2;

mono.D7.CD4ISP.SCF_IL7 = 0;
mono.D7.CD4ISP.IL3_IL7 = 0;
mono.D7.CD4ISP.IL7_TNFa = 0;
mono.D7.CD4ISP.TNFa_CXCL12 = 2;

mono.D7.CD3N.SCF_IL7 = 0;
mono.D7.CD3N.IL3_IL7 = 0;
mono.D7.CD3N.IL7_TNFa = 0;
mono.D7.CD3N.TNFa_CXCL12 = 3;

%Days 14-28
mono.D14.CD3N.SCF_IL7 = 0;
mono.D14.CD3N.IL3_IL7 = 1;
mono.D14.CD3N.IL7_TNFa = 2;

mono.D14.CD3P.SCF_IL7 = 0;
mono.D14.CD3P.IL3_IL7 = 1;
mono.D14.CD3P.IL7_TNFa = 2;

mono.D14.SP.SCF_IL7 = 0;
mono.D14.SP.IL3_IL7 = 1;
mono.D14.SP.IL7_TNFa = 2;


%% Create table

%DAYS 7-14
p.days = 'D7';

sz = [14 4]; % 4 combinations for Days 7-14, 14 parameters max

varTypes = ["double", "double", "double", "double"];
varNames = ["SCF-IL7", "IL3-IL7", "IL7-TNFa","TNFa-CXCL12"];

rowNames = ["E0", "E1", "E2", "E3", "h11", "h12", "h21", "h22", "C11", "C12", "C21", "C22", "alpha12", "alpha21"];

    % PROT CELLS
    p.cell = 'PROT';
        resTable.(p.days).(p.cell) = fun_table(p,interactions, mono, sz,varTypes,varNames,rowNames,fit_parameters);

    %CD4ISP CELLS
    p.cell = "CD4ISP";
        resTable.(p.days).(p.cell) = fun_table(p,interactions, mono, sz,varTypes,varNames,rowNames,fit_parameters);

   %DP (CD3-) CELLS
   p.cell = "CD3N";
       resTable.(p.days).(p.cell) = fun_table(p,interactions, mono, sz,varTypes,varNames,rowNames,fit_parameters);

%DAYS 14-28
p.days = 'D14';

sz = [14 3]; % 4 combinations for Days 7-14, 14 parameters max

varTypes = ["double", "double", "double"];
varNames = ["SCF-IL7", "IL3-IL7", "IL7-TNFa"];

    %DP (CD3-) CELLS
    p.cell = 'CD3N';
        resTable.(p.days).(p.cell) = fun_table(p,interactions, mono, sz,varTypes,varNames,rowNames,fit_parameters);

    %DP (CD3+) CELLS
    p.cell = "CD3P";
        resTable.(p.days).(p.cell) = fun_table(p,interactions, mono, sz,varTypes,varNames,rowNames,fit_parameters);

   %8SP CELLS
   p.cell = "SP";
       resTable.(p.days).(p.cell) = fun_table(p,interactions, mono, sz,varTypes,varNames,rowNames,fit_parameters);

%% Export

filename = "Results/combination_fits_ipscs.xlsx";

writetable(resTable.D7.PROT,filename,'Sheet','D7 PROT','WriteRowNames',true);
writetable(resTable.D7.CD4ISP,filename,'Sheet','D7 CD4ISP','WriteRowNames',true);
writetable(resTable.D7.CD3N,filename,'Sheet','D7 CD3N','WriteRowNames',true);
writetable(resTable.D14.CD3N,filename,'Sheet','D14 CD3N','WriteRowNames',true);
writetable(resTable.D14.CD3P,filename,'Sheet','D14 CD3P','WriteRowNames',true);
writetable(resTable.D14.SP,filename,'Sheet','D14 SP','WriteRowNames',true);