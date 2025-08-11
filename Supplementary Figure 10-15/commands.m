%This code makes the figures for the noise distributions of parameters in
%single fits, the single fit values and the combinations values
%for EC50 and Hill slope parameters

clear variables

load("noiseParameters.mat")

combinationFits.D7.PROT = readtable("combination_fits_ipscs.xlsx",'Sheet','D7 PROT','ReadRowNames',true,'VariableNamingRule','modify');
combinationFits.D7.CD4ISP = readtable("combination_fits_ipscs.xlsx",'Sheet','D7 CD4ISP','ReadRowNames',true,'VariableNamingRule','modify');
combinationFits.D7.CD3N = readtable("combination_fits_ipscs.xlsx",'Sheet','D7 CD3N','ReadRowNames',true,'VariableNamingRule','modify');
combinationFits.D14.CD3N = readtable("combination_fits_ipscs.xlsx",'Sheet','D14 CD3N','ReadRowNames',true,'VariableNamingRule','modify');
combinationFits.D14.CD3P = readtable("combination_fits_ipscs.xlsx",'Sheet','D14 CD3P','ReadRowNames',true,'VariableNamingRule','modify');
combinationFits.D14.SP = readtable("combination_fits_ipscs.xlsx",'Sheet','D14 SP','ReadRowNames',true,'VariableNamingRule','modify');

singleFits.D7.PROT = readtable("single_fits_ipscs.xlsx","Sheet","D7 PROT",'ReadRowNames',true);
singleFits.D7.CD4ISP = readtable("single_fits_ipscs.xlsx","Sheet","D7 CD4ISP",'ReadRowNames',true);
singleFits.D7.CD3N = readtable("single_fits_ipscs.xlsx","Sheet","D7 CD3N",'ReadRowNames',true);
singleFits.D14.CD3N = readtable("single_fits_ipscs.xlsx","Sheet","D14 CD3N",'ReadRowNames',true);
singleFits.D14.CD3P = readtable("single_fits_ipscs.xlsx","Sheet","D14 CD3P",'ReadRowNames',true);
singleFits.D14.SP = readtable("single_fits_ipscs.xlsx","Sheet","D14 SP",'ReadRowNames',true);


%% Days 7-14
p.days = 'D7';

    %Pro-T
    p.cell = "PROT";

        %SCF-IL7
        p.int = "SCF_IL7";
        p.drug1 = "SCF";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL3-IL7
        p.int = "IL3_IL7";
        p.drug1 = "IL3";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL7-TNFa
        p.int = "IL7_TNFa";
        p.drug1 = "IL7";
        p.drug2 = "TNFa";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %TNFa_CXCL12
        p.int = "TNFa_CXCL12";
        p.drug1 = "TNFa";
        p.drug2 = "CXCL12";
        p.drug1_mono = 0;
        p.drug2_mono = 1;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);


   %CD4ISP
    p.cell = "CD4ISP";

        %SCF-IL7
        p.int = "SCF_IL7";
        p.drug1 = "SCF";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL3-IL7
        p.int = "IL3_IL7";
        p.drug1 = "IL3";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL7-TNFa
        p.int = "IL7_TNFa";
        p.drug1 = "IL7";
        p.drug2 = "TNFa";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %TNFa_CXCL12
        p.int = "TNFa_CXCL12";
        p.drug1 = "TNFa";
        p.drug2 = "CXCL12";
        p.drug1_mono = 0;
        p.drug2_mono = 1;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

    %CD3N
    p.cell = "CD3N";

        %SCF-IL7
        p.int = "SCF_IL7";
        p.drug1 = "SCF";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL3-IL7
        p.int = "IL3_IL7";
        p.drug1 = "IL3";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL7-TNFa
        p.int = "IL7_TNFa";
        p.drug1 = "IL7";
        p.drug2 = "TNFa";
        p.drug1_mono = 0;
        p.drug2_mono = 1;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %TNFa_CXCL12
        p.int = "TNFa_CXCL12";
        p.drug1 = "TNFa";
        p.drug2 = "CXCL12";
        p.drug1_mono = 1;
        p.drug2_mono = 1;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

%% Days 14-28
p.days = 'D14';

    %CD3N
    p.cell = "CD3N";

        %SCF-IL7
        p.int = "SCF_IL7";
        p.drug1 = "SCF";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL3-IL7
        p.int = "IL3_IL7";
        p.drug1 = "IL3";
        p.drug2 = "IL7";
        p.drug1_mono = 1;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL7-TNFa
        p.int = "IL7_TNFa";
        p.drug1 = "IL7";
        p.drug2 = "TNFa";
        p.drug1_mono = 0;
        p.drug2_mono = 1;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);


   %CD3P
    p.cell = "CD3P";

        %SCF-IL7
        p.int = "SCF_IL7";
        p.drug1 = "SCF";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL3-IL7
        p.int = "IL3_IL7";
        p.drug1 = "IL3";
        p.drug2 = "IL7";
        p.drug1_mono = 1;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL7-TNFa
        p.int = "IL7_TNFa";
        p.drug1 = "IL7";
        p.drug2 = "TNFa";
        p.drug1_mono = 0;
        p.drug2_mono = 1;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);


    %8SP
    p.cell = "SP";

        %SCF-IL7
        p.int = "SCF_IL7";
        p.drug1 = "SCF";
        p.drug2 = "IL7";
        p.drug1_mono = 0;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL3-IL7
        p.int = "IL3_IL7";
        p.drug1 = "IL3";
        p.drug2 = "IL7";
        p.drug1_mono = 1;
        p.drug2_mono = 0;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);

        %IL7-TNFa
        p.int = "IL7_TNFa";
        p.drug1 = "IL7";
        p.drug2 = "TNFa";
        p.drug1_mono = 0;
        p.drug2_mono = 1;
    
        fun_Dist(p,noiseParam,singleFits,combinationFits);
