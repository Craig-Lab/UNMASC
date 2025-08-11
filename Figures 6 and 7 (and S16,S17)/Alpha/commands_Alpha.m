%% Load data 

clear variables

%Load cytokine concentrations and set names
cytConc = readtable('../../Data/cytConc.xlsx');
cytokines = cytConc.Properties.VariableNames;
cytConc = table2array(cytConc);

%Load surfaces results
load("fit_parameters.mat") %Fit parameters
load("parameters_MC.mat") %Parameters following Monte Carlo simulations
% load("CI.mat") %Confidence intervals for the parameters


%% Alpha values D7-14

%Set alpha values in 'struct' object

%DAYS 7-14 
p.days = "D7";
p.cell = "PROT";

    p.int = "SCF_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL3_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL7_TNFa";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "TNFa_CXCL12";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));


p.cell = "CD4ISP";

    p.int = "SCF_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL3_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL7_TNFa";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "TNFa_CXCL12";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));


p.cell = "CD3N";

    p.int = "SCF_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL3_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL7_TNFa";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "TNFa_CXCL12";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));
    

    data_alpha12_D7 = [alpha12.D7.PROT.SCF_IL7,alpha12.D7.CD4ISP.SCF_IL7,alpha12.D7.CD3N.SCF_IL7,...
            alpha12.D7.PROT.IL3_IL7,alpha12.D7.CD4ISP.IL3_IL7,alpha12.D7.CD3N.IL3_IL7,...
            alpha12.(p.days).PROT.IL7_TNFa, alpha12.(p.days).CD4ISP.IL7_TNFa, alpha12.(p.days).CD3N.IL7_TNFa,...
            alpha12.(p.days).PROT.TNFa_CXCL12, alpha12.(p.days).CD4ISP.TNFa_CXCL12, alpha12.(p.days).CD3N.TNFa_CXCL12];


    data_alpha21_D7 = [alpha21.D7.PROT.SCF_IL7,alpha21.D7.CD4ISP.SCF_IL7,alpha21.D7.CD3N.SCF_IL7,...
        alpha21.D7.PROT.IL3_IL7,alpha21.D7.CD4ISP.IL3_IL7,alpha21.D7.CD3N.IL3_IL7,...
        alpha21.(p.days).PROT.IL7_TNFa, alpha21.(p.days).CD4ISP.IL7_TNFa, alpha21.(p.days).CD3N.IL7_TNFa,...
            alpha21.(p.days).PROT.TNFa_CXCL12, alpha21.(p.days).CD4ISP.TNFa_CXCL12, alpha21.(p.days).CD3N.TNFa_CXCL12];

res.alpha12_d7_mean = mean(log10(data_alpha12_D7));
res.alpha12_d7_std = std(log10(data_alpha12_D7));
res.alpha21_d7_mean = mean(log10(data_alpha21_D7));
res.alpha21_d7_std = std(log10(data_alpha21_D7));

%% Alpha values D14-28 

%Set alpha values in 'struct' object

p.days = "D14";
p.cell = "CD3N";

    p.int = "SCF_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL3_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL7_TNFa";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));    
    

p.cell = "CD3P";

    p.int = "SCF_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL3_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL7_TNFa";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    
p.cell = "SP";

    p.int = "SCF_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL3_IL7";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));

    p.int = "IL7_TNFa";
    alpha12.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2)-1);
    alpha21.(p.days).(p.cell).(p.int) = parameters_MC.(p.days).(p.cell).(p.int)(:,size(parameters_MC.(p.days).(p.cell).(p.int),2));
    

    data_alpha12_D14 = [alpha12.(p.days).CD3N.SCF_IL7,alpha12.(p.days).CD3P.SCF_IL7,alpha12.(p.days).SP.SCF_IL7,...
            alpha12.(p.days).CD3N.IL3_IL7,alpha12.(p.days).CD3P.IL3_IL7,alpha12.(p.days).SP.IL3_IL7,...
            alpha12.(p.days).CD3N.IL7_TNFa, alpha12.(p.days).CD3P.IL7_TNFa, alpha12.(p.days).SP.IL7_TNFa];

    data_alpha21_D14 = [alpha21.(p.days).CD3N.SCF_IL7,alpha21.(p.days).CD3P.SCF_IL7,alpha21.(p.days).SP.SCF_IL7,...
        alpha21.(p.days).CD3N.IL3_IL7,alpha21.(p.days).CD3P.IL3_IL7,alpha21.(p.days).SP.IL3_IL7,...
        alpha21.(p.days).CD3N.IL7_TNFa, alpha21.(p.days).CD3P.IL7_TNFa, alpha21.(p.days).SP.IL7_TNFa];

res.alpha12_d14_mean = mean(log10(data_alpha12_D14));
res.alpha12_d14_std = std(log10(data_alpha12_D14));
res.alpha21_d14_mean = mean(log10(data_alpha21_D14));
res.alpha21_d14_std = std(log10(data_alpha21_D14));

save("alpha_mean_std.mat","res")

%% Create figures for differentiation stage

p.days = "D7";

    %SCF-IL7
    p.int = "SCF_IL7";
    p.cyt1 = "SCF";
    p.cyt2 = "IL7";
    p.start = 1;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 3;%index value in the data_alpha matrix corresponding to this interaction
    p.alpha12 = 1; %for alpha12
    p.ylim = [-3 3];
    
    figureAlpha_fun(p,data_alpha12_D7(:,p.start:p.end))

    p.alpha12 = 0; %for alpha21
    p.ylim = [-2 2];
    figureAlpha_fun(p,data_alpha21_D7(:,p.start:p.end))

    %IL3-IL7
    p.int = "IL3_IL7";
    p.cyt1 = "IL3";
    p.cyt2 = "IL7";
    p.start = 4;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 6;%index value in the data_alpha matrix corresponding to this interaction
    p.alpha12 = 1; %for alpha12
    p.ylim = [-5 5];
    
    figureAlpha_fun(p,data_alpha12_D7(:,p.start:p.end))

    p.alpha12 = 0; %for alpha21
    p.ylim = [-2 2];
    figureAlpha_fun(p,data_alpha21_D7(:,p.start:p.end))

    %IL7-TNFa
    p.int = "IL7_TNFa";
    p.cyt1 = "IL7";
    p.cyt2 = "TNFa";
    p.start = 7;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 9;%index value in the data_alpha matrix corresponding to this interaction
    p.alpha12 = 1; %for alpha12
    p.ylim = [-5 5];
    
    figureAlpha_fun(p,data_alpha12_D7(:,p.start:p.end))

    p.alpha12 = 0; %for alpha21
    p.ylim = [-5 5];
    figureAlpha_fun(p,data_alpha21_D7(:,p.start:p.end))

    %TNFa-CXCL12
    p.int = "TNFa_CXCL12";
    p.cyt1 = "TNFa";
    p.cyt2 = "CXCL12";
    p.start = 10;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 12;%index value in the data_alpha matrix corresponding to this interaction
    p.alpha12 = 1; %for alpha12
    p.ylim = [-4 4];
    
    figureAlpha_fun(p,data_alpha12_D7(:,p.start:p.end))

    p.alpha12 = 0; %for alpha21
    p.ylim = [-3 3];
    figureAlpha_fun(p,data_alpha21_D7(:,p.start:p.end))

%% Create figures for maturation stage

p.days = "D14";

    %SCF-IL7
    p.int = "SCF_IL7";
    p.cyt1 = "SCF";
    p.cyt2 = "IL7";
    p.start = 1;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 3;%index value in the data_alpha matrix corresponding to this interaction
    p.alpha12 = 1; %for alpha12
    p.ylim = [-3 3];
    
    figureAlpha_fun(p,data_alpha12_D14(:,p.start:p.end))

    p.alpha12 = 0; %for alpha21
    p.ylim = [-5 5];
    figureAlpha_fun(p,data_alpha21_D14(:,p.start:p.end))

    %IL3-IL7
    p.int = "IL3_IL7";
    p.cyt1 = "IL3";
    p.cyt2 = "IL7";
    p.start = 4;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 6;%index value in the data_alpha matrix corresponding to this interaction
    p.alpha12 = 1; %for alpha12
    p.ylim = [-5 5];
    
    figureAlpha_fun(p,data_alpha12_D14(:,p.start:p.end))

    p.alpha12 = 0; %for alpha21
    p.ylim = [-4 4];
    figureAlpha_fun(p,data_alpha21_D14(:,p.start:p.end))

    %IL7-TNFa
    p.int = "IL7_TNFa";
    p.cyt1 = "IL7";
    p.cyt2 = "TNFa";
    p.start = 7;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 9;%index value in the data_alpha matrix corresponding to this interaction
    p.alpha12 = 1; %for alpha12
    p.ylim = [-4 4];
    
    figureAlpha_fun(p,data_alpha12_D14(:,p.start:p.end))

    p.alpha12 = 0; %for alpha21
    p.ylim = [-5 5];
    figureAlpha_fun(p,data_alpha21_D14(:,p.start:p.end))
