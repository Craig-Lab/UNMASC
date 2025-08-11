clear variables

%Load cytokine concentrations and set names
cytConc = readtable('../Data/cytConc.xlsx');
cytokines = cytConc.Properties.VariableNames;
cytConc = table2array(cytConc);

%Load coefficients
load("../Data/iPSC_Surface_Coeff.mat") %name of structure loaded is 'coeff'

%Define a "results" structure where we store the results of the single fits
results.D7.PROT = readtable("single_fits_ipscs.xlsx",'Sheet','D7 PROT','ReadRowNames',true);
results.D7.CD4ISP = readtable("single_fits_ipscs.xlsx",'Sheet','D7 CD4ISP','ReadRowNames',true);
results.D7.CD3N = readtable("single_fits_ipscs.xlsx",'Sheet','D7 CD3N','ReadRowNames',true);
results.D14.CD3N = readtable("single_fits_ipscs.xlsx",'Sheet','D14 CD3N','ReadRowNames',true);
results.D14.CD3P = readtable("single_fits_ipscs.xlsx",'Sheet','D14 CD3P','ReadRowNames',true);
results.D14.SP = readtable("single_fits_ipscs.xlsx",'Sheet','D14 SP','ReadRowNames',true);

%Load color for surfaces
p.map = load("ownColor.mat");

% Scaled concentration of cytokines
p.scaled_conc = [-2.366,-1, 0, 1, 2.366];


%% DAYS 7-14 PRO-T CELLS SCF_IL7

p.monotonic = 0;

p.days = 'D7';
p.cell = 'PROT';

p.drug1_name = "SCF";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X) %Displays the interaction to follow where the computation is at

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.9;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);


%% DAYS 7-14 PRO-T CELLS IL3-IL7

p.monotonic = 0;

p.days = 'D7';
p.cell = 'PROT';

p.drug1_name = "IL3";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.9;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 PRO-T CELLS IL7-TNFa

p.monotonic = 0;

p.days = 'D7';
p.cell = 'PROT';

p.drug1_name = "IL7";
p.drug2_name = "TNFa";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.9;
p.E2 = 0.1;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 PRO-T CELLS TNFa-CXCL12

p.monotonic = 2;

p.days = 'D7';
p.cell = 'PROT';

p.drug1_name = "TNFa";
p.drug2_name = "CXCL12";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 0.9;
p.E1 = 1;
p.E2 = 0.8;
p.E3 = 0.5;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD4ISP CELLS SCF-IL7

p.monotonic = 0; 

p.days = 'D7';
p.cell = 'CD4ISP';

p.drug1_name = "SCF";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.7;
p.E2 = 0.1;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD4ISP CELLS IL3-IL7

p.monotonic = 0; 

p.days = 'D7';
p.cell = 'CD4ISP';

p.drug1_name = "IL3";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.1;
p.E2 = 0.9;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD4ISP CELLS IL7-TNFA

p.monotonic = 0; 

p.days = 'D7';
p.cell = 'CD4ISP';

p.drug1_name = "IL7";
p.drug2_name = "TNFa";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD4ISP CELLS TNFA-CXCL12

p.monotonic = 2; 

p.days = 'D7';
p.cell = 'CD4ISP';

p.drug1_name = "TNFa";
p.drug2_name = "CXCL12";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD3N CELLS SCF-IL7

p.monotonic = 0; 

p.days = 'D7';
p.cell = 'CD3N';

p.drug1_name = "SCF";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD3N CELLS IL3-IL7

p.monotonic = 0; 

p.days = 'D7';
p.cell = 'CD3N';

p.drug1_name = "IL3";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD3N CELLS IL7-TNFA

p.monotonic = 0; 

p.days = 'D7';
p.cell = 'CD3N';

p.drug1_name = "IL7";
p.drug2_name = "TNFa";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 7-14 CD3N CELLS TNFA-CXCL12

p.monotonic = 3; %Change this for better fit, altough TNFa is multiphasic alone

p.days = 'D7';
p.cell = 'CD3N';

p.drug1_name = "TNFa";
p.drug2_name = "CXCL12";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 CD3N CELLS SCF-IL7

p.monotonic = 0; 

p.days = 'D14';
p.cell = 'CD3N';

p.drug1_name = "SCF";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 CD3N CELLS IL3-IL7

p.monotonic = 1; 

p.days = 'D14';
p.cell = 'CD3N';

p.drug1_name = "IL3";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 CD3N CELLS IL7-TNFa

p.monotonic = 2; 

p.days = 'D14';
p.cell = 'CD3N';

p.drug1_name = "IL7";
p.drug2_name = "TNFa";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 CD3P CELLS SCF-IL7

p.monotonic = 0; 

p.days = 'D14';
p.cell = 'CD3P';

p.drug1_name = "SCF";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 CD3P CELLS IL3-IL7

p.monotonic = 1; 

p.days = 'D14';
p.cell = 'CD3P';

p.drug1_name = "IL3";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 CD3P CELLS IL7-TNFa

p.monotonic = 2; 

p.days = 'D14';
p.cell = 'CD3P';

p.drug1_name = "IL7";
p.drug2_name = "TNFa";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 SP CELLS SCF-IL7

p.monotonic = 0; 

p.days = 'D14';
p.cell = 'SP';

p.drug1_name = "SCF";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 SP CELLS IL3-IL7

p.monotonic = 1; 

p.days = 'D14';
p.cell = 'SP';

p.drug1_name = "IL3";
p.drug2_name = "IL7";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% DAYS 14-28 SP CELLS IL7-TNFA

p.monotonic = 2; 

p.days = 'D14';
p.cell = 'SP';

p.drug1_name = "IL7";
p.drug2_name = "TNFa";

X = strcat(p.days," ",p.cell," ",p.drug1_name," ",p.drug2_name);
disp(X)

%Guesses for the fit
p.E0 = 1;
p.E1 = 0.5;
p.E2 = 0.5;
p.E3 = 0;

%Initialize data
p = initialize_fun(p,cytokines,cytConc,coeff);

%Perform fit, Monte Carlo simulation and make figure
[fit_parameters.(p.days).(p.cell).(p.int_name),resnorm.(p.days).(p.cell).(p.int_name),...
    parameters_MC.(p.days).(p.cell).(p.int_name),...
    CI.(p.days).(p.cell).(p.int_name),output.(p.days).(p.cell).(p.int_name),...
    costFcnData.(p.days).(p.cell).(p.int_name)] = runCode_fun(p,results);

%% Save results

%Fit parameters
save("Results/fit_parameters.mat","fit_parameters")

%Residual norm
save("Results/resnorm.mat","resnorm")

%Monte Carlo parameters
save("Results/parameters_MC.mat","parameters_MC")

%Confidence intervals Monte Carlo
save("Results/CI.mat","CI")

%Cost function iterations
save("Results/costFcnData.mat","costFcnData")
