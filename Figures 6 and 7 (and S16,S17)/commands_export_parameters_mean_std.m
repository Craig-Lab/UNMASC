%% Load data
clear variables

alpha = load("Alpha/alpha_mean_std.mat");
beta = load("Beta/beta_mean_std.mat");

%% Create table

%DAYS 7-14

sz = [12 6]; % 4 combinations for Days 7-14, 14 parameters max

varTypes = ["double", "double","double", "double","double", "double"];
varNames = ["Mean alpha12", "Standard Deviation alpha12", ...
    "Mean alpha21", "Standard Deviation alpha21", ...
    "Mean beta", "Standard Deviation beta"];

rowNames = ["ProT SCF-IL7", "ProT IL3-IL7", "ProT IL7-TNFa","ProT TNFa-CXCL12", ...
    "CD4ISP SCF-IL7", "CD4ISP IL3-IL7", "CD4ISP IL7-TNFa","CD4ISP TNFa-CXCL12",...
    "CD3N SCF-IL7", "CD3N IL3-IL7", "CD3N IL7-TNFa","CD3N TNFa-CXCL12"];

tb1 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames,'RowNames',rowNames);

tb1(:,1) = array2table(alpha.res.alpha12_d7_mean');
tb1(:,2) = array2table(alpha.res.alpha12_d7_std');
tb1(:,3) = array2table(alpha.res.alpha21_d7_mean'); 
tb1(:,4) = array2table(alpha.res.alpha21_d7_std'); 
tb1(:,5) = array2table(beta.res.beta_d7_mean'); 
tb1(:,6) = array2table(beta.res.beta_d7_std'); 

%DAYS 14-28

sz = [9 6]; % 4 combinations for Days 7-14, 14 parameters max

varTypes = ["double", "double","double", "double","double", "double"];
varNames = ["Mean alpha12", "Standard Deviation alpha12", ...
    "Mean alpha21", "Standard Deviation alpha21", ...
    "Mean beta", "Standard Deviation beta"];

rowNames = ["CD3N SCF-IL7", "CD3N IL3-IL7", "CD3N IL7-TNFa" ...
    "CD3P SCF-IL7", "CD3P IL3-IL7", "CD3P IL7-TNFa"...
    "8SP SCF-IL7", "8SP IL3-IL7", "8SP IL7-TNFa"];

tb2 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames,'RowNames',rowNames);

tb2(:,1) = array2table(alpha.res.alpha12_d14_mean');
tb2(:,2) = array2table(alpha.res.alpha12_d14_std');
tb2(:,3) = array2table(alpha.res.alpha21_d14_mean'); 
tb2(:,4) = array2table(alpha.res.alpha21_d14_std'); 
tb2(:,5) = array2table(beta.res.beta_d14_mean'); 
tb2(:,6) = array2table(beta.res.beta_d14_std'); 

%% Export tables

filename = "alpha_beta_mean_std.xlsx";

writetable(tb1,filename,'Sheet','D7','WriteRowNames',true);
writetable(tb2,filename,'Sheet','D14','WriteRowNames',true);