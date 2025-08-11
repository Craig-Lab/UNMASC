clear variables

load("fit_parameters.mat")
load("parameters_MC.mat")
 
n = 200;

%% Compute beta values : DAYS 7-14

p.days = "D7";

%SCF-IL7
p.int = "SCF_IL7";
    
    %PROT cells
    p.cell = "PROT";
    for i = 1:n
        beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %CD4ISP cells
    p.cell = "CD4ISP";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %DP (CD3-) cells
    p.cell = "CD3N";
    for i = 1:n
        beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 


%IL3-IL7
p.int = "IL3_IL7";
    
    %PROT cells
    p.cell = "PROT";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %CD4ISP cells
    p.cell = "CD4ISP";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %DP (CD3-) cells
    p.cell = "CD3N";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

%IL7-TNFa
p.int = "IL7_TNFa";
    
    %PROT cells
    p.cell = "PROT";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %CD4ISP cells
    p.cell = "CD4ISP";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %DP (CD3-) cells
    p.cell = "CD3N";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

%TNFa-CXCL12
p.int = "TNFa_CXCL12";
    
    %PROT cells
    p.cell = "PROT";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %CD4ISP cells
    p.cell = "CD4ISP";
    for i = 1:n
        beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %DP (CD3-) cells
    p.cell = "CD3N";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    data_beta_D7 = [beta.(p.days).PROT.SCF_IL7',beta.(p.days).CD4ISP.SCF_IL7',beta.(p.days).CD3N.SCF_IL7',...
            beta.(p.days).PROT.IL3_IL7',beta.(p.days).CD4ISP.IL3_IL7',beta.(p.days).CD3N.IL3_IL7',...
            beta.(p.days).PROT.IL7_TNFa',beta.(p.days).CD4ISP.IL7_TNFa',beta.(p.days).CD3N.IL7_TNFa',...
            beta.(p.days).PROT.TNFa_CXCL12', beta.(p.days).CD4ISP.TNFa_CXCL12',beta.(p.days).CD3N.TNFa_CXCL12'];


res.beta_d7_mean = mean(log10(data_beta_D7));
res.beta_d7_std = std(log10(data_beta_D7));

%% Compute beta values : DAYS 14-28

p.days = "D14";

%SCF-IL7
p.int = "SCF_IL7";
    
    %PROT cells
    p.cell = "CD3N";
    for i = 1:n
        beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %CD4ISP cells
    p.cell = "CD3P";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %DP (CD3-) cells
    p.cell = "SP";
    for i = 1:n
        beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 


%IL3-IL7
p.int = "IL3_IL7";
    
    %PROT cells
    p.cell = "CD3N";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %CD4ISP cells
    p.cell = "CD3P";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %DP (CD3-) cells
    p.cell = "SP";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

%IL7-TNFa
p.int = "IL7_TNFa";
    
    %PROT cells
    p.cell = "CD3N";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %CD4ISP cells
    p.cell = "CD3P";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 

    %DP (CD3-) cells
    p.cell = "SP";
    for i = 1:n
         beta.(p.days).(p.cell).(p.int)(i) = beta_fun(parameters_MC.(p.days).(p.cell).(p.int)(i,1:4));
    end 


    data_beta_D14 = [beta.(p.days).CD3N.SCF_IL7',beta.(p.days).CD3P.SCF_IL7',beta.(p.days).SP.SCF_IL7',...
            beta.(p.days).CD3N.IL3_IL7',beta.(p.days).CD3P.IL3_IL7',beta.(p.days).SP.IL3_IL7',...
            beta.(p.days).CD3N.IL7_TNFa',beta.(p.days).CD3P.IL7_TNFa',beta.(p.days).SP.IL7_TNFa'];

res.beta_d14_mean = mean(log10(data_beta_D14));
res.beta_d14_std = std(log10(data_beta_D14));

save("beta_mean_std.mat","res")

%% Create figures for differentiation stage beta

p.days = "D7";

    %SCF-IL7
    p.int = "SCF_IL7";
    p.cyt1 = "SCF";
    p.cyt2 = "IL7";
    p.start = 1;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 3;%index value in the data_alpha matrix corresponding to this interaction
    p.ylim = [-1.5 1.5];
    
    figureBeta_fun(p,data_beta_D7(:,p.start:p.end))

    %IL3-IL7
    p.int = "IL3_IL7";
    p.cyt1 = "IL3";
    p.cyt2 = "IL7";
    p.start = 4;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 6;%index value in the data_alpha matrix corresponding to this interaction
    p.ylim = [-1 1];
    
    figureBeta_fun(p,data_beta_D7(:,p.start:p.end))

    %IL7-TNFa
    p.int = "IL7_TNFa";
    p.cyt1 = "IL7";
    p.cyt2 = "TNFa";
    p.start = 7;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 9;%index value in the data_alpha matrix corresponding to this interaction
    p.ylim = [-3 3];
    
    figureBeta_fun(p,data_beta_D7(:,p.start:p.end))

    %TNFa_CXCL12
    p.int = "TNFa_CXCL12";
    p.cyt1 = "TNFa";
    p.cyt2 = "CXCL12";
    p.start = 10;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 12;%index value in the data_alpha matrix corresponding to this interaction
    p.ylim = [-5 5];
    
    figureBeta_fun(p,data_beta_D7(:,p.start:p.end))

%% Create figures for maturation stage beta

p.days = "D14";

    %SCF-IL7
    p.int = "SCF_IL7";
    p.cyt1 = "SCF";
    p.cyt2 = "IL7";
    p.start = 1;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 3;%index value in the data_alpha matrix corresponding to this interaction
    p.ylim = [-2 2];
    
    figureBeta_fun(p,data_beta_D14(:,p.start:p.end))

    %IL3-IL7
    p.int = "IL3_IL7";
    p.cyt1 = "IL3";
    p.cyt2 = "IL7";
    p.start = 4;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 6;%index value in the data_alpha matrix corresponding to this interaction
    p.ylim = [-5 5];
    
    figureBeta_fun(p,data_beta_D14(:,p.start:p.end))

    %IL7-TNFa
    p.int = "IL7_TNFa";
    p.cyt1 = "IL7";
    p.cyt2 = "TNFa";
    p.start = 7;%index value in the data_alpha matrix corresponding to this interaction
    p.end = 9;%index value in the data_alpha matrix corresponding to this interaction
    p.ylim = [-3 3];
    
    figureBeta_fun(p,data_beta_D14(:,p.start:p.end))

%% Function

function beta = beta_fun(param)

    E0 = param(1);
    E1 = param(2);
    E2 = param(3);
    E3 = param(4);

    beta = (E0-max(E1,E2))./(max(E1,E2)-E3);
end 
