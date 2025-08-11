% This code will generate data points for the single dose-responses as well
% as generate noise for these points based on the standard errors in Yale
% et al., 2022.

clear variables

%Load cytokine concentrations and set names
cytConc = readtable('../../Data/cytConc.xlsx');
cytokines = cytConc.Properties.VariableNames;
cytConc = table2array(cytConc);

%Scaled concentration of cytokines (from Yale et al., 2022)
scaled_conc = [-2.366, -1, 0, 1, 2.366];

%Load polynomial coefficients and standard errors
load('../../Data/iPSC_SingleInt_Coeff.mat')%Polynomial coefficients
load('../../Data/iPSC_SingleInt_Stde.mat')%Standard errors

%% Define function handles for polynomials describing the single dose-responses
%This section defines the function handles for the polynomials

y.D7.PROT.SCF = @(x) coeff.D7.PROT.SCF(1) + coeff.D7.PROT.SCF(2)*x + coeff.D7.PROT.SCF(3)*x.^2;
y.D7.PROT.Flt3L = @(x) coeff.D7.PROT.Flt3L(1) + coeff.D7.PROT.Flt3L(2)*x + coeff.D7.PROT.Flt3L(3)*x.^2;
y.D7.PROT.IL3 = @(x) coeff.D7.PROT.IL3(1) + coeff.D7.PROT.IL3(2)*x + coeff.D7.PROT.IL3(3)*x.^2;
y.D7.PROT.IL7 =  @(x) coeff.D7.PROT.IL7(1) + coeff.D7.PROT.IL7(2)*x + coeff.D7.PROT.IL7(3)*x.^2 + coeff.D7.PROT.IL7(4)*x.^3;
y.D7.PROT.TNFa = @(x) coeff.D7.PROT.TNFa(1) + coeff.D7.PROT.TNFa(2)*x + coeff.D7.PROT.TNFa(3)*x.^2;
y.D7.PROT.CXCL12 = @(x) coeff.D7.PROT.CXCL12(1) + coeff.D7.PROT.CXCL12(2)*x;

y.D7.CD4ISP.SCF = @(x) coeff.D7.CD4ISP.SCF(1) + coeff.D7.CD4ISP.SCF(2)*x + coeff.D7.CD4ISP.SCF(3)*x.^2;
y.D7.CD4ISP.Flt3L = @(x) coeff.D7.CD4ISP.Flt3L(1) + coeff.D7.CD4ISP.Flt3L(2)*x + coeff.D7.CD4ISP.Flt3L(3)*x.^2;
y.D7.CD4ISP.IL3 = @(x) coeff.D7.CD4ISP.IL3(1) + coeff.D7.CD4ISP.IL3(2)*x + coeff.D7.CD4ISP.IL3(3)*x.^2;
y.D7.CD4ISP.IL7 = @(x) coeff.D7.CD4ISP.IL7(1) + coeff.D7.CD4ISP.IL7(2)*x + coeff.D7.CD4ISP.IL7(3)*x.^2 + coeff.D7.CD4ISP.IL7(4)*x.^3;
y.D7.CD4ISP.TNFa = @(x) coeff.D7.CD4ISP.TNFa(1) + coeff.D7.CD4ISP.TNFa(2)*x + coeff.D7.CD4ISP.TNFa(3)*x.^2;
y.D7.CD4ISP.CXCL12 = @(x) coeff.D7.CD4ISP.CXCL12(1) + coeff.D7.CD4ISP.CXCL12(2)*x;

y.D7.CD3N.SCF = @(x) coeff.D7.CD3N.SCF(1) + coeff.D7.CD3N.SCF(2)*x + coeff.D7.CD3N.SCF(3)*x.^2;
y.D7.CD3N.Flt3L = @(x) coeff.D7.CD3N.Flt3L(1) + coeff.D7.CD3N.Flt3L(2)*x + coeff.D7.CD3N.Flt3L(3)*x.^2;
y.D7.CD3N.IL3 = @(x) coeff.D7.CD3N.IL3(1) + coeff.D7.CD3N.IL3(2)*x + coeff.D7.CD3N.IL3(3)*x.^2;
y.D7.CD3N.IL7 = @(x) coeff.D7.CD3N.IL7(1) + coeff.D7.CD3N.IL7(2)*x + coeff.D7.CD3N.IL7(3)*x.^2 + coeff.D7.CD3N.IL7(4)*x.^3;
y.D7.CD3N.TNFa = @(x) coeff.D7.CD3N.TNFa(1) + coeff.D7.CD3N.TNFa(2)*x + coeff.D7.CD3N.TNFa(3)*x.^2;
y.D7.CD3N.CXCL12 = @(x) coeff.D7.CD3N.CXCL12(1) + coeff.D7.CD3N.CXCL12(2)*x;

y.D14.CD3N.SCF = @(x) coeff.D14.CD3N.SCF(1) + coeff.D14.CD3N.SCF(2)*x + coeff.D14.CD3N.SCF(3)*x.^2 + coeff.D14.CD3N.SCF(4)*x.^3;
y.D14.CD3N.Flt3L = @(x) coeff.D14.CD3N.Flt3L(1) + coeff.D14.CD3N.Flt3L(2)*x + coeff.D14.CD3N.Flt3L(3)*x.^2;
y.D14.CD3N.IL3 = @(x) coeff.D14.CD3N.IL3(1) + coeff.D14.CD3N.IL3(2)*x;
y.D14.CD3N.IL7 = @(x) coeff.D14.CD3N.IL7(1) + coeff.D14.CD3N.IL7(2)*x + coeff.D14.CD3N.IL7(3)*x.^2 + coeff.D14.CD3N.IL7(4)*x.^3;
y.D14.CD3N.TNFa = @(x) coeff.D14.CD3N.TNFa(1) + coeff.D14.CD3N.TNFa(2)*x + coeff.D14.CD3N.TNFa(3)*x.^2 + coeff.D14.CD3N.TNFa(4)*x.^3;

y.D14.CD3P.SCF = @(x) coeff.D14.CD3P.SCF(1) + coeff.D14.CD3P.SCF(2)*x + coeff.D14.CD3P.SCF(3)*x.^2 + coeff.D14.CD3P.SCF(4)*x.^3;
y.D14.CD3P.Flt3L = @(x) coeff.D14.CD3P.Flt3L(1) + coeff.D14.CD3P.Flt3L(2)*x + coeff.D14.CD3P.Flt3L(3)*x.^2;
y.D14.CD3P.IL3 = @(x) coeff.D14.CD3P.IL3(1) + coeff.D14.CD3P.IL3(2)*x;
y.D14.CD3P.IL7 = @(x) coeff.D14.CD3P.IL7(1) + coeff.D14.CD3P.IL7(2)*x + coeff.D14.CD3P.IL7(3)*x.^2 + coeff.D14.CD3P.IL7(4)*x.^3;
y.D14.CD3P.TNFa = @(x) coeff.D14.CD3P.TNFa(1) + coeff.D14.CD3P.TNFa(2)*x + coeff.D14.CD3P.TNFa(3)*x.^2 + coeff.D14.CD3P.TNFa(4)*x.^3;

y.D14.SP.SCF = @(x) coeff.D14.SP.SCF(1) + coeff.D14.SP.SCF(2)*x + coeff.D14.SP.SCF(3)*x.^2 + coeff.D14.SP.SCF(4)*x.^3;
y.D14.SP.Flt3L = @(x) coeff.D14.SP.Flt3L(1) + coeff.D14.SP.Flt3L(2)*x + coeff.D14.SP.Flt3L(3)*x.^2;
y.D14.SP.IL3 = @(x) coeff.D14.SP.IL3(1) + coeff.D14.SP.IL3(2)*x;
y.D14.SP.IL7 = @(x) coeff.D14.SP.IL7(1) + coeff.D14.SP.IL7(2)*x + coeff.D14.SP.IL7(3)*x.^2 + coeff.D14.SP.IL7(4)*x.^3;
y.D14.SP.TNFa = @(x) coeff.D14.SP.TNFa(1) + coeff.D14.SP.TNFa(2)*x + coeff.D14.SP.TNFa(3)*x.^2 + coeff.D14.SP.TNFa(4)*x.^3;


%% SCF 
k=1;

p.name = num2str(cell2mat(cytokines(k)));

%DAYS 7-14
p.days = 'D7';

    %PRO-T CELLS
    p.cell = 'PROT';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD4ISP CELLS
    p.cell = 'CD4ISP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3N CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

%DAYS 14-28
p.days = 'D14';

    %CD3- CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3+ CELLS
    p.cell = 'CD3P';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %8SP CELLS
    p.cell = 'SP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

%% FLT3L
k=2;

p.name = num2str(cell2mat(cytokines(k)));

%DAYS 7-14
p.days = 'D7';

    %PRO-T CELLS
    p.cell = 'PROT';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD4ISP CELLS
    p.cell = 'CD4ISP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3N CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

%DAYS 14-28
p.days = 'D14';

    %CD3- CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3+ CELLS
    p.cell = 'CD3P';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %8SP CELLS
    p.cell = 'SP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);
%% IL3
k=3;

p.name = num2str(cell2mat(cytokines(k)));

%DAYS 7-14
p.days = 'D7';

    %PRO-T CELLS
    p.cell = 'PROT';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD4ISP CELLS
    p.cell = 'CD4ISP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3N CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

%DAYS 14-28
p.days = 'D14';

    %CD3- CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3+ CELLS
    p.cell = 'CD3P';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %8SP CELLS
    p.cell = 'SP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);
%% IL7
k=4;

p.name = num2str(cell2mat(cytokines(k)));

%DAYS 7-14
p.days = 'D7';

    %PRO-T CELLS
    p.cell = 'PROT';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD4ISP CELLS
    p.cell = 'CD4ISP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3N CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

%DAYS 14-28
p.days = 'D14';

    %CD3- CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3+ CELLS
    p.cell = 'CD3P';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %8SP CELLS
    p.cell = 'SP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);
%% TNFa
k=5;

p.name = num2str(cell2mat(cytokines(k)));

%DAYS 7-14
p.days = 'D7';

    %PRO-T CELLS
    p.cell = 'PROT';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD4ISP CELLS
    p.cell = 'CD4ISP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3N CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

%DAYS 14-28
p.days = 'D14';

    %CD3- CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3+ CELLS
    p.cell = 'CD3P';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %8SP CELLS
    p.cell = 'SP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);
%% CXCL12
k=6;

p.name = num2str(cell2mat(cytokines(k)));

%DAYS 7-14
p.days = 'D7';

    %PRO-T CELLS
    p.cell = 'PROT';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD4ISP CELLS
    p.cell = 'CD4ISP';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

    %CD3N CELLS
    p.cell = 'CD3N';
    [noise.(p.days).(p.cell).(p.name),data.(p.days).(p.cell).(p.name)] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde);

%% Save results
save("Results/Noise_data.mat","noise")
save("Results/Data_values.mat","data")

