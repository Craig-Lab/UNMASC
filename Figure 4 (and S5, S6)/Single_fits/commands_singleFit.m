%This code will perform the fits of the single dose-responses.

clear variables

%Load cytokine concentrations and set names
cytConc = readtable('../../Data/cytConc.xlsx');
cytokines = cytConc.Properties.VariableNames;
cytConc = table2array(cytConc);

%Scaled concentration of cytokines
scaled_conc = [-2.366, -1, 0, 1, 2.366];

noise.all = load('Noise_data.mat'); % all noise data (not filtered by monoticity)

load("noiseMonotonicDecreasing.mat"); %monoDec
load("noiseMonotonicIncreasing.mat"); %monoInc
load("noiseMultiphasic.mat"); %multi

noise.monoDec = monoDec;
noise.monoInc = monoInc;
noise.multi = multi;

load('Data_values.mat')%data points

% Set up lower and upper bound limits
% E0, C, h, Emax
p.lb1 = [0 -inf 0 0];
p.ub1 = [1.05 inf 5 1.05];

p.lb2 = [0 -inf -5 0];
p.ub2 = [1.05 inf 0 1.05]; 

%% SCF 
p.k = 1; %first cytokine

%Get name of cytokine
p.name = num2str(cell2mat(cytokines(p.k)));

% DAYS 7-14
%Get days 
p.days = 'D7';

    %Get cell type
    p.cell = 'PROT';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 0.5;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2.25;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        %Perform fit
        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD4ISP';

      %Initial guesses for the fit
        p.E01 = 0.5;
        p.C1 = 0.5;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 3;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3N';

      %Initial guesses for the fit
       p.E01 = 0.5;
       p.C1 = 0.5;
       p.h1 = 2;
       p.Emax1 = 1;

       p.E02 = 0.5;
       p.C2 = 2;
       p.h2 = -2;
       p.Emax2 = 1;

       p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
       p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

% DAYS 14-28    
%Get days 
p.days = 'D14';

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 0.5;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2.25;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3P';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 0;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'SP';

        %Initial guesses for the fit
        p.E01 = 0.7;
        p.C1 = 0.7;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0.7;
        p.C2 = 2.25;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

 %% FLT3L
p.k = 2; 

%Get name of cytokine
p.name = num2str(cell2mat(cytokines(p.k)));

% DAYS 7-14
%Get days 
p.days = 'D7';

    %Get cell type
    p.cell = 'PROT';

        %Initial guesses for the fit
        p.E01 = 0.5;
        p.C1 = 0.2;
        p.h1 = 3;
        p.Emax1 = 1;

        p.E02 = 0.5;
        p.C2 = 2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD4ISP';

        %Initial guesses for the fit
        p.E01 = 0.5;
        p.C1 = 0.2;
        p.h1 = 3;
        p.Emax1 = 1;

        p.E02 = 0.5;
        p.C2 = 2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 0.5;
        p.C1 = 0.2;
        p.h1 = 3;
        p.Emax1 = 1;

        p.E02 = 0.5;
        p.C2 = 2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

% DAYS 14-28 
%Get days 
p.days = 'D14';

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 0.5;
        p.C1 = 0.2;
        p.h1 = 3;
        p.Emax1 = 1;

        p.E02 = 0.5;
        p.C2 = 2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3P';

        %Initial guesses for the fit
        p.E01 = 0.3;
        p.C1 = 0.2;
        p.h1 = 3;
        p.Emax1 = 1;

        p.E02 = 0.3;
        p.C2 = 2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'SP';

        %Initial guesses for the fit
        p.E01 = 0.8;
        p.C1 = 0;
        p.h1 = 3;
        p.Emax1 = 1;

        p.E02 = 0.7;
        p.C2 = 2;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

 %% IL3
p.k = 3; 

%Get name of cytokine
p.name = num2str(cell2mat(cytokines(p.k)));

% DAYS 7-14
%Get days 
p.days = 'D7';

    %Get cell type
    p.cell = 'PROT';

        %Initial guesses for the fit
        p.E01 = 0.5;
        p.C1 = -1;
        p.h1 = 1.5;
        p.Emax1 = 1;

        p.E02 = 0.5;
        p.C2 = 1;
        p.h2 = -1.5;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD4ISP';

        %Initial guesses for the fit
        p.E01 = 0.7;
        p.C1 = -1;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0.7;
        p.C2 = 1;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 0.5;
        p.C1 = -1.5;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0.2;
        p.C2 = 1.5;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

% DAYS 14-28
%Get days 
p.days = 'D14';

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit (Monotonic increasing)
        p.E01 = 0.5;
        p.C1 = 0;
        p.h1 = 1;
        p.Emax1 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3P';

        %Initial guesses for the fit (Monotonic increasing)
        p.E01 = 0.5;
        p.C1 = 0;
        p.h1 = 1;
        p.Emax1 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'SP';

        %Initial guesses for the fit (Monotonic decreasing)
        p.E01 = 1;
        p.C1 = 0;
        p.h1 = 1;
        p.Emax1 = 0.9;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

        % % Adjust lower and upper bound limits
        % % E0, C, h, Emax
        % p.lb1 = [0 -inf -5 0];
        % p.ub1 = [1 inf 0 1.05];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

%% IL7
p.k = 4; 

%Get name of cytokine
p.name = num2str(cell2mat(cytokines(p.k)));

%Adjust lower and upper bound limits
% E0, C, h, Emax
p.lb1 = [0 -inf 0 0];
p.ub1 = [1 inf 5 1.05];

% DAYS 7-14
%Get days 
p.days = 'D7';

    %Get cell type
    p.cell = 'PROT';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 1;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2.2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD4ISP';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 1;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2.5;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 1;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2.5;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

% DAYS 14-28
%Get days 
p.days = 'D14';

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 1;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2.2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3P';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 1;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 2.2;
        p.h2 = -3;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'SP';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = 1;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 3;
        p.h2 = -1;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

%% TNFa
p.k=5;

%Get name of cytokine
p.name = num2str(cell2mat(cytokines(p.k)));

% DAYS 7-14
%Get days 
p.days = 'D7';

    %Get cell type
    p.cell = 'PROT';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = -2;
        p.h1 = 1;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 0.5;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD4ISP';

        %Initial guesses for the fit
        p.E01 = 0;
        p.C1 = -2;
        p.h1 = 2;
        p.Emax1 = 1;

        p.E02 = 0;
        p.C2 = 1;
        p.h2 = -2;
        p.Emax2 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [p.E02, p.C2, p.h2, p.Emax2];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 1;
        p.C1 = 0.5;
        p.h1 = 2;
        p.Emax1 = 0;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

% DAYS 14-28
%Get days 
p.days = 'D14';

    %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit (Monotonic decreasing)
        p.E01 = 1;
        p.C1 = -0.25;
        p.h1 = 1;
        p.Emax1 = 0;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD3P';

        %Initial guesses for the fit (Monotonic decreasing)
        p.E01 = 1;
        p.C1 = -1;
        p.h1 = 1;
        p.Emax1 = 0;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'SP';

        %Initial guesses for the fit (Monotonic decreasing)
        p.E01 = 1;
        p.C1 = 0;
        p.h1 = 1;
        p.Emax1 = 0.3;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

        [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

%% CXCL12
p.k=6;

%Get name of cytokine
p.name = num2str(cell2mat(cytokines(p.k)));

% DAYS 7-14 (We only have CXCL12 for these days)
%Get days 
p.days = 'D7';

    %Get cell type
    p.cell = 'PROT';

        %Initial guesses for the fit
        p.E01 = 0.8;
        p.C1 = 0.5;
        p.h1 = 1;
        p.Emax1 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

    %Get cell type
    p.cell = 'CD4ISP';

        %Initial guesses for the fit
        p.E01 = 0.9;
        p.C1 = 0.5;
        p.h1 = 1;
        p.Emax1 = 1;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

     %Get cell type
    p.cell = 'CD3N';

        %Initial guesses for the fit
        p.E01 = 1;
        p.C1 = 1;
        p.h1 = 0.1;
        p.Emax1 = 0.98;

        p.param_guess1 = [p.E01, p.C1, p.h1, p.Emax1];
        p.param_guess2 = [];

       [results.(p.days).(p.cell).(p.name),noiseParam.(p.days).(p.cell).(p.name),resnorm.(p.days).(p.cell).(p.name).one,resnorm.(p.days).(p.cell).(p.name).two] = single_fit(p,noise,data,cytConc);

%% Save noise parameters
save("Results/noiseParameters","noiseParam")

%% Create tables with results
% 6 tables are created
% 3 tables for Days 7-14 and 3 tables for Days 14-28

% DAYS 7 -14
p.days = 'D7';

sz = [24 6]; % 6 cytokines for Days 7-14

varTypes = ["double", "double", "double", "double", "double", "double"];
varNames = ["SCF", "FLT3L", "IL3", "IL7", "TNFa", "CXCL12"];

rowNames = ["E01", "E01 lb", "E01 ub", "C1", "C1 lb", "C1 ub", "h1", "h1 lb", "h1 ub", "Emax1", "Emax1 lb", "Emax1 ub",...
    "E02", "E02 lb", "E02 ub", "C2", "C2 lb", "C2 ub", "h2", "h2 lb", "h2 ub", "Emax2", "Emax2 lb", "Emax2 ub"];

    % PROT CELLS
    p.cell = 'PROT';

    resTable.(p.days).(p.cell) = table_fun(p,cytokines,sz,varTypes,varNames,rowNames,results);

    % CD4ISP CELLS
    p.cell = 'CD4ISP';
    
    resTable.(p.days).(p.cell) = table_fun(p,cytokines,sz,varTypes,varNames,rowNames,results);

    % CD3N CELLS
    p.cell = 'CD3N';
    
    resTable.(p.days).(p.cell) = table_fun(p,cytokines,sz,varTypes,varNames,rowNames,results);

% DAYS 14-28
p.days = 'D14';

sz = [24 5]; 

varTypes = ["double", "double", "double", "double", "double"];
varNames = ["SCF", "FLT3L", "IL3", "IL7", "TNFa"];

rowNames = ["E01", "E01 lb", "E01 ub", "C1", "C1 lb", "C1 ub", "h1", "h1 lb", "h1 ub", "Emax1", "Emax1 lb", "Emax1 ub",...
    "E02", "E02 lb", "E02 ub", "C2", "C2 lb", "C2 ub", "h2", "h2 lb", "h2 ub", "Emax2", "Emax2 lb", "Emax2 ub"];

    % CD3N CELLS
    p.cell = 'CD3N';

    resTable.(p.days).(p.cell) = table_fun(p,cytokines,sz,varTypes,varNames,rowNames,results);

    % CD4ISP CELLS
    p.cell = 'CD3P';
    
    resTable.(p.days).(p.cell) = table_fun(p,cytokines,sz,varTypes,varNames,rowNames,results);

    % SP CELLS
    p.cell = 'SP';
    
    resTable.(p.days).(p.cell) = table_fun(p,cytokines,sz,varTypes,varNames,rowNames,results);

filename = "single_fits_ipscs.xlsx";
writetable(resTable.D7.PROT,filename,'Sheet','D7 PROT','WriteRowNames',true);
writetable(resTable.D7.CD4ISP,filename,'Sheet','D7 CD4ISP','WriteRowNames',true);
writetable(resTable.D7.CD3N,filename,'Sheet','D7 CD3N','WriteRowNames',true);
writetable(resTable.D14.CD3N,filename,'Sheet','D14 CD3N','WriteRowNames',true);
writetable(resTable.D14.CD3P,filename,'Sheet','D14 CD3P','WriteRowNames',true);
writetable(resTable.D14.SP,filename,'Sheet','D14 SP','WriteRowNames',true);

