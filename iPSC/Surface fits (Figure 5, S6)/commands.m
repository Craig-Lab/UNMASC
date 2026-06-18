%% Load data

clear all

% Cytokine data
cytokinesTable = readtable("cytConc.xlsx", "ReadVariableNames", true);
cytokinesNames = cytokinesTable.Properties.VariableNames;
cytokinesConc = table2array(cytokinesTable);

% Load coefficients from polynomials
load("iPSC_Surface_Coeff.mat")

% Load monotonic structure for interactions
load("monotonic.mat")

% Load data values for surfaces
load("dataSurfaces.mat")

% Results from single fits
load("MCMC_params.mat");
load("MCMC_ci.mat");

rng(1)

%% Do surface fits

load("E_guesses.mat")

n = 200; %Number of times to do Monte Carlo

fprintf('Computing surface fits \n');

daysName = ["D7", "D14"];

cellType1 = ["PROT", "CD4ISP", "CD3N"];
cellType2 = ["CD3N", "CD3P", "SP"];
cellIntNames1 = ["SCF_IL7", "IL3_IL7", "IL7_TNFa", "TNFa_CXCL12"];
cellIntNames2 = ["SCF_IL7", "IL3_IL7", "IL7_TNFa"];

for i = 1:numel(daysName)
    days = daysName(i);

    if strcmp(days,"D7") == 1
        cellTypeNames = cellType1;
        cellIntNames = cellIntNames1;
    elseif strcmp(days,"D14") == 1
        cellTypeNames = cellType2;
        cytokinesNames = cytokinesNames(1:5);
        cellIntNames = cellIntNames2; 
    end 

    for j = 1:numel(cellTypeNames)
        cellType = cellTypeNames(j);

        for k = 1:numel(cellIntNames)
            cellInt = cellIntNames(k);

            fprintf('Processing %s | %s | %s\n', days, cellType, cellInt);

            parts = split(cellInt, "_");
            cytokine1 = parts(1);
            cytokine2 = parts(2);
            
            %Get first cytokine concentration
            idx1 = find(strcmp(cytokinesNames, cytokine1));
            x1 = cytokinesConc(:,idx1);

            %Get second cytokine concentration
            idx2 = find(strcmp(cytokinesNames, cytokine2));
            x2 = cytokinesConc(:,idx2);

            %Get data points for surface
            y = data.(days).(cellType).(cellInt)(1:end, 1:end);

            %Get parameters from single fits
            params1 = mcmcParams.(days).(cellType).(cytokine1);
            params2 = mcmcParams.(days).(cellType).(cytokine2);

            %Get confidence intervals from single fits
            ci1 = ci95.(days).(cellType).(cytokine1);
            ci2 = ci95.(days).(cellType).(cytokine2);

            %Assess monotonicity of single fits
            if (length(params1) == 7) && (length(params2) == 7) %both drugs are non-monotonic
                p.monotonic = 0;
            elseif (length(params1) == 5) && (length(params2)==7) %first drug is monotonic
                p.monotonic = 1;
            elseif (length(params1) == 7) && (length(params2)==5) %second drug is monotonic
                p.monotonic = 2;
            elseif (length(params1) == 5) && (length(params2)==5) %both drugs are monotonic
                p.monotonic = 3;
            end 

            p.guesses = E_guesses.(days).(cellType).(cellInt);

            [fit_parameters,resnorm, residual,exitflag,output,iterationData] = ...
                fit_surface(x1,x2,y,params1,params2,ci1,ci2,p);

            surfaceFitsParams.(days).(cellType).(cellInt) = fit_parameters;
            surfaceMonotonicity.(days).(cellType).(cellInt) = p.monotonic;
            resNorm.(days).(cellType).(cellInt) = resnorm;

            N = size(y,1)*size(y,2); 
            resnorm_normalized = resnorm/N;
            rms = sqrt(resnorm/N);

            [CI,parameters_MC] = monteCarlo_fun(x1,x2,p,fit_parameters,rms,n,ci1,ci2);

            surfacesParamsMC.(days).(cellType).(cellInt) = parameters_MC;
            surfacesParamsMC_CI.(days).(cellType).(cellInt) = CI;
        end 
    end 
end 

%% Save results

save("surfacesParamsMC.mat", "surfacesParamsMC")
save("surfacesParamsMC_CI.mat", "surfacesParamsMC_CI")
save("surfaceFitsParams.mat", "surfaceFitsParams")
save("surfaceMonotonicity.mat", "surfaceMonotonicity")
save("surfacesResNorm.mat", "resNorm")

%% Plot surface fits

load("surfacesParamsMC.mat")
load("surfaceFitsParams.mat")
load("surfaceMonotonicity.mat")

daysName = ["D7", "D14"];

cellType1 = ["PROT", "CD4ISP", "CD3N"];
cellType2 = ["CD3N", "CD3P", "SP"];
cellIntNames1 = ["SCF_IL7", "IL3_IL7", "IL7_TNFa", "TNFa_CXCL12"];
cellIntNames2 = ["SCF_IL7", "IL3_IL7", "IL7_TNFa"];

p.map = load("ownColor.mat"); %load own color map for surface

cytokinesNames = cytokinesTable.Properties.VariableNames;


for i = 1:numel(daysName)
    days = daysName(i);

    if strcmp(days,"D7") == 1
        cellTypeNames = cellType1;
        cellIntNames = cellIntNames1;
        cytNames = cytokinesNames;
    elseif strcmp(days,"D14") == 1
        cellTypeNames = cellType2;
        cytNames = cytokinesNames(1:5);
        cellIntNames = cellIntNames2; 
    end 

    for j = 1:numel(cellTypeNames)
        cellType = cellTypeNames(j);

        for k = 1:numel(cellIntNames)
            cellInt = cellIntNames(k);

            fprintf('Processing %s | %s | %s\n', days, cellType, cellInt);

            dataPts = data.(days).(cellType).(cellInt)(1:end, 1:end);

            parts = split(cellInt, "_");
            cytokine1 = parts(1);
            cytokine2 = parts(2);

            p.drug1_name = cytokine1; 
            p.drug2_name = cytokine2; 

            % Original drug concentrations
            idx1 = find(strcmp(cytNames, cytokine1));
            idx2 = find(strcmp(cytNames, cytokine2));
            
            drug1 = cytokinesConc(:,idx1);
            drug2 = cytokinesConc(:,idx2);
            
            param = surfaceFitsParams.(days).(cellType).(cellInt);
            p.monotonic = surfaceMonotonicity.(days).(cellType).(cellInt);

            [fig,E] = figure_fun(drug1, drug2, dataPts, param, p);

            %Define file name
            file_name = strcat("Figures/",days,"_", cellType,"_", cellInt,".pdf");

            exportgraphics(fig,file_name)
        end 
    end 
end 


