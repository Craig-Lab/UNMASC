%% Load data

clear all

% Cytokine data
cytokinesTable = readtable("cytokines.xlsx", "ReadRowNames", true);
cytokinesNames = cytokinesTable.Properties.RowNames;
cytokinesConc = table2array(cytokinesTable);
logConc = log10(cytokinesConc);

% Results from single fits
load("MCMC_params.mat");
load("MCMC_ci.mat");

% Data points of surfaces
vials = 2:4;

files = {
    'SCF_IL7', 'CD7pCD5n', 'SCF-IL7_CD7+CD5-_vial%d.xlsx'
    'SCF_IL7', 'CD7pCD5p', 'SCF-IL7_CD7+CD5+_vial%d.xlsx'
    'IL3_IL7', 'CD7pCD5n', 'IL3-IL7_CD7+CD5-_vial%d.xlsx'
    'IL3_IL7', 'CD7pCD5p', 'IL3-IL7_CD7+CD5+_vial%d.xlsx'
};

for v = vials
    vialName = sprintf('vial%d', v);
    
    for i = 1:size(files,1)
        conditionName = files{i,1};
        cellTypeName = files{i,2};
        fileName = sprintf(files{i,3}, v);

        filePath = fullfile("Data", fileName);
        
        dataNotNorm = table2array(readtable(filePath));
        data.(vialName).(conditionName).(cellTypeName) = ...    
            dataNotNorm./(max(dataNotNorm(:)));
    end
end

%Remove outlier from data
data.("vial2").("IL3_IL7").("CD7pCD5n")(1,2) = NaN;
data.("vial2").("IL3_IL7").("CD7pCD5p")(2,2) = NaN;

rng(1)

%% Do surface fits

n = 200; %Number of times to do Monte Carlo

fprintf('Computing surface fits \n');
nbVials = 4;
cellInt = ["SCF_IL7", "IL3_IL7"];
cellType = ["CD7pCD5n", "CD7pCD5p"];

for i = 2:nbVials
    vialName = sprintf('vial%d', i);

    for j = 1:numel(cellInt)
        conditionName = cellInt(j);

        for k = 1:numel(cellType)
            cellTypeName = cellType(k);
            fprintf('Processing %s | %s | %s\n', vialName, conditionName, cellTypeName);

            parts = split(conditionName, "_");
            cytokine1 = parts(1);
            cytokine2 = parts(2);
            
            %Get first cytokine concentration
            idx1 = find(strcmp(cytokinesNames, cytokine1));
            x1 = cytokinesConc(idx1,1:end)';


            %Get second cytokine concentration
            idx2 = find(strcmp(cytokinesNames, cytokine2));
            x2 = cytokinesConc(idx2,1:end)';

            %Get data points for surface
            y = data.(vialName).(conditionName).(cellTypeName)(1:end,1:end);

            %Get parameters from single fits
            params1 = mcmcParams.(vialName).(conditionName).(cellTypeName).(cytokine1);
            params2 = mcmcParams.(vialName).(conditionName).(cellTypeName).(cytokine2);

            %Get confidence intervals from single fits
            ci1 = ci95.(vialName).(conditionName).(cellTypeName).(cytokine1);
            ci2 = ci95.(vialName).(conditionName).(cellTypeName).(cytokine2);

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

            [fit_parameters,resnorm, residual,exitflag,output,iterationData] = ...
                fit_surface(x1,x2,y,params1,params2,ci1,ci2,p);

            surfaceFitsParams.(vialName).(conditionName).(cellTypeName) = fit_parameters;
            surfaceMonotonicity.(vialName).(conditionName).(cellTypeName) = p.monotonic;
            resNorm.(vialName).(conditionName).(cellTypeName) = resnorm;

            N = size(y,1)*size(y,2); 
            resnorm_normalized = resnorm/N;
            rms = sqrt(resnorm/N);

            [CI,parameters_MC] = monteCarlo_fun(x1,x2,p,fit_parameters,rms,n,ci1,ci2);

            surfacesParamsMC.(vialName).(conditionName).(cellTypeName) = parameters_MC;
            surfacesParamsMC_CI.(vialName).(conditionName).(cellTypeName) = CI;
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

nbVials = 4;
cellInt = ["SCF_IL7", "IL3_IL7"];
cellType = ["CD7pCD5n", "CD7pCD5p"];

load("surfaceFitsParams.mat")
load("surfaceMonotonicity.mat")


p.map = load("ownColor.mat"); %load own color map for surface

for i = 2:nbVials
    vialName = sprintf('vial%d', i);

    for j = 1:numel(cellInt)
        conditionName = cellInt(j);

        for k = 1:numel(cellType)
            cellTypeName = cellType(k);
            fprintf('Plotting surface %s | %s | %s\n', vialName, conditionName, cellTypeName);

             % Get data points
            dataPts = data.(vialName).(conditionName).(cellTypeName)(2:end, 2:end);

            parts = split(conditionName, "_");
            cytokine1 = parts(1);
            cytokine2 = parts(2);

            p.drug1_name = cytokine1; 
            p.drug2_name = cytokine2; 

            % Original drug concentrations
            idx1 = find(strcmp(cytokinesNames, cytokine1));
            idx2 = find(strcmp(cytokinesNames, cytokine2));
            
            drug1 = cytokinesConc(idx1,2:end)';
            drug2 = cytokinesConc(idx2,2:end)';

            param = surfaceFitsParams.(vialName).(conditionName).(cellTypeName);
            p.monotonic = surfaceMonotonicity.(vialName).(conditionName).(cellTypeName);

            
            [fig,E] = figure_fun(drug1, drug2, dataPts, param, p);

            %Define file name
            file_name = strcat("Figures/",vialName,"_", conditionName,"_", cellTypeName,".pdf");

            exportgraphics(fig,file_name)
        end 
    end 
end 



