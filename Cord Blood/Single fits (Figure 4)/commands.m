clear all 

%% Prepare data

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
        
        data.(vialName).(conditionName).(cellTypeName) = ...
            table2array(readtable(filePath));
    end
end

% Cytokine data
cytokinesTable = readtable("cytokines.xlsx", "ReadRowNames", true);
cytokinesNames = cytokinesTable.Properties.RowNames;
cytokinesConc = table2array(cytokinesTable);
cytokinesConc(:,1) = 1e-6;
logConc = log10(cytokinesConc);

rng(1)


%% Perform single fits

fprintf('Computing MCMC for single fits \n');
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

            % Fit and MCMC for 1st cytokine
            idx1 = find(strcmp(cytokinesNames, cytokine1));
            x1 = cytokinesConc(idx1,:)';
            y1 = data.(vialName).(conditionName).(cellTypeName)(:,1);
            y1 = y1./max(y1);

            [samples.(vialName).(conditionName).(cellTypeName).(cytokine1),...
                ci95.(vialName).(conditionName).(cellTypeName).(cytokine1)] = single_fit(x1,y1);


            % Fit and MCMC for 2nd cytokine
            idx2 = find(strcmp(cytokinesNames, cytokine2));
            x2 = cytokinesConc(idx2,:);
            y2 = data.(vialName).(conditionName).(cellTypeName)(1,:);
            y2 = y2./max(y2);

            [samples.(vialName).(conditionName).(cellTypeName).(cytokine2),...
               ci95.(vialName).(conditionName).(cellTypeName).(cytokine2)] = single_fit(x2,y2);

        end
    end
end

save("MCMC_samples.mat", "samples")

%% Save results of fit

for i = 2:nbVials
    vialName = sprintf('vial%d', i);

    for j = 1:numel(cellInt)
        conditionName = cellInt(j);

        for k = 1:numel(cellType)
            cellTypeName = cellType(k);

            parts = split(conditionName, "_");
            cytokineList = [parts(1), parts(2)];

            for c = 1:numel(cytokineList)
                cytokineName = cytokineList(c);

                params = median(samples.(vialName).(conditionName).(cellTypeName).(cytokineName));

                mcmcParams.(vialName).(conditionName).(cellTypeName).(cytokineName) = params;
            end 
        end 
    end 
end 

save('MCMC_params.mat', 'mcmcParams'); %save parameters
save('MCMC_ci.mat','ci95') %save confidence intervals

%% Plot results

load("MCMC_params.mat")
load("MCMC_ci.mat")
load("MCMC_samples.mat")

fprintf('Making figures\n');

nbVials = 4;
cellInt = ["SCF_IL7", "IL3_IL7"];
cellType = ["CD7pCD5n", "CD7pCD5p"];

marker = "#de8a7d";
line = "#555e7d";
markerSize = 50;
figureSize = [0 0 3 3]; 
figureUnits = "centimeters";

nbins = 100;

for i = 2:nbVials
    vialName = sprintf('vial%d', i);

    for j = 1:numel(cellInt)
        conditionName = cellInt(j);

        for k = 1:numel(cellType)
            cellTypeName = cellType(k);
            fprintf('Plotting %s | %s | %s\n', vialName, conditionName, cellTypeName);

            parts = split(conditionName, "_");
            cytokineList = [parts(1), parts(2)];

            for c = 1:numel(cytokineList)
                cytokineName = cytokineList(c);

                idx = find(strcmp(cytokinesNames, cytokineName));
                x = log10(cytokinesConc(idx,2:end)');

                if c ==1 
                    y = data.(vialName).(conditionName).(cellTypeName)(:,1);
                    y = y ./ max(y);
                    y = y(2:end);
                    mcmc = samples.(vialName).(conditionName).(cellTypeName).(cytokineName);
                elseif c==2
                    y = data.(vialName).(conditionName).(cellTypeName)(1,:);
                    y = y ./ max(y);
                    y = y(2:end);
                    mcmc = samples.(vialName).(conditionName).(cellTypeName).(cytokineName);
                end 


                if size(mcmc,2) == 5   % monotonic
                    param = median(mcmc);
                    behavior = 'monotonic';

                    x_vec = linspace(x(1)-0.5, x(end)+0.5);
                    y_vec = hill(x_vec, param);

                    fig = figure("Visible","off");
                    f.Units = figureUnits; 
                    f.Position = figureSize; 

                    plot(x_vec, y_vec, "LineWidth", 3,"Color",line)
                    hold on
                    scatter(x, y, markerSize, "MarkerFaceColor",marker, "MarkerEdgeColor",marker)
                    hold off
                    xlabel(['log_{10}(' char(cytokineName) ')'])
                    ylabel("Effect")
                    fontsize(16,"points")
                    leg=legend("Fit","Data");
                    leg.Location = "northeastoutside";
                    ylim([0 1.1])
                    axis square


                    file_name = strcat('Figures/', vialName, '_', conditionName, '_', ...
                                       cellTypeName, '_', cytokineName, '.pdf');
                    exportgraphics(fig, file_name, 'ContentType', 'vector');


                elseif size(mcmc,2) == 7   % non-monotonic

                    behavior = 'non-monotonic';
                    
                    param1 = median(mcmc);

                    x_vec = linspace(x(1)-0.5, x(end)+0.5);
                    y_vec = hillNonmonotonic(x_vec, param1);

                    fig = figure("Visible","off");
                    f.Units = figureUnits; 
                    f.Position = figureSize; 

                    plot(x_vec, y_vec, "LineWidth", 3, "Color",line)
                    hold on
                    scatter(x, y, markerSize, "MarkerFaceColor",marker, "MarkerEdgeColor",marker)
                    hold off
                    xlabel(['log_{10}(' char(cytokineName) ')'])
                    ylabel("Effect")
                    fontsize(16,"points")
                    leg=legend("Fit","Data");
                    leg.Location = "northeastoutside";
                    ylim([0 1.1])
                    axis square

                    file_name = strcat('Figures/', vialName, '_', conditionName, '_', ...
                                       cellTypeName, '_', cytokineName, '.pdf');
                    exportgraphics(fig, file_name, 'ContentType', 'vector');
                end
            end
        end
    end
end

%% hill function
function y = hill(x,param)

    E01 = param(1);
    C1 = param(2);
    h1 = param(3);
    Emax1 = param(4);

    E = @(d) E01 + (Emax1-E01)./(1+10.^((C1-d)*h1));
    y = E(x);

end 

function y = hillNonmonotonic(x, param)

    E0 = param(1);
    E1 = param(2);
    C1 = param(3);
    h1 = param(4);   % negative
    C2 = param(5);
    h2 = param(6);   % positive

    H1 = 1 ./ (1 + 10.^((x - C1) .* h1));
    H2 = 1 ./ (1 + 10.^((x - C2) .* h2));

    y = E1 + (E0 - E1) .* H1 .* H2;
end

