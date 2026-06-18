%% Load data

clear all

% Cytokine data
cytokinesTable = readtable("cytConc.xlsx", "ReadVariableNames", true);
cytokinesNames = cytokinesTable.Properties.VariableNames;
cytokinesConc = table2array(cytokinesTable);
logConc = log10(cytokinesConc);

% Load data points
load("Data_values.mat")

rng(1)

%% Computing MCMC for single fits

fprintf('Computing MCMC for single fits \n');

daysName = ["D7", "D14"];

cellType1 = ["PROT", "CD4ISP", "CD3N"];
cellType2 = ["CD3N", "CD3P", "SP"];

for i = 1:numel(daysName)
    days = daysName(i);

    if strcmp(days,"D7") == 1
        cellTypeNames = cellType1;
    elseif strcmp(days,"D14") == 1
        cellTypeNames = cellType2;
        cytokinesNames = cytokinesNames(1:5);
    end 

    for j = 1:numel(cellTypeNames)
        cellType = cellTypeNames(j);

        for k = 1:length(cytokinesNames)
            cytokine = cytokinesNames{k};
            fprintf('Processing %s | %s | %s\n', days, cellType, cytokine);

            idx = find(strcmp(cytokinesNames, cytokine));
            x = cytokinesConc(:,idx)';
            y = data.(days).(cellType).(cytokine)(:);

            [samples.(days).(cellType).(cytokine),...
                ci95.(days).(cellType).(cytokine)] = single_fit(x',y);
        end 
    end 
end 

%% Save results of fit

cytokinesNames = cytokinesTable.Properties.VariableNames;

for i = 1:numel(daysName)
    days = daysName(i);

    if strcmp(days,"D7") == 1
        cellTypeNames = cellType1;
    elseif strcmp(days,"D14") == 1
        cellTypeNames = cellType2;
        cytokinesNames = cytokinesNames(1:5);
    end 

    for j = 1:numel(cellTypeNames)
        cellType = cellTypeNames(j);

        for k = 1:length(cytokinesNames)
            cytokine = cytokinesNames{k};

            params = median(samples.(days).(cellType).(cytokine));

            mcmcParams.(days).(cellType).(cytokine) = params;
        end 
    end 
end 


save('MCMC_params.mat', 'mcmcParams'); %save parameters
save('MCMC_ci.mat','ci95') %save confidence intervals
save("MCMC_samples.mat", "samples")

%% Create tables with MCMC parameter estimates and confidence intervals

load("MCMC_params.mat")   % load MCMC_params
load("MCMC_ci.mat")       % load MCMC_ci

resTable = struct();

rowNames = ["E0", "E0 lb", "E0 ub", ...
            "Emax1", "Emax1 lb", "Emax1 ub", ...
            "C1", "C1 lb", "C1 ub", ...
            "h1", "h1 lb", "h1 ub", ...
            "C2", "C2 lb", "C2 ub", ...
            "h2", "h2 lb", "h2 ub"];


% DAYS 7-14
p.days = 'D7';

cytokines_D7 = ["SCF", "Flt3L", "IL3", "IL7", "TNFa", "CXCL12"];

varNames = cytokines_D7;
varTypes = repmat("double", 1, numel(varNames));


sz = [numel(rowNames), numel(varNames)];

cells_D7 = ["PROT", "CD4ISP", "CD3N"];

for c = 1:numel(cells_D7)
    p.cell = cells_D7(c);

    resTable.(p.days).(p.cell) = table_fun_mcmc( ...
        p, cytokines_D7, sz, varTypes, varNames, rowNames, ...,
            mcmcParams, ci95);
end


% DAYS 14-28
p.days = 'D14';

cytokines_D14 = ["SCF", "Flt3L", "IL3", "IL7", "TNFa"];

varNames = cytokines_D14;
varTypes = repmat("double", 1, numel(varNames));

sz = [numel(rowNames), numel(varNames)];

cells_D14 = ["CD3N", "CD3P", "SP"];

for c = 1:numel(cells_D14)
    p.cell = cells_D14(c);

    resTable.(p.days).(p.cell) = table_fun_mcmc( ...
        p, cytokines_D14, sz, varTypes, varNames, rowNames, ...
            mcmcParams, ci95);
end


filename = "single_fits_ipscs_MCMC.xlsx";

writetable(resTable.D7.PROT,   filename, 'Sheet', 'D7 PROT',   'WriteRowNames', true);
writetable(resTable.D7.CD4ISP, filename, 'Sheet', 'D7 CD4ISP', 'WriteRowNames', true);
writetable(resTable.D7.CD3N,   filename, 'Sheet', 'D7 CD3N',   'WriteRowNames', true);

writetable(resTable.D14.CD3N, filename, 'Sheet', 'D14 CD3N', 'WriteRowNames', true);
writetable(resTable.D14.CD3P, filename, 'Sheet', 'D14 CD3P', 'WriteRowNames', true);
writetable(resTable.D14.SP,   filename, 'Sheet', 'D14 SP',   'WriteRowNames', true);


%% Plot results

load("MCMC_params.mat")
load("MCMC_ci.mat")
load("MCMC_samples.mat")

fprintf('Making figures\n');

marker = "#de8a7d";
line = "#555e7d";
markerSize = 50;
figureSize = [0 0 2 2]; 
figureUnits = "centimeters";

nbins = 100;

daysName = ["D7", "D14"];

cellType1 = ["PROT", "CD4ISP", "CD3N"];
cellType2 = ["CD3N", "CD3P", "SP"];

for i = 1:numel(daysName)
    days = daysName(i);

    if strcmp(days,"D7") == 1
        cellTypeNames = cellType1;
    elseif strcmp(days,"D14") == 1
        cellTypeNames = cellType2;
        cytokinesNames = cytokinesNames(1:5);
    end 

    for j = 1:numel(cellTypeNames)
        cellType = cellTypeNames(j);

        for k = 1:length(cytokinesNames)
            cytokine = cytokinesNames{k};
            

            idx = find(strcmp(cytokinesNames, cytokine));
            x = log10(cytokinesConc(:,idx)');
            y = data.(days).(cellType).(cytokine)(:);

            mcmc = samples.(days).(cellType).(cytokine);

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
                xlabel(['log_{10}(' char(cytokine) ')'])
                ylabel("Effect")
                fontsize(16,"points")
                leg=legend("Fit","Data");
                leg.Location = "northeastoutside";
                ylim([0 1.1])
                axis square


                file_name = strcat('Figures/', days, '_', ...
                                   cellType, '_', cytokine, '.pdf');
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
                xlabel(['log_{10}(' char(cytokine) ')'])
                ylabel("Effect")
                fontsize(16,"points")
                leg=legend("Fit","Data");
                leg.Location = "northeastoutside";
                ylim([0 1.1])
                axis square

                file_name = strcat('Figures/', days, '_', ...
                                   cellType, '_', cytokine, '.pdf');
                exportgraphics(fig, file_name, 'ContentType', 'vector');
            end 
        end 
    end 
end 

%% Functions
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


function tb = table_fun_mcmc(p, cytokines, sz, varTypes, varNames, rowNames, ...
                             MCMC_params, MCMC_ci)

    tb = table('Size', sz, ...
               'VariableTypes', varTypes, ...
               'VariableNames', varNames, ...
               'RowNames', rowNames);

    for i = 1:numel(cytokines)

        p.name = char(cytokines(i));

        valsRaw = MCMC_params.(p.days).(p.cell).(p.name);
        cisRaw  = MCMC_ci.(p.days).(p.cell).(p.name);

        valsRaw = valsRaw(:);

        valsRaw = valsRaw(1:end-1);

        if size(cisRaw, 1) == 2
            cisRaw = cisRaw(:, 1:end-1);
            cisStd = cisRaw;

        elseif size(cisRaw, 2) == 2
            cisRaw = cisRaw(1:end-1, :);
            cisStd = cisRaw';

        else
            error("For %s %s %s: confidence interval has unexpected size.", ...
                  p.days, p.cell, p.name);
        end

        nRaw = numel(valsRaw);

        % [E0, Emax1, C1, h1, C2, h2]
        vals = nan(6, 1);
        cis  = nan(2, 6);

        if nRaw == 6
            %[E0, Emax1, C1, h1, C2, h2]

            vals(:) = valsRaw(:);
            cis(:, :) = cisStd(:, 1:6);

        elseif nRaw == 4
            % Monotonic single fit:
            % [E0, C1, h1, Emax1]
            %
            % Non-monotonic single fit:
            % [E0, Emax1, C1, h1, C2, h2]

            vals(1) = valsRaw(1);   % E0
            vals(2) = valsRaw(4);   % Emax1
            vals(3) = valsRaw(2);   % C1
            vals(4) = valsRaw(3);   % h1
            vals(5) = NaN;          % C2 missing
            vals(6) = NaN;          % h2 missing

            cis(:,1) = cisStd(:,1); % E0
            cis(:,2) = cisStd(:,4); % Emax1
            cis(:,3) = cisStd(:,2); % C1
            cis(:,4) = cisStd(:,3); % h1
            cis(:,5) = [NaN; NaN];  % C2 missing
            cis(:,6) = [NaN; NaN];  % h2 missing

        end

        new_vec = nan(3*6, 1);

        k = 1;
        for j = 1:6

            new_vec(k)   = vals(j);
            new_vec(k+1) = cis(1,j);
            new_vec(k+2) = cis(2,j);

            k = k + 3;

        end

        tb(:, i) = array2table(new_vec);

    end

end


