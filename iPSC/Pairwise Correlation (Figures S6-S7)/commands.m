%% Load data

clear all

% Cytokine data
cytokinesTable = readtable("cytConc.xlsx", "ReadVariableNames", true);
cytokinesNames = cytokinesTable.Properties.VariableNames;
cytokinesConc = table2array(cytokinesTable);
logConc = log10(cytokinesConc);


load("MCMC_samples.mat")



%% Calculate posterior parameter correlations

corrFolder = fullfile("Figures", "PosteriorCorrelations");
if ~exist(corrFolder, "dir")
    mkdir(corrFolder)
end

posteriorCorr = struct();


hexColors = [
    "#292f40"
    "#3e475d"
    "#556177"
    "#6d7d91"
    "#8699ab"
    "#a1b7c4"
    "#c0d5dd"
    "#f6ebe3"
    "#eec6bd"
    "#dda49a"
    "#c6857a"
    "#ab675e"
    "#8f4c45"
    "#72332e"
    "#551b19"
];

anchorColors = zeros(numel(hexColors), 3);

for q = 1:numel(hexColors)
    hexValue = char(extractAfter(hexColors(q), "#"));

    anchorColors(q,:) = [
        hex2dec(hexValue(1:2)), ...
        hex2dec(hexValue(3:4)), ...
        hex2dec(hexValue(5:6))
    ] / 255;
end

cmap = interp1( ...
    linspace(-1, 1, size(anchorColors,1)), ...
    anchorColors, ...
    linspace(-1, 1, 256), ...
    "linear");

daysName = ["D7", "D14"];

cellType1 = ["PROT", "CD4ISP", "CD3N"];
cellType2 = ["CD3N", "CD3P", "SP"];


allCytokineNames = cytokinesNames;

color5Parameters = "#a781a2";  % Purple background
color7Parameters = "#ec9f79";  % Orange background

for i = 1:numel(daysName)

    days = daysName(i);
    daysField = char(days);

    if days == "D7"
        cellTypeNames = cellType1;
    
        cellTypeDisplayNames = [
            "ProT"
            "CD4ISP (CD3−)"
            "DP (CD3−)"
        ];
    
        cytokineList = allCytokineNames;
        dayTitle = "Days 0–7";
    else
        cellTypeNames = cellType2;
    
        cellTypeDisplayNames = [
            "DP (CD3−)"
            "DP (CD3+)"
            "8SP (CD3+)"
        ];
    
        cytokineList = allCytokineNames(1:5);
        dayTitle = "Days 7–21 ";
    end

    nRows = numel(cytokineList);
    nCols = numel(cellTypeNames);
    
    figWidth = 22;
    figHeight = 5*nRows + 5;

    fig = figure( ...
        "Visible", "off", ...
        "Units", "centimeters", ...
        "Position", [0, 0, figWidth, figHeight]);    
    tl = tiledlayout(fig, nRows, nCols, ...
        "TileSpacing", "tight", ...
        "Padding", "compact");

    colormap(fig, cmap)

    axLast = [];

    for j = 1:numel(cellTypeNames)

        cellTypeName = cellTypeNames(j);
        cellTypeField = char(cellTypeName);

        cellTypeDisplayName = cellTypeDisplayNames(j);

        for c = 1:numel(cytokineList)

            cytokineName = cytokineList{c};

            mcmc = samples.(daysField)...
                .(cellTypeField).(cytokineName);

            R = corr(mcmc, ...
                "Type", "Pearson", ...
                "Rows", "complete");

            posteriorCorr.(daysField)...
                .(cellTypeField).(cytokineName) = R;

            nParameters = size(mcmc,2);
            if nParameters == 5
                paramNames = {
                    '$\mathbf{E_0}$', ...
                    '$\mathbf{EC_{50}}$', ...
                    '$\mathbf{h}$', ...
                    '$\mathbf{E_1}$', ...
                    '$\mathbf{\sigma}$'
                };
                matrixBackground = color5Parameters;
            elseif nParameters == 7
                paramNames = {
                    '$\mathbf{E_0}$', ...
                    '$\mathbf{E_1}$', ...
                    '$\mathbf{EC_{50,1}}$', ...
                    '$\mathbf{h_1}$', ...
                    '$\mathbf{EC_{50,2}}$', ...
                    '$\mathbf{h_2}$', ...
                    '$\mathbf{\sigma}$'
                };
                matrixBackground = color7Parameters;
            end


            tileNumber = (c-1)*nCols + j;
            ax = nexttile(tl, tileNumber);
            axLast = ax;

            imagesc(ax, R)

            ax.Color = matrixBackground;
            
            backgroundMargin = 0.18;
            
            xlim(ax, [
                0.5 - backgroundMargin, ...
                nParameters + 0.5 + backgroundMargin
            ]);
            
            ylim(ax, [
                0.5 - backgroundMargin, ...
                nParameters + 0.5 + backgroundMargin
            ]);
            
            axis(ax, "square")
            caxis(ax, [-1 1])
            ax.Layer = "top";

            xticks(ax, 1:numel(paramNames))
            yticks(ax, 1:numel(paramNames))
            xticklabels(ax, paramNames)
            yticklabels(ax, paramNames)

            ax.TickLabelInterpreter = "latex";
            ax.XTickLabelRotation = 45;
            ax.FontSize = 10;
            
            ax.XTickLabelRotation = 45;
    
            if c == 1
                title(ax, ...
                    {char(cellTypeDisplayName)}, ...
                    "Interpreter", "none", ...
                    "FontSize", 12, ...
                    "FontWeight", "bold");
            end
            
            if j == 1
                yl = ylabel(ax, cytokineName, ...
                    "Interpreter", "none", ...
                    "FontSize", 12, ...
                    "FontWeight", "bold");
            
                yl.Rotation = 0;
                yl.HorizontalAlignment = "right";
                yl.VerticalAlignment = "middle";
            end
                end
    end

    title(tl, dayTitle, ...
    "FontSize", 14, ...
    "FontWeight", "bold");

    cb = colorbar(axLast);
    cb.Layout.Tile = "south";
    cb.Orientation = "horizontal";
    cb.Label.String = "Posterior correlation";
    cb.FontSize = 9;
    cb.Ticks = -1:0.2:1;


    fileName = sprintf("%s_posterior_correlations.pdf", days);

    exportgraphics(fig, ...
        fullfile(corrFolder, fileName), ...
        "ContentType", "vector");

    close(fig)
end

save("MCMC_correlations.mat", "posteriorCorr")
