%% Load data

clear all

% Cytokine data
cytokinesTable = readtable("cytokines.xlsx", "ReadRowNames", true);
cytokinesNames = cytokinesTable.Properties.RowNames;


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



vialNames = ["vial2", "vial3", "vial4"];

combinationNames = [
    "SCF_IL7"
    "IL3_IL7"
];

combinationDisplayNames = [
    "SCF-IL7"
    "IL3-IL7"
];

cellTypeNames = [
    "CD7pCD5n"
    "CD7pCD5p"
];

cellTypeDisplayNames = [
    "CD7^+CD5^-"
    "CD7^+CD5^+"
];

cytokinesByCombination = {
    ["SCF", "IL7"]
    ["IL3", "IL7"]
};

cytokineDisplayNames = {
    ["SCF", "IL7"]
    ["IL3", "IL7"]
};

color5Parameters = "#a781a2";  % Purple background
color7Parameters = "#ec9f79";  % Orange background



for v = 1:numel(vialNames)

    vialName  = vialNames(v);
    vialField = char(vialName);

    for k = 1:numel(combinationNames)

        combinationName  = combinationNames(k);
        combinationField = char(combinationName);
        combinationTitle = combinationDisplayNames(k);

        cytokineNames = cytokinesByCombination{k};
        cytokineLabels = cytokineDisplayNames{k};

        nRows = numel(cytokineNames);
        nCols = numel(cellTypeNames);

        fig = figure( ...
            "Visible", "off", ...
            "Units", "centimeters", ...
            "Position", [0 0 16 17]);

        tl = tiledlayout(fig, nRows, nCols, ...
            "TileSpacing", "compact", ...
            "Padding", "compact");

        colormap(fig, cmap)

        axLast = [];

        for c = 1:numel(cytokineNames)

            cytokineName  = cytokineNames(c);
            cytokineField = char(cytokineName);
            cytokineLabel = cytokineLabels(c);

            for j = 1:numel(cellTypeNames)

                cellTypeName  = cellTypeNames(j);
                cellTypeField = char(cellTypeName);
                cellTypeLabel = cellTypeDisplayNames(j);

                mcmc = samples.(vialField) ...
                    .(combinationField) ...
                    .(cellTypeField) ...
                    .(cytokineField);

                R = corr(mcmc, ...
                    "Type", "Pearson", ...
                    "Rows", "complete");


                posteriorCorr.(vialField) ...
                    .(combinationField) ...
                    .(cellTypeField) ...
                    .(cytokineField) = R;

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
            else
                error( ...
                    "Unexpected number of MCMC parameters for " + ...
                    vialName + ", " + combinationName + ", " + ...
                    cellTypeName + ", " + cytokineName + ".");
            end

                tileNumber = (c - 1)*nCols + j;
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
                ax.FontSize = 12;

                if c == 1
                    title(ax, ['$\mathrm{' char(cellTypeLabel) '}$'], ...
                        "Interpreter", "latex", ...
                        "FontSize", 11, ...
                        "FontWeight", "bold");
                end

                if j == 1
                    yl = ylabel(ax, cytokineLabel, ...
                        "Interpreter", "none", ...
                        "FontSize", 11, ...
                        "FontWeight", "bold");

                    yl.Rotation = 0;
                    yl.HorizontalAlignment = "right";
                    yl.VerticalAlignment = "middle";
                end
            end
        end

        donorNumber = v;

        title(tl, ...
            sprintf("Donor %d: %s", donorNumber, combinationTitle), ...
            "Interpreter", "none", ...
            "FontSize", 14, ...
            "FontWeight", "bold");

        cb = colorbar(axLast);
        cb.Layout.Tile = "south";
        cb.Orientation = "horizontal";
        cb.Label.String = "Posterior correlation";
        cb.FontSize = 9;
        cb.Ticks = -1:0.2:1;

        % Export figure
        fileName = sprintf( ...
            "%s_%s_posterior_correlations.pdf", ...
            vialName, combinationName);

        exportgraphics(fig, ...
            fullfile(corrFolder, fileName), ...
            "ContentType", "vector");

        close(fig)
    end
end

save("MCMC_correlations.mat", "posteriorCorr")
