function [] = fun_Dist_SingleFits(p, noiseParam, singleFits, combinationFits)

    alphaVal = 0.5;

    colorSingle = [1.00 0.00 0.00];   % red
    colorComb   = [0.10 0.48 0.13];   % green
    colorNoise  = [0.41 0.39 0.78];   % purple/blue
    colorInt    = [0.75 0.75 0.75];   % grey

    fs = 10;
    mainTitleFS = 10;
    lineWidth = 1.2;
    axisLW = 0.9;

    outDir = "Figures";
    if ~exist(outDir, "dir")
        mkdir(outDir)
    end

    fileTag = string(p.days) + "_" + string(p.cell) + "_" + string(p.int);
    fileTag = replace(fileTag, "/", "_");

    filePDF = fullfile(outDir, fileTag + "_allParams.pdf");

    % Non-monotonic samples:
    % E0, E1, EC50_1, h1, EC50_2, h2, ...
    %
    % Monotonic samples:
    % E0, EC50_1, h1, E1, ...

    if p.drug1_mono == 0
        drug1_C1_col = 3;
        drug1_h1_col = 4;
        drug1_C2_col = 5;
        drug1_h2_col = 6;
    else
        drug1_C1_col = 2;
        drug1_h1_col = 3;
        drug1_C2_col = NaN;
        drug1_h2_col = NaN;
    end

    if p.drug2_mono == 0
        drug2_C1_col = 3;
        drug2_h1_col = 4;
        drug2_C2_col = 5;
        drug2_h2_col = 6;
    else
        drug2_C1_col = 2;
        drug2_h1_col = 3;
        drug2_C2_col = NaN;
        drug2_h2_col = NaN;
    end

    plotInfo = {};

    plotInfo(end+1,:) = {"C11", p.drug1, drug1_C1_col, ...
        "C1 lb", "C1 ub", "C1", "C11", false, "$\log_{10}(C_{11})$", true};

    if p.drug1_mono == 0
        plotInfo(end+1,:) = {"C12", p.drug1, drug1_C2_col, ...
            "C2 lb", "C2 ub", "C2", "C12", false, "$\log_{10}(C_{12})$", true};
    end

    plotInfo(end+1,:) = {"C21", p.drug2, drug2_C1_col, ...
        "C1 lb", "C1 ub", "C1", "C21", false, "$\log_{10}(C_{21})$", true};

    if p.drug2_mono == 0
        plotInfo(end+1,:) = {"C22", p.drug2, drug2_C2_col, ...
            "C2 lb", "C2 ub", "C2", "C22", false, "$\log_{10}(C_{22})$", true};
    end

    if p.drug1_mono == 1
        flip_h11 = true;
    else
        flip_h11 = false;
    end

    plotInfo(end+1,:) = {"h11", p.drug1, drug1_h1_col, ...
        "h1 lb", "h1 ub", "h1", "h11", flip_h11, "$h_{11}$", true};

    if p.drug1_mono == 0
        plotInfo(end+1,:) = {"h12", p.drug1, drug1_h2_col, ...
            "h2 lb", "h2 ub", "h2", "h12", false, "$h_{12}$", true};
    end

    if p.drug2_mono == 1
        flip_h21 = true;
    else
        flip_h21 = false;
    end

    plotInfo(end+1,:) = {"h21", p.drug2, drug2_h1_col, ...
        "h1 lb", "h1 ub", "h1", "h21", flip_h21, "$h_{21}$", true};

    if p.drug2_mono == 0
        plotInfo(end+1,:) = {"h22", p.drug2, drug2_h2_col, ...
            "h2 lb", "h2 ub", "h2", "h22", false, "$h_{22}$", true};
    end

    while size(plotInfo, 1) < 8
        plotInfo(end+1,:) = {"", "", NaN, "", "", "", "", false, "", false};
    end

    nParams = 8;

    % Figure

    fig = figure("Visible", "off", "Color", "w");

    figWidth = 11.0;
    figHeight = 2.3;

    fig.Units = "inches";
    fig.Position = [1, 1, figWidth, figHeight];

    fig.PaperUnits = "inches";
    fig.PaperPosition = [0, 0, figWidth, figHeight];
    fig.PaperSize = [figWidth, figHeight];


    annotation(fig, "textbox", [0.20 0.91 0.60 0.07], ...
        "String", string(p.int), ...
        "Interpreter", "none", ...
        "HorizontalAlignment", "center", ...
        "VerticalAlignment", "middle", ...
        "FontWeight", "bold", ...
        "FontSize", mainTitleFS, ...
        "EdgeColor", "none");


    leftMargin = 0.04;
    rightMargin = 0.01;
    bottomMargin = 0.22;
    axHeight = 0.62;
    axGap = 0.028;

    axWidth = (1 - leftMargin - rightMargin - (nParams - 1)*axGap) / nParams;

    for k = 1:nParams

        paramName    = plotInfo{k,1};
        drugName     = plotInfo{k,2};
        sampleCol    = plotInfo{k,3};
        lbRow        = plotInfo{k,4};
        ubRow        = plotInfo{k,5};
        singleRow    = plotInfo{k,6};
        combRow      = plotInfo{k,7};
        flipCombSign = plotInfo{k,8};
        xlab         = plotInfo{k,9};
        doPlot       = plotInfo{k,10};

        axLeft = leftMargin + (k-1)*(axWidth + axGap);

        ax = axes(fig, "Position", [axLeft bottomMargin axWidth axHeight]);
        hold(ax, "on")

        if ~doPlot
            axis(ax, "off")
            continue
        end


        sampleMat = noiseParam.(p.days).(p.cell).(drugName);

        if sampleCol > size(sampleMat, 2)
            error("For %s %s %s %s: sample column %d does not exist. Sample matrix has only %d columns.", ...
                string(p.days), string(p.cell), string(p.int), string(paramName), ...
                sampleCol, size(sampleMat, 2));
        end

        samples = sampleMat(:, sampleCol);
        samples = samples(:);
        samples = samples(~isnan(samples));

        lb = singleFits.(p.days).(p.cell){lbRow, drugName};
        ub = singleFits.(p.days).(p.cell){ubRow, drugName};

        singleVal = singleFits.(p.days).(p.cell){singleRow, drugName};

        combValRaw = combinationFits.(p.days).(p.cell){combRow, p.int};
        combVal = combValRaw;

        if flipCombSign
            combVal = -combVal;
        end


        lb2 = min(lb, ub);
        ub2 = max(lb, ub);


        if numel(samples) >= 10
            q = quantile(samples, [0.01 0.99]);
        else
            q = [min(samples), max(samples)];
        end

        xMin = min([q(1), lb2, ub2, singleVal, combVal]);
        xMax = max([q(2), lb2, ub2, singleVal, combVal]);

        pad = 0.10 * (xMax - xMin);

        if isnan(pad) || pad == 0
            pad = 0.1;
        end

        xMin = xMin - pad;
        xMax = xMax + pad;


        hNoise = histogram(ax, samples, ...
            "BinLimits", [xMin xMax], ...
            "NumBins", 35, ...
            "FaceColor", colorNoise, ...
            "EdgeColor", [0.25 0.25 0.45], ...
            "FaceAlpha", 0.8);

        yMax = max(hNoise.Values);

        if isempty(yMax) || isnan(yMax) || yMax == 0
            yMax = 1;
        end


        if yMax >= 1000
            yExp = floor(log10(yMax));
            yScale = 10^yExp;

            yTop = ceil(yMax / yScale) * yScale;

            ylim(ax, [0, yTop])
            ax.YTick = [0, yTop/2, yTop];

            ax.YAxis.Exponent = yExp;
            ytickformat(ax, "%.1f")
        else
            yTop = ceil(yMax / 50) * 50;

            if yTop == 0
                yTop = 1;
            end

            ylim(ax, [0, yTop])
            ax.YTick = [0, yTop/2, yTop];

            ax.YAxis.Exponent = 0;
            ytickformat(ax, "%.0f")
        end

        xlim(ax, [xMin, xMax])


        hInt = patch(ax, ...
            [lb2 ub2 ub2 lb2], ...
            [0 0 yTop yTop], ...
            colorInt, ...
            "FaceAlpha", alphaVal, ...
            "EdgeColor", "none");

        uistack(hInt, "bottom")


        xline(ax, singleVal, ...
            "LineWidth", lineWidth, ...
            "Color", colorSingle);

        xline(ax, combVal, ...
            "LineWidth", lineWidth, ...
            "Color", colorComb);


        xlabel(ax, xlab, ...
            "Interpreter", "latex", ...
            "FontSize", fs)

        ax.FontSize = fs;
        ax.LineWidth = axisLW;
        ax.TickDir = "out";
        ax.TickLength = [0.018 0.018];

        xt = linspace(xMin, xMax, 5);
        ax.XTick = xt([2 4]);
        ax.XTickLabel = compose("%.1f", ax.XTick);

        pbaspect(ax, [1 1 1])

        box(ax, "off")
        ax.XAxisLocation = "bottom";
        ax.YAxisLocation = "left";
        ax.XColor = "k";
        ax.YColor = "k";
        ax.LineWidth = axisLW;

        plot(ax, [xMin xMax], [yTop yTop], ...
            "k-", "LineWidth", axisLW, "Clipping", "on");

        plot(ax, [xMax xMax], [0 yTop], ...
            "k-", "LineWidth", axisLW, "Clipping", "on");

    end

    drawnow
    exportgraphics(fig, filePDF, "ContentType", "vector")
    close(fig)

end