%% Load data

% Cytokine data
cytokinesTable = readtable("cytConc.xlsx", "ReadVariableNames", true);
cytokinesNames = cytokinesTable.Properties.VariableNames;
cytokinesConc = table2array(cytokinesTable);

load("surfacesParamsMC.mat")

daysName = ["D7", "D14"];

cellType1 = ["PROT", "CD4ISP", "CD3N"];
cellType2 = ["CD3N", "CD3P", "SP"];

cellIntNames1 = ["SCF_IL7", "IL3_IL7", "IL7_TNFa", "TNFa_CXCL12"];
cellIntNames2 = ["SCF_IL7", "IL3_IL7", "IL7_TNFa"];

alpha_type = ["alpha12", "alpha21"];

%% Alpha

for i = 1:numel(daysName)

    days = daysName(i);

    if strcmp(days, "D7") == 1
        cellTypeNames = cellType1;
        cellIntNames = cellIntNames1;
        cytNames = cytokinesNames;
    elseif strcmp(days, "D14") == 1
        cellTypeNames = cellType2;
        cytNames = cytokinesNames(1:5);
        cellIntNames = cellIntNames2; 
    end 

    for j = 1:numel(cellTypeNames)

        cellType = cellTypeNames(j);

        for k = 1:numel(cellIntNames)

            cellInt = cellIntNames(k);
                
            w = surfacesParamsMC.(days).(cellType).(cellInt);

            data_alpha12.(days).(cellInt).(cellType) = w(:, size(w, 2)-1); 
            data_alpha21.(days).(cellInt).(cellType) = w(:, size(w, 2)); 

        end 

    end 

end 

save("Alpha/alpha12_iPSC.mat", "data_alpha12")
save("Alpha/alpha21_iPSC.mat", "data_alpha21")

for i = 1:numel(daysName)

    days = daysName(i);

    if strcmp(days, "D7") == 1
        cellTypeNames = cellType1;
        cellIntNames = cellIntNames1;
        cytNames = cytokinesNames;
    elseif strcmp(days, "D14") == 1
        cellTypeNames = cellType2;
        cytNames = cytokinesNames(1:5);
        cellIntNames = cellIntNames2; 
    end 

    for m = 1:numel(alpha_type)

        alpha_name = alpha_type(m);
    
        for j = 1:numel(cellIntNames)

            conditionName = cellIntNames(j);
    
            parts = split(conditionName, "_");
            cytokine1 = parts(1);
            cytokine2 = parts(2);

            if strcmp(alpha_name, "alpha12") == 1

                y_violin2 = data_alpha12.(days).(conditionName).(cellTypeNames(1)); 
                y_violin3 = data_alpha12.(days).(conditionName).(cellTypeNames(2)); 
                y_violin4 = data_alpha12.(days).(conditionName).(cellTypeNames(3)); 

            elseif strcmp(alpha_name, "alpha21") == 1

                y_violin2 = data_alpha21.(days).(conditionName).(cellTypeNames(1)); 
                y_violin3 = data_alpha21.(days).(conditionName).(cellTypeNames(2)); 
                y_violin4 = data_alpha21.(days).(conditionName).(cellTypeNames(3)); 

            end 

            y_violin = [y_violin2; y_violin3; y_violin4];

            xgroupdata = [ones(numel(y_violin2), 1); ...
                          2*ones(numel(y_violin3), 1); ...
                          3*ones(numel(y_violin4), 1)];

            xgroupdata = categorical(xgroupdata, [1 2 3], {'vial2', 'vial3', 'vial4'});

            x0 = 0;
            y0 = 0;
            width = 7;
            height = 5;
            alpha = 1;
            fs = 9;
            densityWidth = 0.5;

            color1 = "#8B5058";
            color2 = "#525B7A";
            color3 = "#EC9F79";
            color4 = "#A781A2";
            color5 = "#9FC2CC";
            color6 = "#d17466";

            fig = figure("Visible", "off");
            set(fig, 'units', 'centimeters', 'position', [x0, y0, width, height])
            hold on

            v1 = violinplot(ones(numel(y_violin2), 1), log10(y_violin2), ...
                FaceAlpha=alpha, DensityWidth=densityWidth);

            v2 = violinplot(2*ones(numel(y_violin3), 1), log10(y_violin3), ...
                FaceAlpha=alpha, DensityWidth=densityWidth);

            v3 = violinplot(3*ones(numel(y_violin4), 1), log10(y_violin4), ...
                FaceAlpha=alpha, DensityWidth=densityWidth);

            xlim([0.5 3.5])

            yl = ylim;
            ylim_val = max(abs(yl));

            ylim_val = ceil(ylim_val);

            if ylim_val == 0
                ylim_val = 1;
            end

            ylim([-ylim_val ylim_val])
            yticks([-ylim_val 0 ylim_val])

            yline(0, "LineWidth", 1, "LineStyle", "--", "Color", "#000000")

            if strcmp(alpha_name, "alpha12")
                ylabel(["Log-fold change of", "potency of cytokine 2"])
            else
                ylabel(["Log-fold change of", "potency of cytokine 1"])
            end

            if strcmp(days, "D7") == 1 

                v1.FaceColor = color1;
                v2.FaceColor = color2;
                v3.FaceColor = color3;

            elseif strcmp(days, "D7") == 0

                v1.FaceColor = color3;
                v2.FaceColor = color4;
                v3.FaceColor = color5;

            end 

            leg = legend(cellTypeNames(1), cellTypeNames(2), cellTypeNames(3));
            leg.Location = "bestoutside";

            set(gca, 'TickLength', [0 0])

            xlabel(sprintf('Cytokine 1: %s\nCytokine 2: %s', cytokine1, cytokine2))

            fontsize(fs, "points")

            H = gca;
            H.LineWidth = 1;

            if strcmp(alpha_name, "alpha12") == 1

                file_name = strcat("Alpha/", days, "_", conditionName, "_alpha12.pdf");
                exportgraphics(fig, file_name)

            elseif strcmp(alpha_name, "alpha12") == 0

                file_name = strcat("Alpha/", days, "_", conditionName, "_alpha21.pdf");
                exportgraphics(fig, file_name)

            end 

        end 

    end 
end 

%% Code for violin plots of beta MC parameters

for i = 1:numel(daysName)
    days = daysName(i);

    if strcmp(days,"D7")
        cellTypeNames = cellType1;
        cellIntNames = cellIntNames1;
    else
        cellTypeNames = cellType2;
        cellIntNames = cellIntNames2;
    end

    for j = 1:numel(cellTypeNames)
        cellType = cellTypeNames(j);

        for k = 1:numel(cellIntNames)
            cellInt = cellIntNames(k);

            w = surfacesParamsMC.(days).(cellType).(cellInt)(:,1:4);

            beta_vals = zeros(size(w,1),1);

            for m = 1:size(w,1)
                beta_vals(m) = beta_fun(w(m,:));
            end
            
            data_beta.(days).(cellInt).(cellType) = beta_vals;

        end
    end
end

save("Beta/beta_iPSC.mat", "data_beta")

for i = 1:numel(daysName)
    days = daysName(i);

    if strcmp(days,"D7")
        cellTypeNames = cellType1;
        cellIntNames = cellIntNames1;
    else
        cellTypeNames = cellType2;
        cellIntNames = cellIntNames2;
    end

    for j = 1:numel(cellIntNames)
        conditionName = cellIntNames(j);

        parts = split(conditionName, "_");
        cytokine1 = parts(1);
        cytokine2 = parts(2);

        y_violin2 = data_beta.(days).(conditionName).(cellTypeNames(1));
        y_violin3 = data_beta.(days).(conditionName).(cellTypeNames(2));
        y_violin4 = data_beta.(days).(conditionName).(cellTypeNames(3));

        y_violin = [y_violin2; y_violin3; y_violin4];

        x0=0;
        y0=0;
        width=7;
        height=5;
        alpha = 1;
        fs = 9; 
        densityWidth = 0.3;

        color1 = "#8B5058";
        color2 = "#525B7A";
        color3 = "#EC9F79";
        color4 = "#A781A2";
        color5 = "#9FC2CC";
        color6 = "#d17466";

        fig = figure("Visible","off");
            set(fig,'units','centimeters','position',[x0,y0,width,height])
            hold on

        v1 = violinplot(ones(numel(y_violin2),1), y_violin2, ...
            FaceAlpha=alpha, DensityWidth=densityWidth);
        v2 = violinplot(2*ones(numel(y_violin3),1), y_violin3, ...
            FaceAlpha=alpha, DensityWidth=densityWidth);
        v3 = violinplot(3*ones(numel(y_violin4),1), y_violin4, ...
            FaceAlpha=alpha, DensityWidth=densityWidth);

        xlim([0.5 3.5])
        
        yl = ylim;
        ylim_val = max(abs(yl));
        
        ylim_val = ceil(ylim_val);
        
        if ylim_val == 0
            ylim_val = 1;
        end
        
        ylim([-ylim_val ylim_val])
        yticks([-ylim_val 0 ylim_val])
        
        yline(0, "LineWidth", 1, "LineStyle", "--", "Color", "#000000")

        ylabel('$\beta$','Interpreter','latex')

        if strcmp(days,"D7") == 1 
            v1.FaceColor = color1;
            v2.FaceColor = color2;
            v3.FaceColor = color3;
        elseif strcmp(days,"D7") == 0
            v1.FaceColor = color3;
            v2.FaceColor = color4;
            v3.FaceColor = color5;
        end 

        leg = legend(cellTypeNames(1), cellTypeNames(2), cellTypeNames(3));
        leg.Location = "bestoutside";

        xlabel(sprintf('Cytokine 1: %s\nCytokine 2: %s', cytokine1, cytokine2))

        set(gca,'TickLength',[0 0])
        H = gca;
        H.LineWidth = 1;

        fontsize(9,"points")

        file_name = strcat("Beta/",days,"_",conditionName,"_beta.pdf");
        exportgraphics(fig,file_name)
    end
end


%% -----------------------------------------------------------------------

function beta = beta_fun(param)

    E0 = param(1);
    E1 = param(2);
    E2 = param(3);
    E3 = param(4);

    beta = (E0-max(E1,E2))./(max(E1,E2)-E3);
end 