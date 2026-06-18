%% Load data

% Cytokine data
cytokinesTable = readtable("cytokines.xlsx", "ReadRowNames", true);
cytokinesNames = cytokinesTable.Properties.RowNames;
cytokinesConc = table2array(cytokinesTable);
logConc = log10(cytokinesConc);

alpha_type = ["alpha12", "alpha21"];
nbVials = 4;
cellInt = ["SCF_IL7", "IL3_IL7"];
cellType = ["CD7pCD5n", "CD7pCD5p"];

load("surfacesParamsMC.mat")

%%  Alpha
for j = 1:numel(cellInt)
    conditionName = cellInt(j);

    for k = 1:numel(cellType)
        cellTypeName = cellType(k);

        for i = 2:nbVials
            vialName = sprintf('vial%d', i);
                
            w = surfacesParamsMC.(vialName).(conditionName).(cellTypeName);

            data_alpha12.(conditionName).(cellTypeName).(vialName) = w(:,size(w,2)-1); 
            data_alpha21.(conditionName).(cellTypeName).(vialName) = w(:,size(w,2)); 

        end 
    end 
end 

save("Alpha/alpha12.mat", "data_alpha12")
save("Alpha/alpha21.mat", "data_alpha21")


for m = 1:numel(alpha_type)
    alpha_name = alpha_type(m);

    for j = 1:numel(cellInt)
        conditionName = cellInt(j);
    
        for k= 1:numel(cellType)
            cellTypeName = cellType(k);
    
            parts = split(conditionName, "_");
            cytokine1 = parts(1);
            cytokine2 = parts(2);

            if strcmp(alpha_name,"alpha12")==1
                y_violin2 = data_alpha12.(conditionName).(cellTypeName).("vial2"); 
                y_violin3 = data_alpha12.(conditionName).(cellTypeName).("vial3"); 
                y_violin4 = data_alpha12.(conditionName).(cellTypeName).("vial4"); 
            elseif strcmp(alpha_name,"alpha12")==0
                y_violin2 = data_alpha21.(conditionName).(cellTypeName).("vial2"); 
                y_violin3 = data_alpha21.(conditionName).(cellTypeName).("vial3"); 
                y_violin4 = data_alpha21.(conditionName).(cellTypeName).("vial4"); 
            end 
    
            y_violin = [y_violin2; y_violin3; y_violin4];
            xgroupdata = [ones(numel(y_violin2),1); 2*ones(numel(y_violin3),1); 3*ones(numel(y_violin4),1)];

            xgroupdata = categorical(xgroupdata, [1 2 3], {'vial2','vial3','vial4'});
    
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
            
            v1 = violinplot(ones(numel(y_violin2),1), log10(y_violin2), ...
                FaceAlpha=alpha, DensityWidth=densityWidth);
            
            v2 = violinplot(2*ones(numel(y_violin3),1), log10(y_violin3), ...
                FaceAlpha=alpha, DensityWidth=densityWidth);
            
            v3 = violinplot(3*ones(numel(y_violin4),1), log10(y_violin4), ...
                FaceAlpha=alpha, DensityWidth=densityWidth);
            
            xlim([0.5 3.5])
            
            yl = ylim;
            ylim_val = max(abs(yl));
            
            ylim_val = 1.10 * ylim_val;
            ylim_val = ceil(ylim_val);
            
            if ylim_val == 0
                ylim_val = 1;
            end
            
            ylim([-ylim_val ylim_val])
            yticks([-ylim_val 0 ylim_val])
            
            yline(0, "LineWidth", 1, "LineStyle", "--", "Color", "#000000")
            
            if strcmp(alpha_name,"alpha12")
                ylabel(["Log-fold change of", "potency of cytokine 2"])
            else
                ylabel(["Log-fold change of", "potency of cytokine 1"])
            end

            if strcmp(cellTypeName,"CD7pCD5n") == 1 
                v1.FaceColor = color1;
                v2.FaceColor = color2;
                v3.FaceColor = color3;
            elseif strcmp(cellTypeName,"CD7pCD5n") == 0
                v1.FaceColor = color4;
                v2.FaceColor = color5;
                v3.FaceColor = color6;
            end 

            leg = legend("Donor 1", "Donor 2", "Donor 3");
            leg.Location = "bestoutside";
            
            xticks([1 2 3])
            xticklabels([])
            set(gca,'TickLength',[0 0])
            
            xlabel(sprintf('Cytokine 1: %s\nCytokine 2: %s', cytokine1, cytokine2))
            
            fontsize(fs,"points")
            H = gca;
            H.LineWidth = 1;
    
            if strcmp(alpha_name,"alpha12")==1
                file_name= strcat("Alpha/",conditionName,"_",cellTypeName,"_alpha12.pdf");
                exportgraphics(fig,file_name)
            elseif strcmp(alpha_name,"alpha12")==0
                file_name= strcat("Alpha/",conditionName,"_",cellTypeName,"_alpha21.pdf");
                exportgraphics(fig,file_name)
            end 
      
        end 
    end 
end 



%% Beta 

for j = 1:numel(cellInt)
    conditionName = cellInt(j);

    for k = 1:numel(cellType)
        cellTypeName = cellType(k);

        for i = 2:nbVials
            vialName = sprintf('vial%d', i);
                
            w = surfacesParamsMC.(vialName).(conditionName).(cellTypeName);

            beta_vals = zeros(size(w,1),1);

            for m = 1:size(w,1)
                beta_vals(m) = beta_fun(w(m,:));
            end

            data_beta.(conditionName).(cellTypeName).(vialName) = beta_vals;
        end 
    end 
end 

save("Beta/beta.mat", "data_beta")


for j = 1:numel(cellInt)
    conditionName = cellInt(j);

    for k= 1:numel(cellType)
        cellTypeName = cellType(k);

        parts = split(conditionName, "_");
        cytokine1 = parts(1);
        cytokine2 = parts(2);

        raw2 = data_beta.(conditionName).(cellTypeName).("vial2");
        y_violin2 = raw2(raw2 >= -5 & raw2 <= 5);
        
        raw3 = data_beta.(conditionName).(cellTypeName).("vial3");
        y_violin3 = raw3(raw3 >= -5 & raw3 <= 5);
        
        raw4 = data_beta.(conditionName).(cellTypeName).("vial4");
        y_violin4 = raw4(raw4 >= -5 & raw4 <= 5);
   

        y_violin = [y_violin2; y_violin3; y_violin4];

        xgroupdata = [ones(numel(y_violin2),1); 2*ones(numel(y_violin3),1); 3*ones(numel(y_violin4),1)];

        xgroupdata = categorical(xgroupdata, [1 2 3], {'vial2','vial3','vial4'});

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

        yy = y_violin;
        yy = yy(isfinite(yy));
        
        y_violin2 = y_violin2(isfinite(y_violin2));
        y_violin3 = y_violin3(isfinite(y_violin3));
        y_violin4 = y_violin4(isfinite(y_violin4));
            
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
        
        ylim_val = 1.10 * ylim_val;
        ylim_val = ceil(ylim_val);
        
        if ylim_val == 0
            ylim_val = 1;
        end
        
        ylim([-ylim_val ylim_val])
        yticks([-ylim_val 0 ylim_val])
        
        yline(0, "LineWidth", 1, "LineStyle", "--", "Color", "#000000")

       ylabel("\beta",'Interpreter','tex','FontWeight','bold')
    
        if strcmp(cellTypeName,"CD7pCD5n") == 1 
            v1.FaceColor = color1;
            v2.FaceColor = color2;
            v3.FaceColor = color3;
        elseif strcmp(cellTypeName,"CD7pCD5n") == 0
            v1.FaceColor = color4;
            v2.FaceColor = color5;
            v3.FaceColor = color6;
        end 

        leg = legend("Vial 2", "Vial 3", "Vial 4");
        leg.Location = "bestoutside";
        
        xticks([1 2 3])
        xticklabels([])
        set(gca,'TickLength',[0 0])
        
        xlabel(sprintf('Cytokine 1: %s\nCytokine 2: %s', cytokine1, cytokine2))
        
        fontsize(fs,"points")
        H = gca;
        H.LineWidth = 1;
      
        file_name= strcat("Beta/",conditionName,"_",cellTypeName,".pdf");
        exportgraphics(fig,file_name)
    end 
end 

function beta = beta_fun(param)

    E0 = param(1);
    E1 = param(2);
    E2 = param(3);
    E3 = param(4);

    beta = (E0-max(E1,E2))./(max(E1,E2)-E3);
end 
 
