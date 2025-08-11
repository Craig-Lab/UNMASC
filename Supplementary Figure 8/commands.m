clear variables

%Load cytokine concentrations and set names
cytConc = readtable('../Data/cytConc.xlsx');
cytokines = cytConc.Properties.VariableNames;
cytConc = table2array(cytConc);

%Load color for surfaces
p.map = load("ownColor.mat");

%Load surfaces results
load("fit_parameters.mat")
load("parameters_MC.mat")
load("CI.mat")


%% TNFA_CXCL12 DAYS 7-14 PRO-T CELLS

p.days = "D7";
p.cells = "PROT";
p.int = "TNFa_CXCL12";
p.drug1_name = "TNFa";
p.drug2_name = "CXCL12";

p.monotonic = 2;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% TNFA_CXCL12 DAYS 7-14 CD4ISP CELLS

p.days = "D7";
p.cells = "CD4ISP";
p.int = "TNFa_CXCL12";
p.drug1_name = "TNFa";
p.drug2_name = "CXCL12";

p.monotonic = 2;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% TNFA_CXCL12 DAYS 7-14 CD3N CELLS

p.days = "D7";
p.cells = "CD3N";
p.int = "TNFa_CXCL12";
p.drug1_name = "TNFa";
p.drug2_name = "CXCL12";

p.monotonic = 3;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% IL3-IL7 DAYS 14-28 CD3N CELLS

p.days = "D14";
p.cells = "CD3N";
p.int = "IL3_IL7";
p.drug1_name = "IL3";
p.drug2_name = "IL7";

p.monotonic = 1;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% IL7-TNFa DAYS 14-28 CD3N CELLS 

p.days = "D14";
p.cells = "CD3N";
p.int = "IL7_TNFa";
p.drug1_name = "IL7";
p.drug2_name = "TNFa";

p.monotonic = 2;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% IL3-IL7 DAYS 14-28 CD3P CELLS

p.days = "D14";
p.cells = "CD3P";
p.int = "IL3_IL7";
p.drug1_name = "IL3";
p.drug2_name = "IL7";

p.monotonic = 1;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% IL7-TNFa DAYS 14-28 CD3N CELLS 

p.days = "D14";
p.cells = "CD3P";
p.int = "IL7_TNFa";
p.drug1_name = "IL7";
p.drug2_name = "TNFa";

p.monotonic = 2;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% IL3-IL7 DAYS 14-28 SP CELLS

p.days = "D14";
p.cells = "SP";
p.int = "IL3_IL7";
p.drug1_name = "IL3";
p.drug2_name = "IL7";

p.monotonic = 1;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% IL7-TNFa DAYS 14-28 CD3N CELLS 

p.days = "D14";
p.cells = "SP";
p.int = "IL7_TNFa";
p.drug1_name = "IL7";
p.drug2_name = "TNFa";

p.monotonic = 2;

param = fit_parameters.(p.days).(p.cells).(p.int);

fig = new_figure_fun(p,param);

%% ---------------------------------------

function fig = new_figure_fun(p,param)
        
        alphaVal = 0.7;

        drug1_new = logspace(log10(0.001),log10(1000));
        drug2_new = logspace(log10(0.001),log10(1000));

        E = surface_fun(drug1_new,drug2_new,param,p);

        [X,Y] = meshgrid(log10(drug1_new), log10(drug2_new));   
        
        fig = figure('Visible','off','OuterPosition',[100 100 500 500]);
        x0=10;
        y0=10;
        width=400;
        height=300;
        set(fig,'position',[x0,y0,width,height])

        surf(X,Y,E,'EdgeColor','none','FaceAlpha',alphaVal);
        grid off
        colormap(fig,p.map.map)
        clim([0 1])
        cb = colorbar;   
        cb.Limits = [0 1];
        cb.Ticks = [0 0.2 0.4 0.6 0.8 1];
        cb.Label.String = "Effect";
        cb.Label.VerticalAlignment = "middle";
        cb.Label.Rotation = 0;
        cb.Label.Position = [0.5 1+0.1];
        set(cb,'position',[0.9 0.2 .02 .5])
       
        zlim([0 1.2])
        zticks([0 0.5 1])
        xlh = xlabel(strcat("Log_{10}(",p.drug1_name,")"),"Rotation",12.5);
        xlh.Position(1) = xlh.Position(1); 
        xlh.Position(2) = xlh.Position(2) + 0.35;
        ylh = ylabel(strcat("Log_{10}(",p.drug2_name,")"),"Rotation",-27.5);
        ylh.Position(2) = ylh.Position(2) - 0.40;
        ylh.Position(1) = ylh.Position(1) + 0.25;
        zlabel("Effect")
        fontsize(11,"points")
        axis square
        view(-33.56,17.79)
        
        %Set Line width
        H=gca;
        H.LineWidth=1.5;

        %Set alpha for whole figure
        alpha(.7)

        file_name = strcat("Extended/",p.days,"_",p.cells,"_",p.drug1_name,"_",p.drug2_name,".pdf");

        exportgraphics(fig,file_name)
end 
