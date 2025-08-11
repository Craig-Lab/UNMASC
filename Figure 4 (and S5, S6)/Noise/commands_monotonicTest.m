%This code separates the noisy data based on monotonicity.

clear variables

%Load cytokine concentrations and set names
cytConc = readtable('../../Data/cytConc.xlsx');
cytokines = cytConc.Properties.VariableNames;
cytConc = table2array(cytConc);

%Scaled concentration of cytokines
scaled_conc = [-2.366, -1, 0, 1, 2.366];

load('Results/Noise_data.mat') %noise data
load('Results/Data_values.mat')%data points

%% DAYS 7-14
%Cytokines for days 7-14 are SCF, Flt3L, IL3, IL7, TNFa, CXCL12

p.days = 'D7';

%PROT
p.cell = 'PROT';
[sol.(p.days).(p.cell),monoInc.(p.days).(p.cell),monoDec.(p.days).(p.cell),multi.(p.days).(p.cell)] = monotonicTest(p,cytokines,noise);

%CD4ISP
p.cell = 'CD4ISP';
[sol.(p.days).(p.cell),monoInc.(p.days).(p.cell),monoDec.(p.days).(p.cell),multi.(p.days).(p.cell)] = monotonicTest(p,cytokines,noise);

%CD3N
p.cell = 'CD3N';
[sol.(p.days).(p.cell),monoInc.(p.days).(p.cell),monoDec.(p.days).(p.cell),multi.(p.days).(p.cell)] = monotonicTest(p,cytokines,noise);

%% DAYS 14-28
%Cytokines for days 14-28 are SCF, Flt3L, IL3, IL7, TNFa
cytokinesD14 = cytokines(1:5); %We remove CXCL12

p.days = 'D14';

%CD3N
p.cell = 'CD3N';
[sol.(p.days).(p.cell),monoInc.(p.days).(p.cell),monoDec.(p.days).(p.cell),multi.(p.days).(p.cell)] = monotonicTest(p,cytokinesD14,noise);

%CD3P
p.cell = 'CD3P';
[sol.(p.days).(p.cell),monoInc.(p.days).(p.cell),monoDec.(p.days).(p.cell),multi.(p.days).(p.cell)] = monotonicTest(p,cytokinesD14,noise);

%SP
p.cell = 'SP';
[sol.(p.days).(p.cell),monoInc.(p.days).(p.cell),monoDec.(p.days).(p.cell),multi.(p.days).(p.cell)] = monotonicTest(p,cytokinesD14,noise);

%% Save results

save("noiseMonotonicIncreasing.mat")
save("noiseMonotonicDecreasing.mat")
save("noiseMultiphasic.mat")

%% Display figures in a tiled layout
set(0,'DefaultFigureVisible','on');

h1 = openfig("Results/Figures_D7_PROT_Bar.fig");
ax1 = gca;

h2 = openfig("Results/Figures_D7_CD4ISP_Bar.fig");
ax2 = gca;

h3 = openfig("Results/Figures_D7_CD3N_Bar.fig");
ax3 = gca;

h4 = openfig("Results/Figures_D14_CD3N_Bar.fig");
ax4 = gca;

h5 = openfig("Results/Figures_D14_CD3P_Bar.fig");
ax5 = gca;

h6 = openfig("Results/Figures_D14_SP_Bar.fig");
ax6 = gca;

x0=10;
y0=10;
width=800;
height=1000;

fig = figure();
set(fig,'position',[x0,y0,width,height])
tl = tiledlayout(1,2);

tl1 = tiledlayout(tl,3,1);
tl1.Layout.Tile = 1;
   
ax1c = copyobj(ax1, tl1);
ax1c.Layout.Tile = 1;
xticks(ax1c,[])
title(ax1c,"PRO-T cells")

ax2c = copyobj(ax2, tl1);
ax2c.Layout.Tile = 2;
xticks(ax2c,[])
title(ax2c,"CD4ISP cells")

ax3c = copyobj(ax3, tl1);
ax3c.Layout.Tile = 3;
title(ax3c,"CD3^- cells")

tl2 = tiledlayout(tl,3,1);
tl2.Layout.Tile = 2;

ax4c = copyobj(ax4, tl2);
ax4c.Layout.Tile = 1;
xticks(ax4c,[])
yticks(ax4c, [])
title(ax4c,"CD3^- cells")

ax5c = copyobj(ax5, tl2);
ax5c.Layout.Tile = 2;
xticks(ax5c,[])
yticks(ax5c, [])
title(ax5c,"CD3^+ cells")

ax6c = copyobj(ax6, tl2);
ax6c.Layout.Tile = 3;
yticks(ax6c, [])
title(ax6c,"SP cells")


leg = legend("Monotonic Increasing","Monotonic Decreasing","Multiphasic");
leg.Layout.Tile = "east";

linkaxes([ax1c ax2c ax3c ax4c ax5c ax6c],'y')
linkaxes([ax1c ax2c ax3c],'x')
linkaxes([ax4c ax5c ax6c],'x')
fontsize(16,"points")

title(tl1,"Days 7-14")
title(tl2,"Days 14-28")
title(tl,"Noise distribution")

saveas(fig,"Results/monotonicityBarGraph.pdf")

