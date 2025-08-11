%This code generates the figures displaying the single fits and the noise
%data generated.

%% Make figure Days 7-14


set(0,'DefaultFigureVisible','on');

%PRO-T CELLS
openfig("Results/D7/PROT/SCF_Fit.fig");
ax1 = gca;
xticks(ax1,[])
title(ax1,"")

openfig("Results/D7/PROT/Flt3L_Fit.fig");
ax2 = gca;
yticks(ax2,[])
xticks(ax2,[])
title(ax2,"")

openfig("Results/D7/PROT/IL3_Fit.fig");
ax3 = gca;
yticks(ax3,[])
xticks(ax3,[])
title(ax3,"")

openfig("Results/D7/PROT/IL7_Fit.fig");
ax4 = gca;
yticks(ax4,[])
xticks(ax4,[])
title(ax4,"")

openfig("Results/D7/PROT/TNFa_Fit.fig");
ax5 = gca;
yticks(ax5,[])
xticks(ax5,[])
title(ax5,"")

openfig("Results/D7/PROT/CXCL12_Fit.fig");
ax6 = gca;
yticks(ax6,[])
xticks(ax6,[])
title(ax6,"")

%CD4ISP cells
openfig("Results/D7/CD4ISP/SCF_Fit.fig");
ax7 = gca;
title(ax7,"")
xticks(ax7,[])

openfig("Results/D7/CD4ISP/Flt3L_Fit.fig");
ax8 = gca;
title(ax8,"")
yticks(ax8,[])
xticks(ax8,[])

openfig("Results/D7/CD4ISP/IL3_Fit.fig");
ax9 = gca;
title(ax9,"")
yticks(ax9,[])
xticks(ax9,[])

openfig("Results/D7/CD4ISP/IL7_Fit.fig");
ax10 = gca;
title(ax10,"")
yticks(ax10,[])
xticks(ax10,[])

openfig("Results/D7/CD4ISP/TNFa_Fit.fig");
ax11 = gca;
title(ax11,"")
yticks(ax11,[])
xticks(ax11,[])

openfig("Results/D7/CD4ISP/CXCL12_Fit.fig");
ax12 = gca;
title(ax12,"")
yticks(ax12,[])
xticks(ax12,[])

%CD3N cells
openfig("Results/D7/CD3N/SCF_Fit.fig");
ax13 = gca;
title(ax13,"")
xlabel(ax13,"SCF")

openfig("Results/D7/CD3N/Flt3L_Fit.fig");
ax14 = gca;
title(ax14,"")
yticks(ax14,[])
xlabel(ax14,"Flt3L")

openfig("Results/D7/CD3N/IL3_Fit.fig");
ax15 = gca;
title(ax15,"")
yticks(ax15,[])
xlabel(ax15,"IL-3")

openfig("Results/D7/CD3N/IL7_Fit.fig");
ax16 = gca;
title(ax16,"")
yticks(ax16,[])
xlabel(ax16,"IL-7")

openfig("Results/D7/CD3N/TNFa_Fit.fig");
ax17 = gca;
title(ax17,"")
yticks(ax17,[])
xlabel(ax17,"TNFa")

openfig("Results/D7/CD3N/CXCL12_Fit.fig");
ax18 = gca;
title(ax18,"")
yticks(ax18,[])
xlabel(ax18,"CXCL12")

figure("Units","centimeters","Position",[10 10 20 10]);

T = tiledlayout(3,1);
T.TileSpacing = "tight";
T.Padding = "tight";

tl1 = tiledlayout(T,1,6);
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [1 1];
title(tl1,"ProT","fontweight","bold")

tl2 = tiledlayout(T,1,6);
tl2.Layout.Tile = 2;
tl2.Layout.TileSpan = [1 1];
title(tl2,"CD4ISP (CD3-)","fontweight","bold")

tl3 = tiledlayout(T,1,6);
tl3.Layout.Tile = 3;
tl3.Layout.TileSpan = [1 1];
title(tl3,"DP (CD3-)","fontweight","bold")

ax1c = copyobj(ax1, tl1);
ax1c.Layout.Tile = 1;

ax2c = copyobj(ax2, tl1);
ax2c.Layout.Tile = 2;

ax3c = copyobj(ax3, tl1);
ax3c.Layout.Tile = 3;

ax4c = copyobj(ax4, tl1);
ax4c.Layout.Tile = 4;

ax5c = copyobj(ax5, tl1);
ax5c.Layout.Tile = 5;

ax6c = copyobj(ax6, tl1);
test1 = ax6c.Children(4);
test2 = ax6c.Children(1);
test3 = ax6c.Children(3);
ax6c.Layout.Tile = 6;

leg = legend([test1 test2 test3],{"Noise","Data + Error", "Fit"});
leg.Location = "bestoutside";

ax7c = copyobj(ax7, tl2);
ax7c.Layout.Tile = 1;

ax8c = copyobj(ax8, tl2);
ax8c.Layout.Tile = 2;

ax9c = copyobj(ax9, tl2);
ax9c.Layout.Tile = 3;

ax10c = copyobj(ax10, tl2);
ax10c.Layout.Tile = 4;

ax11c = copyobj(ax11, tl2);
ax11c.Layout.Tile = 5;

ax12c = copyobj(ax12, tl2);
ax12c.Layout.Tile = 6;

ax13c = copyobj(ax13, tl3);
ax13c.Layout.Tile = 1;

ax14c = copyobj(ax14, tl3);
ax14c.Layout.Tile = 2;

ax15c = copyobj(ax15, tl3);
ax15c.Layout.Tile = 3;

ax16c = copyobj(ax16, tl3);
ax16c.Layout.Tile = 4;

ax17c = copyobj(ax17, tl3);
ax17c.Layout.Tile = 5;

ax18c = copyobj(ax18, tl3);
ax18c.Layout.Tile = 6;

fontsize(10,"points")
ylabel(T,"Effect","fontweight","bold")
xlabel(T,"Log_{10}(Concentration (ng/mL))",'fontweight','bold')

%% Make figure Days 14-28

set(0,'DefaultFigureVisible','on');

%CD3N CELLS
h1 = openfig("Results/D14/CD3N/SCF_Fit.fig");
ax1 = gca;
xticks(ax1,[])
yticks(ax1,[0 0.5 1 1.5])
ylim(ax1,[0 2])
title(ax1,"")

h2 = openfig("Results/D14/CD3N/Flt3L_Fit.fig");
ax2 = gca;
yticks(ax2,[])
xticks(ax2,[])
ylim(ax2,[0 2])
title(ax2,"")

h3 = openfig("Results/D14/CD3N/IL3_Fit.fig");
ax3 = gca;
yticks(ax3,[])
xticks(ax3,[])
ylim(ax3,[0 2])
title(ax3,"")

h4 = openfig("Results/D14/CD3N/IL7_Fit.fig");
ax4 = gca;
yticks(ax4,[])
xticks(ax4,[])
ylim(ax4,[0 2])
title(ax4,"")

h5 = openfig("Results/D14/CD3N/TNFa_Fit.fig");
ax5 = gca;
yticks(ax5,[])
xticks(ax5,[])
ylim(ax5,[0 2])
title(ax5,"")

%CD4ISP cells
h7 = openfig("Results/D14/CD3P/SCF_Fit.fig");
ax7 = gca;
title(ax7,"")
ylim(ax7,[0 2])
yticks(ax7,[0 0.5 1 1.5])
xticks(ax7,[])

h8 = openfig("Results/D14/CD3P/Flt3L_Fit.fig");
ax8 = gca;
title(ax8,"")
yticks(ax8,[])
xticks(ax8,[])
ylim(ax8,[0 2])

h9 = openfig("Results/D14/CD3P/IL3_Fit.fig");
ax9 = gca;
title(ax9,"")
yticks(ax9,[])
xticks(ax9,[])
ylim(ax9,[0 2])

h10 = openfig("Results/D14/CD3P/IL7_Fit.fig");
ax10 = gca;
title(ax10,"")
yticks(ax10,[])
xticks(ax10,[])
ylim(ax10,[0 2])

h11 = openfig("Results/D14/CD3P/TNFa_Fit.fig");
ax11 = gca;
title(ax11,"")
yticks(ax11,[])
xticks(ax11,[])
ylim(ax11,[0 2])

%CD3N cells
h13 = openfig("Results/D14/SP/SCF_Fit.fig");
ax13 = gca;
title(ax13,"")
xlabel(ax13,"SCF")
ylim(ax13,[0 2])
yticks(ax13,[0 0.5 1 1.5])

h14 = openfig("Results/D14/SP/Flt3L_Fit.fig");
ax14 = gca;
title(ax14,"")
yticks(ax14,[])
xlabel(ax14,"Flt3L")
ylim(ax14,[0 2])

h15 = openfig("Results/D14/SP/IL3_Fit.fig");
ax15 = gca;
title(ax15,"")
yticks(ax15,[])
xlabel(ax15,"IL-3")
ylim(ax15,[0 2])

h16 = openfig("Results/D14/SP/IL7_Fit.fig");
ax16 = gca;
title(ax16,"")
yticks(ax16,[])
xlabel(ax16,"IL-7")
ylim(ax16,[0 2])

h17 = openfig("Results/D14/SP/TNFa_Fit.fig");
ax17 = gca;
title(ax17,"")
yticks(ax17,[])
xlabel(ax17,"TNFa")
ylim(ax17,[0 2])

figure("Units","centimeters","Position",[10 10 20 10]);
T = tiledlayout(3,1);
T.TileSpacing = "tight";
T.Padding = "tight";

tl1 = tiledlayout(T,1,5);
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [1 1];
title(tl1,"DP (CD3-)","fontweight","bold")

tl2 = tiledlayout(T,1,5);
tl2.Layout.Tile = 2;
tl2.Layout.TileSpan = [1 1];
title(tl2,"DP (CD3+)","fontweight","bold")

tl3 = tiledlayout(T,1,5);
tl3.Layout.Tile = 3;
tl3.Layout.TileSpan = [1 1];
title(tl3,"8SP (CD3+)","fontweight","bold")

ax1c = copyobj(ax1, tl1);
ax1c.Layout.Tile = 1;

ax2c = copyobj(ax2, tl1);
ax2c.Layout.Tile = 2;

ax3c = copyobj(ax3, tl1);
ax3c.Layout.Tile = 3;

ax4c = copyobj(ax4, tl1);
ax4c.Layout.Tile = 4;

ax5c = copyobj(ax5, tl1);

test1 = ax5c.Children(4);
test2 = ax5c.Children(1);
test3 = ax5c.Children(3);
ax5c.Layout.Tile = 5;

leg = legend([test1 test2 test3],{"Noise","Data + Error", "Fit"});
leg.Location = "bestoutside";


ax7c = copyobj(ax7, tl2);
ax7c.Layout.Tile = 1;

ax8c = copyobj(ax8, tl2);
ax8c.Layout.Tile = 2;

ax9c = copyobj(ax9, tl2);
ax9c.Layout.Tile = 3;

ax10c = copyobj(ax10, tl2);
ax10c.Layout.Tile = 4;

ax11c = copyobj(ax11, tl2);
ax11c.Layout.Tile = 5;


ax13c = copyobj(ax13, tl3);
ax13c.Layout.Tile = 1;

ax14c = copyobj(ax14, tl3);
ax14c.Layout.Tile = 2;

ax15c = copyobj(ax15, tl3);
ax15c.Layout.Tile = 3;

ax16c = copyobj(ax16, tl3);
ax16c.Layout.Tile = 4;

ax17c = copyobj(ax17, tl3);
ax17c.Layout.Tile = 5;


fontsize(10,"points")
ylabel(T,"Effect","fontweight","bold")
xlabel(T,"Log_{10}(Concentration (ng/mL))",'fontweight','bold')
