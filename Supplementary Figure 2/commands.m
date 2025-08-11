%This code displays how different values of the multiphasic curve impact
%the shape of the resulting effect. 

%% Test multiphasic EC501 effect

E01 = 1;
E02 = 1;
E1 = 0;
E2 = 0;
h1 = -2;
h2 = 2;
C1 = [0.25 0.5 0.75 1];
C2 = 2;

drug = log10(logspace(log10(0.1),log10(400),500));

param1 = [E01, h1, C1(1), E1, E02, h2, C2, E2];
param2 = [E01, h1, C1(2), E1, E02, h2, C2, E2];
param3 = [E01, h1, C1(3), E1, E02, h2, C2, E2];
param4 = [E01, h1, C1(4), E1, E02, h2, C2, E2];

[s11,s12,curve1] = multi(drug,param1);
[s21,s22,curve2]= multi(drug, param2);
[s31,s32,curve3] = multi(drug, param3);
[s41,s42,curve4] = multi(drug,param4);

fig = figure('Visible','on');
x0=10;
y0=10;
width=800;
height=200;
set(fig,'position',[x0,y0,width,height])

tl = tiledlayout(1,4,"TileSpacing","tight");
lw1 = 1.5; %line width
lw2 = 3;
colorCurve = "#7a5186";
colors1 = "#e68a00";
colors2 = "#588477";
colors3 = "#ad862c";

nexttile;
plot(drug,curve1,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s11,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s12,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("C_1 = ", mat2str(C1(1))))
ylim([0 1.1])
% ylabel("Effect")
axis square
H=gca;
H.LineWidth=1.5;

nexttile;
plot(drug,curve2,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s21,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s22,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("C_1 = ", mat2str(C1(2))))
ylim([0 1.1])
yticks([])
axis square
H=gca;
H.LineWidth=1.5;

nexttile;
plot(drug,curve3,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s31,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s32,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("C_1 = ", mat2str(C1(3))))
ylim([0 1.1])
yticks([])
axis square
H=gca;
H.LineWidth=1.5;

nexttile;
plot(drug,curve4,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s41,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s42,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("C_1 = ", mat2str(C1(4))))
ylim([0 1.1])
yticks([])
axis square
H=gca;
H.LineWidth=1.5;

leg = legend('Multiphasic','Phase 1','Phase 2');
leg.Location = "northeastoutside";
fontsize(14,"points")
% tl.XLabel.String = "Log(drug)";

file_name = "C1_effect.pdf";
exportgraphics(fig,file_name)

%% Test multiphasic h1 effect

E01 = 1;
E02 = 1;
E1 = 0;
E2 = 0;
h1 = [-2, -1.5, -1, -0.5];
h2 = 2;
C1 = 0.25;
C2 = 2;

drug = log10(logspace(log10(0.1),log10(400),500));

param1 = [E01, h1(1), C1, E1, E02, h2, C2, E2];
param2 = [E01, h1(2), C1, E1, E02, h2, C2, E2];
param3 = [E01, h1(3), C1, E1, E02, h2, C2, E2];
param4 = [E01, h1(4), C1, E1, E02, h2, C2, E2];

[s11,s12,curve1] = multi(drug,param1);
[s21,s22,curve2]= multi(drug, param2);
[s31,s32,curve3] = multi(drug, param3);
[s41,s42,curve4] = multi(drug,param4);

fig = figure('Visible','on');
x0=10;
y0=10;
width=800;
height=200;
set(fig,'position',[x0,y0,width,height])

tl = tiledlayout(1,4,'TileSpacing','tight');
lw1 = 1.5; %line width
lw2 = 3;
colorCurve = "#7a5186";
colors1 = "#e68a00";
colors2 = "#588477";

nexttile;
plot(drug,curve1,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s11,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s12,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("h_1 = ", mat2str(h1(1))))
ylim([0 1.1])
axis square
H=gca;
H.LineWidth=1.5;

nexttile;
plot(drug,curve2,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s21,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s22,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("h_1 = ", mat2str(h1(2))))
ylim([0 1.1])
yticks([])
axis square
H=gca;
H.LineWidth=1.5;

nexttile;
plot(drug,curve3,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s31,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s32,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("h_1 = ", mat2str(h1(3))))
ylim([0 1.1])
yticks([])
axis square
H=gca;
H.LineWidth=1.5;

nexttile;
plot(drug,curve4,"LineWidth",lw1,"Color",colorCurve)
hold on
plot(drug,s41,'LineWidth',lw2,"Color",colors3,'LineStyle','--')
hold on
plot(drug,s42,"LineWidth",lw2,"Color",colors2,"LineStyle","--")
hold off
title(strcat("h_1 = ", mat2str(h1(4))))
ylim([0 1.1])
yticks([])
axis square
H=gca;
H.LineWidth=1.5;

leg = legend('Multiphasic','Phase 1','Phase 2');
leg.Location = "northeastoutside";
fontsize(14,"points")

file_name = "h1_effect.pdf";
exportgraphics(fig,file_name)



%% Functions
function [s1,s2,d] = multi(drug,param)
    E01 = param(1);
    h1 = param(2);
    C1 = param(3);
    E1 = param(4);

    E02 = param(5);
    h2 = param(6);
    C2 = param(7);
    E2 = param(8);

    single1 = @(d) E01 + (E1-E01)./(1+10.^((C1-d)*h1));
    single2 = @(d) E02 + (E2-E02)./(1+10.^((C2-d)*h2));
    double = @(d) single1(d).*single2(d);
    
    s1 = single1(drug);
    s2 = single2(drug);
    d = double(drug);
end 

function s = single(drug,param)
    E0 = param(1);
    h = param(2);
    C = param(3);
    E1 = param(4);
    
    single = @(d) E0 + (E1-E0)./(1+10.^((C-d)*h));

    s = single(drug);
end 



