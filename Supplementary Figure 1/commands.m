%% Hill vs multiphasic

drug = log10(logspace(log10(0.1),log10(1000),500));

E0 = 0.2;
h = 2;
C = 0.5;
E1 = 1;

param_hill = [E0, h, C, E1];

E01 = 0.2;
h1 = 2;
C1 = 0.1;
Emax1 = 1;

E02 = 0.4;
h2 = -3;
C2 = 2;
Emax2 = 1;

param_phase1 = [E01, h1, C1, Emax1];
param_phase2 = [E02, h2, C2, Emax2];

param_multi = [E01, h1, C1, Emax1, E02, h2, C2, Emax2];

s = single(drug,param_hill);
[phase1,phase2,m] = multi(drug,param_multi);

colorCurve = "#7a5186";
colors1 = "#847852";
colors2 = "#588477";
colors3 = "#2c53ad";

fig = figure('Visible','on');
x0=10;
y0=10;
width=1000;
height=500;
set(fig,'position',[x0,y0,width,height])

tl = tiledlayout(1,2,"TileSpacing","loose");

nexttile;
line([log10(0.1) log10(1000)], [E1 E1],"LineStyle","--","LineWidth",3,"Color",colors1)
hold on
plot(drug,s,'LineWidth',3,"Color",colorCurve)
hold on
plot(C,0.6,'o',"MarkerFaceColor",colors3,"MarkerEdgeColor",colors3,"MarkerSize",8)
hold on
line([log10(0.1) C], [0.6 0.6],"LineWidth",3,"Color",colors3,"LineStyle","--")
hold on
line([C C], [0 0.6],"LineWidth",3,"Color",colors3,"LineStyle","--")
hold off
ylim([0 1.1])
xlim([log10(0.1) log10(1000)])
yticks([E0 E1])
yticklabels(["E_0" "E_1"])
xticks(C)
xticklabels("C")
axis square
xlabel("Log(agent)")
%Set Line width
H=gca;
H.LineWidth=1.5;

nexttile;
line([log10(0.1) log10(1000)], [E02 E02],"LineStyle","--","LineWidth",3,"Color",colors1)
hold on
line([log10(0.1) 1.15], [Emax1 Emax1],"LineStyle","--","LineWidth",3,"Color",colors1)
hold on
plot(drug,m,"LineWidth",3,"LineStyle","-","Color",colorCurve)
hold on
plot(C1,0.6,'o',"MarkerFaceColor",colors3,"MarkerEdgeColor",colors3,"MarkerSize",8)
hold on
line([log10(0.1) C1], [0.6 0.6], "LineWidth",3,"Color",colors3,"LineStyle","--" )
hold on
line([C1 C1], [0 0.6],"LineWidth",3,"Color",colors3,"LineStyle","--" )
hold on
plot(C2,0.7,'o',"MarkerFaceColor",colors3,"MarkerEdgeColor",colors3,"MarkerSize",8)
hold on
line([log10(0.1) C2], [0.7 0.7], "LineWidth",3,"Color",colors3,"LineStyle","--" )
hold on
line([C2 C2], [0 0.7],"LineWidth",3,"Color",colors3,"LineStyle","--" )
ylim([0 1.1])
xlim([log10(0.1) log10(1000)])
yticks([E01 E02 Emax1])
yticklabels(["E_{01}" "E_{02}" "E_1\cdotE_2"])
xticks([C1 C2])
xticklabels(["C1" "C2"])
axis square
xlabel("Log(agent)")

tl.YLabel.String = "Effect";
fontsize(20,"points") 

%Set Line width
H=gca;
H.LineWidth=1.5;

file_name = "hill_vs_multiphasic.pdf";
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
