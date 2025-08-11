%This code creates figures that display how the alpha parameter can impact
%the dose-response curve

drug1 = log10(logspace(log10(0.001),log10(100000),1500));
drug2 = log10(logspace(log10(0.001),log10(100000),1500));

E01 = 1;
h1 = -2;
C1 = -0.5;
E1 = 0;

E02 = 1;
h2 = 2;
C2 = 2.5;
E2 = 0;

alpha21 = [0.1, 1, 10];

param1 = [E01, h1, C1, E1, E02, h2, C2, E2, alpha21(1)];
param2 = [E01, h1, C1, E1, E02, h2, C2, E2, alpha21(2)];
param3 = [E01, h1, C1, E1, E02, h2, C2, E2, alpha21(3)];

param4 = [E01, h1, C1, E1, E02, h2, C2, E2, alpha21(1)];
param5 = [E01, h1, C1, E1, E02, h2, C2, E2, alpha21(2)];
param6 = [E01, h1, C1, E1, E02, h2, C2, E2, alpha21(3)];

[~,~,curve1] = multi(drug1,param1);
[~,~,curve2] = multi(drug1,param2);
[~,~,curve3] = multi(drug1,param3);

[~,~,curve4] = multi(drug1,param4);
[~,~,curve5] = multi(drug1,param5);
[~,~,curve6] = multi(drug1,param6);

%set parameters for figure
lineWidth = 4;
lineColor1 = "#08519c";
lineColor2 = "#6baed6";
lineColor3 = "#c6dbef";

fig = figure('Visible','off');
x0=10;
y0=10;
width=1500;
height=500;
set(fig,'position',[x0,y0,width,height])

tl = tiledlayout(1,2,"TileSpacing","tight");

ax1= nexttile;
plot(drug1,curve1,'LineWidth',lineWidth,'Color',lineColor1)
hold on
plot(drug1,curve2,'LineWidth',lineWidth,'Color',lineColor2)
hold on
plot(drug1,curve3,'LineWidth',lineWidth,'Color',lineColor3)
hold off
xlabel("Log(agent 1)")
ylabel("Effect")
name1 = strcat("\alpha_{21} = ", num2str(alpha21(1)));
name2 = strcat("\alpha_{21} = ", num2str(alpha21(2)));
name3 = strcat("\alpha_{21} = ", num2str(alpha21(3)));

leg1 = legend(name1, name2, name3);
leg1.Location = "northeastoutside";
ylim([0 1.1])
xlim([-3 5])
axis square
ax1.LineWidth = 1.5;

ax2 = nexttile;
plot(drug1,1-curve4,'LineWidth',lineWidth,'Color',lineColor1)
hold on
plot(drug1,1-curve5,'LineWidth',lineWidth,'Color',lineColor2)
hold on
plot(drug1,1-curve6,'LineWidth',lineWidth,'Color',lineColor3)
hold off
xlh2 = xlabel("Log(agent 1)");
ylabel("Effect")
name1 = strcat("\alpha_{21} = ", num2str(alpha21(1)));
name2 = strcat("\alpha_{21} = ", num2str(alpha21(2)));
name3 = strcat("\alpha_{21} = ", num2str(alpha21(3)));
ax2.LineWidth = 1.5;

leg2 = legend(name1, name2, name3);
leg2.Location = "northeastoutside";
ylim([0 1.1])
xlim([-3 5])
axis square

fontsize(18,"points")

file_name = "alpha_example.pdf";
 exportgraphics(fig,file_name)

%--------------------------------------------------------------------------
function [s1,s2,d] = multi(drug,param)
    E01 = param(1);
    h1 = param(2);
    C1 = param(3);
    E1 = param(4);

    E02 = param(5);
    h2 = param(6);
    C2 = param(7);
    E2 = param(8);

    alpha21 = param(9);

    single1 = @(d) E01 + (E1-E01)./(1+10.^((C1-d-log10(alpha21))*h1));
    single2 = @(d) E02 + (E2-E02)./(1+10.^((C2-d-log10(alpha21))*h2));
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
