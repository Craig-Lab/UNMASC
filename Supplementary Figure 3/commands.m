%% Initialize
clear
%Load color for surfaces
p.map = load("ownColor.mat");

%Define vectors of drug concentration
drug1 = logspace(log10(0.0001),log10(10000),500);
drug2 = logspace(log10(0.0001),log10(10000),500);

lineColor1 = "#9ebcda";
lineColor2 = "#8c6bb1";
lineColor3 = "#08519c";
lineGrey = "#9c9c9c";
lineWidth = 4;

x0=10;
y0=10;
width = 1500;
height= 300;

%% Both drugs are multiphasic

p.monotonic = 0;

E0 = 1;
E1 = 0.9;
E2 = 0.3;
E3 = 0;

h11 = -1;
h12 = 1;
h21 = -1.5;
h22 = 1.5;

C11 = -2;
C12 = 2;
C21 = -2;
C22 = 1;

alpha12 = 1;
alpha21 = 1;

param = [E0, E1, E2, E3, h11, h12, h21, h22, C11, C12, C21, C22, alpha12, alpha21];

Eff1 = surface_fun(drug1,drug2,param,p);

[X2,Y2] = meshgrid(log10(drug1), log10(drug2));

fig = figure('Visible','off');
set(fig,'position',[x0,y0,width,height])
tiledlayout(1,4,'TileSpacing','loose');

ax1 = nexttile;
surf(X2,Y2,Eff1,'EdgeColor','none','FaceAlpha',0.5);
grid off
colormap(fig,p.map.map)
clim([0 1])
cb = colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1]);
cb.Location = 'westoutside';
cb.Label.String = 'Effect';
cb.Label.VerticalAlignment = "middle";
cb.Label.Rotation = 0;
cb.Label.Position = [0.5 1+0.1];
set(cb,'position',[0.03 0.3 .01 .5])
xticks([])
yticks([])
zticks([])
xlim([log10(drug1(1)) log10(drug1(end))])
ylim([log10(drug2(1)) log10(drug2(end))])
zlim([0 1])
xlh = xlabel("Log(Agent 1)","Rotation",22.5);
xlh.Position(1) = xlh.Position(1)-2; 
xlh.Position(2) = xlh.Position(2)-1;
ylh = ylabel("Log(Agent 2)","Rotation",-32.5);
ylh.Position(1) = ylh.Position(1)-1;
ylh.Position(2) = ylh.Position(2)-2;
zlabel("Effect")
axis square
ax1.LineWidth = 2;

ax2 = nexttile;
ax2.LineWidth = 2;
line([C11 C11], [0 E2/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([C12 C12], [0 E2/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C11], [E2/2 E2/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C12], [E2/2 E2/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug1),Eff1(1,:),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([0 E2+0.1])
xticks([C11 C12])
xticklabels({"C_{11}","C_{12}"})
yticks([E3 0.5*E2 E2])
xlabel("Log(Agent 1)")
ylabel("Effect")
axis square

ax3 = nexttile;
ax3.LineWidth = 2;
line([C21 C21], [0 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([C22 C22], [0 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C21], [E1/2 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C22], [E1/2 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug2),Eff1(:,1),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([0 E1+0.1])
xticks([C21 C22])
xticklabels({"C_{21}","C_{22}"})
yticks([E3 0.5*E1 E1])
xlabel("Log(Agent 2)")
ylabel("Effect")
axis square

nexttile;
textStr = ["Effect parameters";strcat('E0 = '," ",num2str(E0)); strcat('E1 = '," ",num2str(E1));strcat('E2 = '," ",num2str(E2));strcat('E3 = '," ",num2str(E3))];
text(0,0.5, textStr,'HorizontalAlignment','left',...
    'VerticalAlignment','middle')
axis off

fontsize(14,"points")

file_name1 = "bothBiphasicID.pdf";
exportgraphics(fig,file_name1)

%% First drug is monotonic and second drug is biphasic

p.monotonic = 1;

E0 = 1;
E1 = 0.5;
E2 = 0.8;
E3 = 0;

h11 = -1.5;
h21 = -2;
h22 = 1.5;

C11 = 0.3;
C21 = -2;
C22 = 1;

alpha12 = 1;
alpha21 = 1;

param = [E0, E1, E2, E3, h11, h21, h22, C11, C21, C22, alpha12, alpha21];

Eff2 = surface_fun(drug1,drug2,param,p);

[X2,Y2] = meshgrid(log10(drug1), log10(drug2));

fig = figure('Visible','off');
set(fig,'position',[x0,y0,width,height])
tiledlayout(1,4,'TileSpacing','loose');

ax1 = nexttile;
surf(X2,Y2,Eff2,'EdgeColor','none','FaceAlpha',0.5);
grid off
colormap(fig,p.map.map)
clim([0 1])
cb = colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1]);
cb.Location = 'westoutside';
cb.Label.String = 'Effect';
cb.Label.VerticalAlignment = "middle";
cb.Label.Rotation = 0;
cb.Label.Position = [0.5 1+0.1];
set(cb,'position',[0.03 0.3 .01 .5])
xticks([])
yticks([])
zticks([])
xlim([log10(drug1(1)) log10(drug1(end))])
ylim([log10(drug2(1)) log10(drug2(end))])
zlim([0 1])
xlh = xlabel("Log(Agent 1)","Rotation",22.5);
xlh.Position(1) = xlh.Position(1)-2; 
xlh.Position(2) = xlh.Position(2)-1;
ylh = ylabel("Log(Agent 2)","Rotation",-32.5);
ylh.Position(1) = ylh.Position(1)-1;
ylh.Position(2) = ylh.Position(2)-2;
zlabel("Effect")
axis square
ax1.LineWidth = 2;

ax2 = nexttile;
ax2.LineWidth = 2;
line([C11 C11], [0 E2/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C11], [E2/2 E2/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug1),Eff2(1,:),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([0 E2+0.1])
xticks(C11)
xticklabels({"C_{11}","C_{12}"})
yticks([E3 0.5*E2 E2])
xlabel("Log(Agent 1)")
ylabel("Effect")
axis square

ax3 = nexttile;
ax3.LineWidth = 2;
line([C21 C21], [0 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([C22 C22], [0 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C21], [E1/2 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C22], [E1/2 E1/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug2),Eff2(:,1),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([0 E1+0.1])
xticks([C21 C22])
xticklabels({"C_{21}","C_{22}"})
yticks([E3 0.5*E1 E1])
xlabel("Log(Agent 2)")
ylabel("Effect")
axis square

nexttile;
textStr = ["Effect parameters";strcat('E0 = '," ",num2str(E0)); strcat('E1 = '," ",num2str(E1));strcat('E2 = '," ",num2str(E2));strcat('E3 = '," ",num2str(E3))];
text(0,0.5, textStr,'HorizontalAlignment','left',...
    'VerticalAlignment','middle')
axis off

fontsize(14,"points")

file_name2 = "monotonicBiphasic.pdf";
exportgraphics(fig,file_name2)


%% Dose-responses are both decreasing and then increasing
p.monotonic = 0;

E0 = 0;
E1 = 0.4;
E2 = 0.1;
E3 = 1;

h11 = -1;
h12 = 2;
h21 = -1;
h22 = 2;

C11 = -1;
C12 = 2;
C21 = -1;
C22 = 2;

alpha12 = 1;
alpha21 = 1;

param = [E0, E1, E2, E3, h11, h12, h21, h22, C11, C12, C21, C22, alpha12, alpha21];

Eff3 = surface_fun(drug1,drug2,param,p);

[X2,Y2] = meshgrid(log10(drug1), log10(drug2));

fig = figure('Visible','off');
set(fig,'position',[x0,y0,width,height])
tiledlayout(1,4,'TileSpacing','loose');

ax1 = nexttile;
surf(X2,Y2,Eff3,'EdgeColor','none','FaceAlpha',0.5);
grid off
colormap(fig,p.map.map)
clim([0 1])
cb = colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1]);
cb.Location = 'westoutside';
cb.Label.String = 'Effect';
cb.Label.VerticalAlignment = "middle";
cb.Label.Rotation = 0;
cb.Label.Position = [0.5 1+0.1];
set(cb,'position',[0.03 0.3 .01 .5])
xticks([])
yticks([])
zticks([])
xlim([log10(drug1(1)) log10(drug1(end))])
ylim([log10(drug2(1)) log10(drug2(end))])
zlim([0 1])
xlh = xlabel("Log(Agent 1)","Rotation",22.5);
xlh.Position(1) = xlh.Position(1)-2; 
xlh.Position(2) = xlh.Position(2)-1;
ylh = ylabel("Log(Agent 2)","Rotation",-32.5);
ylh.Position(1) = ylh.Position(1)-1;
ylh.Position(2) = ylh.Position(2)-2;
zlabel("Effect")
axis square
ax1.LineWidth = 2;

ax2 = nexttile;
line([C11 C11], [0 1-(1-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([C12 C12], [0 1-(1-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C11], [1-(1-E2)/2 1-(1-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C12], [1-(1-E2)/2 1-(1-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug1),Eff3(1,:),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([0 1])
xticks([C11 C12])
xticklabels({"C_{11}","C_{12}"})
yticks([E2 1-(1-E2)/2 E3])
xlabel("Log(Agent 1)")
ylabel("Effect")
axis square
ax2.LineWidth = 2;

ax3 = nexttile;
line([C21 C21], [0 1-(1-E1)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([C22 C22], [0 1-(1-E1)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C21], [1-(1-E1)/2 1-(1-E1)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C22], [1-(1-E1)/2 1-(1-E1)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug2),Eff3(:,1),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([0 1])
xticks([C21 C22])
xticklabels({"C_{21}","C_{22}"})
yticks([E1 1-(1-E1)/2 E3])
xlabel("Log(Agent 2)")
ylabel("Effect")
axis square
ax3.LineWidth = 2;

nexttile;
textStr = ["Effect parameters";strcat('E0 = '," ",num2str(E0)); strcat('E1 = '," ",num2str(E1));strcat('E2 = '," ",num2str(E2));strcat('E3 = '," ",num2str(E3))];
text(0,0.5, textStr,'HorizontalAlignment','left',...
    'VerticalAlignment','middle')
axis off

fontsize(14,"points")

file_name3 = "bothBiphasicDI.pdf";
exportgraphics(fig,file_name3)


%% Dose-responses have different increasing and decreasing phases

p.monotonic = 0;

E0 = 1;
E1 = 0.8;
E2 = 0.1;
E3 = 0.3;

h11 = -1;
h12 = 2;
h21 = -2;
h22 = 1;

C11 = -1;
C12 = 2;
C21 = -2;
C22 = 2;

alpha12 = 1;
alpha21 = 1;

param = [E0, E1, E2, E3, h11, h12, h21, h22, C11, C12, C21, C22, alpha12, alpha21];

Eff4 = surface_fun(drug1,drug2,param,p);

[X2,Y2] = meshgrid(log10(drug1), log10(drug2));

fig = figure('Visible','off');
set(fig,'position',[x0,y0,width,height])
tiledlayout(1,4,'TileSpacing','loose');

ax1 = nexttile;
surf(X2,Y2,Eff4,'EdgeColor','none','FaceAlpha',0.5);
grid off
colormap(fig,p.map.map)
clim([0 1])
cb = colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1]);
cb.Location = 'westoutside';
cb.Label.String = 'Effect';
cb.Label.VerticalAlignment = "middle";
cb.Label.Rotation = 0;
cb.Label.Position = [0.5 1+0.1];
set(cb,'position',[0.03 0.3 .01 .5])
xticks([])
yticks([])
zticks([])
xlim([log10(drug1(1)) log10(drug1(end))])
ylim([log10(drug2(1)) log10(drug2(end))])
zlim([0 1])
xlh = xlabel("Log(Agent 1)","Rotation",22.5);
xlh.Position(1) = xlh.Position(1)-2; 
xlh.Position(2) = xlh.Position(2)-1;
ylh = ylabel("Log(Agent 2)","Rotation",-32.5);
ylh.Position(1) = ylh.Position(1)-1;
ylh.Position(2) = ylh.Position(2)-2;
zlabel("Effect")
axis square
ax1.LineWidth = 2;

ax2 = nexttile;
line([C11 C11], [0 E3-(E3-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([C12 C12], [0 E3-(E3-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C11], [E3-(E3-E2)/2 E3-(E3-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C12], [E3-(E3-E2)/2 E3-(E3-E2)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug1),Eff4(1,:),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([E2-0.1 E3+0.1])
xticks([C11 C12])
xticklabels({"C_{11}","C_{12}"})
yticks([E2 E3-(E3-E2)/2 E3])
xlabel("Log(Agent 1)")
ylabel("Effect")
axis square
ax2.LineWidth = 2;

ax3 = nexttile;
line([C21 C21], [0 E3+(E1-E3)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([C22 C22], [0 E3+(E1-E3)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C21], [E3+(E1-E3)/2 E3+(E1-E3)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
line([-4 C22], [E3+(E1-E3)/2 E3+(E1-E3)/2],"LineWidth",lineWidth,'Color',lineGrey,'LineStyle',"--")
hold on
plot(log10(drug2),Eff4(:,1),'Color',lineColor3,'LineWidth',lineWidth)
xlim([-4 4])
ylim([0 E1+0.1])
xticks([C21 C22])
xticklabels({"C_{21}","C_{22}"})
yticks([E3 E3+(E1-E3)/2 E1])
xlabel("Log(Agent 2)")
ylabel("Effect")
axis square
ax3.LineWidth = 2;

nexttile;
textStr = ["Effect parameters";strcat('E0 = '," ",num2str(E0)); strcat('E1 = '," ",num2str(E1));strcat('E2 = '," ",num2str(E2));strcat('E3 = '," ",num2str(E3))];
text(0,0.5, textStr,'HorizontalAlignment','left',...
    'VerticalAlignment','middle')
axis off

fontsize(14,"points")

file_name4 = "differentBehaviours.pdf";
exportgraphics(fig,file_name4)
