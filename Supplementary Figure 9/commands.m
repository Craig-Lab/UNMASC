clear variables

load("resnorm.mat")

color1 = "#8B5058";
color2 = "#525B7A";
color3 = "#EC9F79";
color4 = "#A781A2";
color5 = "#9FC2CC";

x0=10;
y0=10;
width=10;
height=6;


%% D7-14

y = zeros(4,3);

y(1,1) = resnorm.D7.PROT.SCF_IL7;
y(1,2) = resnorm.D7.CD4ISP.SCF_IL7;
y(1,3) = resnorm.D7.CD3N.SCF_IL7;

y(2,1) = resnorm.D7.PROT.IL3_IL7;
y(2,2) = resnorm.D7.CD4ISP.IL3_IL7;
y(2,3) = resnorm.D7.CD3N.IL3_IL7;

y(3,1) = resnorm.D7.PROT.IL7_TNFa;
y(3,2) = resnorm.D7.CD4ISP.IL7_TNFa;
y(3,3) = resnorm.D7.CD3N.IL7_TNFa;

y(4,1) = resnorm.D7.PROT.TNFa_CXCL12;
y(4,2) = resnorm.D7.CD4ISP.TNFa_CXCL12;
y(4,3) = resnorm.D7.CD3N.TNFa_CXCL12;

x = [1 2 3 4];

fig1 = figure(1);
set(fig1,"units","centimeters",'position',[x0,y0,width,height])

b = bar(x,y,'BarWidth',0.8,'GroupWidth',0.4);
% xticklabels({"SCF-IL7", "IL3-IL7","IL7-TNFa","TNFa-CXCL12"})
row1 = {"SCF" "IL-3" "IL7" "TNFa"};
row2 = {"IL-7" "IL-7" "TNFa" "CXCL12"};
labelArray = [row1; row2];
% labelArray = strjust(pad(labelArray),'left'); 
xticklabels(strtrim(sprintf('%s\\newline%s\n', labelArray{:})))

b(1).FaceColor = color1;
b(2).FaceColor = color2;
b(3).FaceColor = color3;

fontsize(12,"points")
% leg = legend("PROT","CD4ISP","CD3-");
% leg.Location = "bestoutside";
 ylabel(["Normalized", "residual norm"])

ax = gca;
ax.YAxis.Exponent = -3;

%% D14-28

y = zeros(3,3);

y(1,1) = resnorm.D14.CD3N.SCF_IL7;
y(1,2) = resnorm.D14.CD3P.SCF_IL7;
y(1,3) = resnorm.D14.SP.SCF_IL7;

y(2,1) = resnorm.D14.CD3N.IL3_IL7;
y(2,2) = resnorm.D14.CD3P.IL3_IL7;
y(2,3) = resnorm.D14.SP.IL3_IL7;

y(3,1) = resnorm.D14.CD3N.IL7_TNFa;
y(3,2) = resnorm.D14.CD3P.IL7_TNFa;
y(3,3) = resnorm.D14.SP.IL7_TNFa;

x = [1 2 3];

fig2 = figure(2);
set(fig2,"units","centimeters",'position',[x0,y0,width,height])
b = bar(x,y,'BarWidth',0.8,'GroupWidth',0.4);
% xticklabels({"SCF-IL7", "IL3-IL7","IL7-TNFa"})
row1 = {"SCF" "IL-3" "IL7"};
row2 = {"IL-7" "IL-7" "TNFa"};
labelArray = [row1; row2];
% labelArray = strjust(pad(labelArray),'left'); 
xticklabels(strtrim(sprintf('%s\\newline%s\n', labelArray{:})))

b(1).FaceColor = color3;
b(2).FaceColor = color4;
b(3).FaceColor = color5;

fontsize(12,"points")
% leg = legend("CD3-","CD3+","8SP");
% leg.Location = "bestoutside";
 ylabel(["Normalized", "residual norm"])

ax = gca;
ax.YAxis.Exponent = -3;

%% Export

file_name1 = "Results/resnorm_D7.pdf";
file_name2 = "Results/resnorm_D14.pdf";

exportgraphics(fig1,file_name1)
exportgraphics(fig2,file_name2)
