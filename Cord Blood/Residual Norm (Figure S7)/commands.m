clear variables

load("surfacesResNorm.mat")

color1 = "#8B5058";
color2 = "#525B7A";
color3 = "#EC9F79";
color4 = "#A781A2";
color5 = "#9FC2CC";
color6 = "#d17466";

x0=10;
y0=10;
width=10;
height=6;

%% CD7+CD5-

y = zeros(2,3);

y(1,1) = resNorm.vial2.SCF_IL7.CD7pCD5n;
y(1,2) = resNorm.vial3.SCF_IL7.CD7pCD5n;
y(1,3) = resNorm.vial4.SCF_IL7.CD7pCD5n;

y(2,1) = resNorm.vial2.IL3_IL7.CD7pCD5n;
y(2,2) = resNorm.vial3.IL3_IL7.CD7pCD5n;
y(2,3) = resNorm.vial4.IL3_IL7.CD7pCD5n;


x = [1 2];

fig1 = figure(1);
set(fig1,"units","centimeters",'position',[x0,y0,width,height])

b = bar(x,y,'BarWidth',0.8,'GroupWidth',0.4);
row1 = {"SCF" "IL-3"};
row2 = {"IL-7" "IL-7"};
labelArray = [row1; row2];
xticklabels(strtrim(sprintf('%s\\newline%s\n', labelArray{:})))

b(1).FaceColor = color1;
b(2).FaceColor = color2;
b(3).FaceColor = color3;

fontsize(12,"points")
ylabel(["Normalized", "residual norm"])

ax = gca;
ax.YAxis.Exponent = -1;

%% CD7+CD5+

y = zeros(2,3);

y(1,1) = resNorm.vial2.SCF_IL7.CD7pCD5p;
y(1,2) = resNorm.vial3.SCF_IL7.CD7pCD5p;
y(1,3) = resNorm.vial4.SCF_IL7.CD7pCD5p;

y(2,1) = resNorm.vial2.IL3_IL7.CD7pCD5p;
y(2,2) = resNorm.vial3.IL3_IL7.CD7pCD5p;
y(2,3) = resNorm.vial4.IL3_IL7.CD7pCD5p;

x = [1 2];

fig2 = figure(2);
set(fig2,"units","centimeters",'position',[x0,y0,width,height])

b = bar(x,y,'BarWidth',0.8,'GroupWidth',0.4);
row1 = {"SCF" "IL-3"};
row2 = {"IL-7" "IL-7"};
labelArray = [row1; row2];
xticklabels(strtrim(sprintf('%s\\newline%s\n', labelArray{:})))

b(1).FaceColor = color4;
b(2).FaceColor = color5;
b(3).FaceColor = color6;

fontsize(12,"points")
 ylabel(["Normalized", "residual norm"])

ax = gca;
ax.YAxis.Exponent = -1;

%% Export

file_name1 = "resNorm_CD7pCD5n.pdf";
file_name2 = "resNorm_CD7pCD5p.pdf";

exportgraphics(fig1,file_name1)
exportgraphics(fig2,file_name2)