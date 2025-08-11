function [] = fun_Dist(p,noiseParam,singleFits,combinationFits)

    %Set parameters for the figures
    alphaVal = 0.5;
    color1 = "r";
    color2 = "#1a7a20";
    color3 = "#6863c7";
    color4 = "#8c8c8c"; %grey
    fs = 11; %font size

    x0=10;
    y0=10;
    width=150;
    height=150;
    maxVal = 1000;

    %Generate figures

    %For the first EC50 of drug 1
    fig1 = figure('Visible','off');
        set(fig1,'position',[x0,y0,width,height])
        area([singleFits.(p.days).(p.cell){"C1 lb",p.drug1} singleFits.(p.days).(p.cell){"C1 ub",p.drug1}],...
            [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
        hold on
        hist = histogram(noiseParam.(p.days).(p.cell).(p.drug1)(:,2),'FaceColor',color3,'FaceAlpha',0.8);
        hold on
        xline(singleFits.(p.days).(p.cell){"C1",p.drug1},'LineWidth',3,'Color',color1)
        hold on
        xline(combinationFits.(p.days).(p.cell){"C11",p.int},'LineWidth',3,'Color',color2)
        hold off
        fontsize(fs,"points")
        xlabel("log(C_{11})")       
        H=gca;
        H.LineWidth=1.5;
        ylim([0 max(hist.Values)])

        file_name1 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/C11",".pdf");
        exportgraphics(fig1,file_name1)

    %For the second EC50 of drug 1
    if (p.drug1_mono == 0) %if drug 1 is non-monotonic
        fig2 = figure('Visible','off');
            set(fig2,'position',[x0,y0,width,height])    
            area([singleFits.(p.days).(p.cell){"C2 lb",p.drug1} singleFits.(p.days).(p.cell){"C2 ub",p.drug1}],...
                [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
            hold on
            hist = histogram(noiseParam.(p.days).(p.cell).(p.drug1)(:,6),'FaceColor',color3,'FaceAlpha',0.8);
            hold on
            xline(singleFits.(p.days).(p.cell){"C2",p.drug1},'LineWidth',3,'Color',color1)
            hold on
            xline(combinationFits.(p.days).(p.cell){"C12",p.int},'LineWidth',3,'Color',color2)
            hold off
            fontsize(fs,"points")
            xlabel("log(C_{12})")
            H=gca;
            H.LineWidth=1.5;
            ylim([0 max(hist.Values)])

            file_name2 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/C12",".pdf");
            exportgraphics(fig2,file_name2)
    end 
    
    %For the first EC50 of drug 2
    fig3 = figure('Visible','off');
        set(fig3,'position',[x0,y0,width,height])
        area([singleFits.(p.days).(p.cell){"C1 lb",p.drug2} singleFits.(p.days).(p.cell){"C1 ub",p.drug2}],...
            [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
        hold on
        hist = histogram(noiseParam.(p.days).(p.cell).(p.drug2)(:,2),'FaceColor',color3,'FaceAlpha',0.8);
        hold on
        xline(singleFits.(p.days).(p.cell){"C1",p.drug2},'LineWidth',3,'Color',color1)
        hold on
        xline(combinationFits.(p.days).(p.cell){"C21",p.int},'LineWidth',3,'Color',color2)
        hold off
        fontsize(fs,"points")
        xlabel("log(C_{21})")
        H=gca;
        H.LineWidth=1.5;
        ylim([0 max(hist.Values)])

        file_name3 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/C21",".pdf");
        exportgraphics(fig3,file_name3)
    
    %For the second EC50 of drug 2
    if (p.drug2_mono == 0) %if drug 2 is non-monotonic 
        fig4 = figure('Visible','off');
            set(fig4,'position',[x0,y0,width,height])
            area([singleFits.(p.days).(p.cell){"C2 lb",p.drug2} singleFits.(p.days).(p.cell){"C2 ub",p.drug2}],...
                [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
            hold on
            hist = histogram(noiseParam.(p.days).(p.cell).(p.drug2)(:,6),'FaceColor',color3,'FaceAlpha',0.8);
            hold on
            xline(singleFits.(p.days).(p.cell){"C2",p.drug2},'LineWidth',3,'Color',color1)
            hold on
            xline(combinationFits.(p.days).(p.cell){"C22",p.int},'LineWidth',3,'Color',color2)
            hold off
            fontsize(fs,"points")
            xlabel("log(C_{22})")
            H=gca;
            H.LineWidth=1.5;
            ylim([0 max(hist.Values)])

            file_name4 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/C22",".pdf");
            exportgraphics(fig4,file_name4)
    end 

    %For the first h of drug 1
    fig5 = figure('Visible','off');
        set(fig5,'position',[x0,y0,width,height])
        area([singleFits.(p.days).(p.cell){"h1 lb",p.drug1} singleFits.(p.days).(p.cell){"h1 ub",p.drug1}],...
            [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
        hold on
        hist = histogram(noiseParam.(p.days).(p.cell).(p.drug1)(:,3),'FaceColor',color3,'FaceAlpha',0.8);
        hold on
        xline(singleFits.(p.days).(p.cell){"h1",p.drug1},'LineWidth',3,'Color',color1)
        hold on
        xline(-combinationFits.(p.days).(p.cell){"h11",p.int},'LineWidth',3,'Color',color2)
        hold off
        fontsize(fs,"points")
        xlabel("h_{11}")       
        H=gca;
        H.LineWidth=1.5;
        ylim([0 max(hist.Values)])

        file_name5 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/h11",".pdf");
        exportgraphics(fig5,file_name5)

   %For the second h of drug 1
   if (p.drug1_mono == 0) %if drug 1 is non-monotonic
    fig6 = figure('Visible','off');
        set(fig6,'position',[x0,y0,width,height])
        area([singleFits.(p.days).(p.cell){"h2 lb",p.drug1} singleFits.(p.days).(p.cell){"h2 ub",p.drug1}],...
            [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
        hold on
        hist = histogram(noiseParam.(p.days).(p.cell).(p.drug1)(:,7),'FaceColor',color3,'FaceAlpha',0.8);
        hold on
        xline(singleFits.(p.days).(p.cell){"h2",p.drug1},'LineWidth',3,'Color',color1)
        hold on
        xline(-combinationFits.(p.days).(p.cell){"h12",p.int},'LineWidth',3,'Color',color2)
        hold off
        fontsize(fs,"points")
        xlabel("h_{12}")       
        H=gca;
        H.LineWidth=1.5;
        ylim([0 max(hist.Values)])

        file_name6 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/h12",".pdf");
        exportgraphics(fig6,file_name6)
   end

    %For the first h of drug 2
    fig7 = figure('Visible','off');
        set(fig7,'position',[x0,y0,width,height])
        area([singleFits.(p.days).(p.cell){"h1 lb",p.drug2} singleFits.(p.days).(p.cell){"h1 ub",p.drug2}],...
            [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
        hold on
        hist = histogram(noiseParam.(p.days).(p.cell).(p.drug2)(:,3),'FaceColor',color3,'FaceAlpha',0.8);
        hold on
        xline(singleFits.(p.days).(p.cell){"h1",p.drug2},'LineWidth',3,'Color',color1)
        hold on
        xline(-combinationFits.(p.days).(p.cell){"h21",p.int},'LineWidth',3,'Color',color2)
        hold off
        fontsize(fs,"points")
        xlabel("h_{21}")       
        H=gca;
        H.LineWidth=1.5;
        ylim([0 max(hist.Values)])

        file_name7 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/h21",".pdf");
        exportgraphics(fig7,file_name7)

    %For the second h of drug 2
    if (p.drug2_mono == 0) %if drug 1 is non-monotonic
    fig8 = figure('Visible','off');
        set(fig8,'position',[x0,y0,width,height])
        area([singleFits.(p.days).(p.cell){"h2 lb",p.drug2} singleFits.(p.days).(p.cell){"h2 ub",p.drug2}],...
            [maxVal maxVal], 'FaceAlpha',alphaVal,'FaceColor',color4)
        hold on
        hist = histogram(noiseParam.(p.days).(p.cell).(p.drug2)(:,7),'FaceColor',color3,'FaceAlpha',0.8);
        hold on
        xline(singleFits.(p.days).(p.cell){"h2",p.drug2},'LineWidth',3,'Color',color1)
        hold on
        xline(-combinationFits.(p.days).(p.cell){"h22",p.int},'LineWidth',3,'Color',color2)
        hold off
        fontsize(fs,"points")
        xlabel("h_{22}")       
        H=gca;
        H.LineWidth=1.5;
        ylim([0 max(hist.Values)])

        file_name8 = strcat("Results/",p.days,"/",p.cell,"/",p.int,"/h22",".pdf");
        exportgraphics(fig8,file_name8)
   end