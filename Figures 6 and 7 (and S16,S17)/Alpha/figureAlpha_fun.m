function [] = figureAlpha_fun(p,data)

    n = 200; %number of data points for each distribution

    xgroupdata = [1*ones(n,1), 3*ones(n,1), 5*ones(n,1)]; %xvalues for the violin plot

    x0=0;
    y0=0;
    width=7;
    height=5;
    alpha = 1;
    fs = 9; %font size
    densityWidth = 1;
    

    color1 = "#8B5058";
    color2 = "#525B7A";
    color3 = "#EC9F79";
    color4 = "#A781A2";
    color5 = "#9FC2CC";

    %Create figure

    fig = figure("Visible","off");
    set(fig,'units','centimeters','position',[x0,y0,width,height])
    v = violinplot(xgroupdata,log10(data),FaceAlpha=alpha,DensityWidth=densityWidth);
    hold on
    yline(0,"LineWidth",1,"LineStyle","--","Color","#000000")
    
    if p.alpha12 == 1
        ylabel(["Log-fold change of", "potency of cytokine 2"])
    elseif p.alpha12 == 0
        ylabel(["Log-fold change of", "potency of cytokine 1"])
    end 

    if strcmp(p.days,'D7')==1
        v(1).FaceColor = color1;
        v(2).FaceColor = color2;
        v(3).FaceColor = color3;
        leg = legend("PROT","CD4ISP","CD3-");
        leg.Location = "bestoutside";
        stage = "Differentiation";

    elseif strcmp(p.days,'D14')==1
        v(1).FaceColor = color3;
        v(2).FaceColor = color4;
        v(3).FaceColor = color5;
        leg = legend("CD3-","CD3+","8SP");
        leg.Location = "bestoutside";
        stage = "Maturation";
    end 
    
    xlim([0 6])
    ylim(p.ylim)

    xticks(3)
    yticks([p.ylim(1) 0 p.ylim(2)])
    set(gca,'TickLength',[0 0])
    
    row1 = strcat('Cytokine 1:'," ", p.cyt1);
    row2 = strcat('Cytokine 2:'," ", p.cyt2);
    labelArray = [row1; row2];
    labelArray = strjust(pad(labelArray),'left'); 
    xticklabels(strtrim(sprintf('%s\\newline%s\n', labelArray{:})))
   
    
    fontsize(fs,"points")
    % H=gca;
    % H.LineWidth=1;
    
    
    if p.alpha12 == 1  
        file_name= strcat("Figures/",stage,"_",p.int,"_alpha12.pdf");
        exportgraphics(fig,file_name)
    elseif p.alpha12 == 0
        file_name= strcat("Figures/",stage,"_",p.int,"_alpha21.pdf");
        exportgraphics(fig,file_name)
    end 
end 