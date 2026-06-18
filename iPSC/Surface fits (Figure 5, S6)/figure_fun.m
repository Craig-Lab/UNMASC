
function [fig,E] = figure_fun(drug1, drug2, data,fit_parameters , p)

    param = fit_parameters;

    [X,Y] = meshgrid(log10(drug1), log10(drug2));

    nGrid = 100;

    drug1_grid = logspace(log10(drug1(1)), log10(drug1(end)), nGrid);
    drug2_grid = logspace(log10(drug2(1)), log10(drug2(end)), nGrid);

    [X1, Y1] = meshgrid(log10(drug1_grid), log10(drug2_grid));

    E = surface_fun(drug1_grid,drug2_grid,param ,p);

    
        alphaVal = 0.7;

        fig = figure('Visible','off');
        x0=10;
        y0=10;
        width=400;
        height=300;
        set(fig,'position',[x0,y0,width,height])

        surf(X1,Y1,E,'EdgeColor','none','FaceAlpha',alphaVal);
        grid off
        colormap(fig,p.map.map)
        clim([0 1])
        cb = colorbar;   
        cb.Limits = [0 1];
        cb.Ticks = [0 0.2 0.4 0.6 0.8 1];
        cb.Label.String = "Effect";
        cb.Label.VerticalAlignment = "middle";
        cb.Label.Rotation = 0;
        cb.Label.Position = [0.5 1+0.1];
        set(cb,'position',[0.8 0.1 .02 .5])
        hold on
        scatter3(X,Y,data,20,"o",'filled','MarkerEdgeColor','#000000','MarkerFaceColor','#000000','LineWidth',1)
        hold off
        zlim([0 1.2])
        zticks([0 0.5 1])
        xlh = xlabel(strcat("Log_{10}(",p.drug1_name,")"),"Rotation",12.5);
        xlh.Position(1) = xlh.Position(1) - 0.40; 
        ylh = ylabel(strcat("Log_{10}(",p.drug2_name,")"),"Rotation",-27.5);
        ylh.Position(2) = ylh.Position(2) - 0.40;
        zlabel("Effect")
        fontsize(11,"points")
        leg = legend('Surface fit','Data');
        leg.Location = "bestoutside";
        axis square
        view(-33.56,17.79)
        
        H=gca;
        H.LineWidth=1.5;

        alpha(.7)


end 