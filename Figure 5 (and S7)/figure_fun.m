% Function to create a figure of the surface fit

%INPUTS:
    % p: 'struct' object containing
        % drug1_conc and drug2_conc: concentrations of cytokines 1 and 2
        % respectively (in ng/mL)
        % map.map: own color for the surface display
        % monotonic: 
            % if monotonic == 0: drugs 1 and 2 are non-monotonic
            % if monotonic == 1: drug 1 is monotonic and drug 2 is
            % non-monotonic
            % if monotonic == 2: drug 1 is non-monotonic and drug 2 is
            % monotonic
            % if monotonic ==3: drugs 1 and 2 are monotonic
        % data: data points used for the fit
        % fit_parameters: best-fit parameters for the surface

%OUTPUTS: 
    % fig: figure of the surface
    % E: matrix containing the points corresponding to the surface
    % evaluated at each data point concentrations 


function [fig,E] = figure_fun(p,data,fit_parameters)

    param = fit_parameters;

    %Define drug concentration vectors of length 500 for surface
    drug1 = linspace(p.drug1_conc(1),p.drug1_conc(end),500);
    drug2 = linspace(p.drug2_conc(1),p.drug2_conc(end),500);

    %Evaluate points on surface
    E = surface_fun(drug1,drug2,param,p);

    %Define meshgrid for the data points
    [X,Y] = meshgrid(log10(p.drug1_conc),log10(p.drug2_conc));

    %Define meshgrid for the surface
    [X2,Y2] = meshgrid(log10(drug1), log10(drug2));
    
        alphaVal = 0.7;

        fig = figure('Visible','off');
        x0=10;
        y0=10;
        width=400;
        height=300;
        set(fig,'position',[x0,y0,width,height])

        surf(X2,Y2,E,'EdgeColor','none','FaceAlpha',alphaVal);
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
        
        %Set Line width
        H=gca;
        H.LineWidth=1.5;

        %Set alpha for whole figure
        alpha(.7)
        
        pause(2)
        %Define file name
        file_name = strcat("Results/",p.days,"/",p.cell,"/",p.drug1_name,"_",p.drug2_name,".pdf");

        exportgraphics(fig,file_name)

end 