%This function is to perform the single fit on the data as well as the
%single fits on all the noise data and create figures. The figure is saved
%as a 'fig' file. 

%INPUTS: 
    % p: 'struct' object that contains
        %k: number attributed to cytokine
        %name: name of the cytokine
        %days: days of the experiment
        %cell: cell type
        %E01,C1,h1,Emax1,E02,C2,h2,Emax2: initial guesses for the fit
        %p.param_guess1: initial guesses for the first Hill function
        %p.param_guess2: initial guesses for the second Hill function
   % noise: noisy data previously generated
   % data: data points to fit
   % conc: concentration of the cytokines (in ng/mL)

%OUTPUTS:
    % results: parameters with the 95% confidence intervals generated using
    % the noisy data.
    % noise_param: all parameters of the fits for the noisy data.
    % resnorm1: residual norm for the first fit (first Hill function)
    % resnorm2: residual norm for the second fit (second Hill function)

    function [results,noise_param,resnorm1,resnorm2] = single_fit(p,noise,data,conc)

%% Perform fits
        Y = data.(p.days).(p.cell).(p.name);
        X = log10(conc(:,p.k))';
        %noise_val = noise.(p.days).(p.cell).(p.name);
        
        %Perform single fit on data
        [param,fit1_results,fit2_results] = fit_fun(X,Y,p);
        
        resnorm1 = fit1_results.resnorm;
        if isempty(fit2_results)
            resnorm2 = [];
        else 
            resnorm2 = fit2_results.resnorm;
        end 

        %Determine if data is monotonic or not
        %If data is monotonic increasing, we select the noisy data that is
        %monotically increasing
        %If data is monotonic decreasing, we select the noisy data that is
        %monotically decreasing
        %If data is multiphasic, we select the noisy data that is
        %multiphasic
        idx = find(Y==max(Y));
        if idx == length(Y) %Data is monotonic increasing
            monotonic = 1;
            noise_val = noise.monoInc.(p.days).(p.cell).(p.name);
            noise_param = zeros(size(noise_val,1),4);
        elseif idx == 1 %Data is monotonic decreasing
            monotonic = 1;
            noise_val = noise.monoDec.(p.days).(p.cell).(p.name);
            noise_param = zeros(size(noise_val,1),4);
        else 
            monotonic = 0;
            noise_val = noise.multi.(p.days).(p.cell).(p.name);
            noise_param = zeros(size(noise_val,1),8);
        end 

       % Perform single fit on all the noisy data

            for i = 1:size(noise_val,1)
                [param_noise] = ...
                    fit_fun(X,noise_val(i,:),p);

                noise_param(i,:) = param_noise;
            end 
        
        CI = zeros(2,size(noise_param,2));

        %Build intervals for parameters
        for i = 1:size(noise_param,2)
            S = std(noise_param(:,i));
            CI(1,i) = param(i) -2*S;
            CI(2,i) = param(i) +2*S;
        end 

        results = [param; CI];
        file_name_results = strcat('Results/',p.days,'/',p.cell,'/',p.name); %Define file name
        save(file_name_results,"results")

%% Make figures

        new_X = linspace(X(1)-1,X(end)+1,100);

        color = "#1468a2";
        colorPts = "#A2142F";
        x0 = 10;
        y0 = 10;
        width = 3;
        height = 3;

        all_noise = noise.(p.days).(p.cell).(p.name);

        err_neg = Y - min(all_noise);
        err_pos = max(all_noise) - Y;

        if monotonic == 1 %if dose-response curve is monotonic
            fig1 = figure("units","centimeters","position",[x0,y0,width,height],'Visible','off');
            plot(X, all_noise,'Color','#d5d5d5','LineWidth',1)
            hold on
            plot(new_X, monotonic_log(param,new_X),'LineWidth',1,'Color',color)
            hold on
            scatter(X,Y,5,"o",'filled','MarkerEdgeColor',colorPts,'MarkerFaceColor',colorPts,'LineWidth',2)
            hold on
            errorbar(X,Y,err_neg,err_pos,'o','MarkerSize',0.5,'MarkerFaceColor',...
                 "#A2142F",'MarkerEdgeColor',"#A2142F",'LineWidth',1,'Color',colorPts)
            hold off
            ylim([0 1.5])
            yticks([0 0.5 1])
            xlim([X(1)-0.01 X(end)+0.01])
            axis square
            title(p.name)
            fontsize(10,"points")


        elseif monotonic == 0
            fig1 = figure("units","centimeters","position",[x0,y0,width,height],'Visible','off');
            plot(X, all_noise,'Color','#d5d5d5','LineWidth',1)
            hold on
            plot(new_X,multiphasic_log(param,new_X),'LineWidth',1,'Color', color)
            hold on
            scatter(X,Y,5,"o",'filled','MarkerEdgeColor',colorPts,'MarkerFaceColor',colorPts,'LineWidth',2)
            hold on
            errorbar(X,Y,err_neg,err_pos,'o','MarkerSize',0.5,'MarkerFaceColor',...
                 "#A2142F",'MarkerEdgeColor',"#A2142F",'LineWidth',1,'Color',colorPts)
            hold off
            ylim([0 1.5])
            yticks([0 0.5 1])
            xlim([X(1)-0.01 X(end)+0.01])
            axis square
            title(p.name)
            fontsize(10,"points")

        end 

        fontsize(10,"points")

        file_name1 = strcat('Results/',p.days,'/',p.cell,'/',p.name,'_Fit','.fig'); %Define file name for fig
        file_name2 = strcat('Results/',p.days,'/',p.cell,'/',p.name,'_Fit','.pdf'); %Define file name for PDF

        savefig(fig1,file_name1); % Saves the figure as fig file
        exportgraphics(fig1,file_name2) %Saves the figure as PDF file

        function E = monotonic_log(param,conc)
            E01 = param(1);
            C1 = param(2);
            h1 = param(3);
            Emax1 = param(4);
        
            Eff = @(d) E01 + (Emax1-E01)./(1+10.^((C1-d)*h1));
            E = Eff(conc);
        end 

        function E = multiphasic_log(param,conc)
            E01 = param(1);
            C1 = param(2);
            h1 = param(3);
            Emax1 = param(4);
        
            E02 = param(5);
            C2 = param(6);
            h2 = param(7);
            Emax2 = param(8);
        
            Eff = @(d) (E01 + (Emax1-E01)./(1+10.^((C1-d)*h1))).*(E02 + (Emax2-E02)./(1+10.^((C2-d)*h2)));
            E = Eff(conc);
        end 

    end 