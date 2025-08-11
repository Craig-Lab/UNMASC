%This function performs the Monte Carlo simulations for the surfaces

%INPUTS:
    % p: 'struct' object containing
        % drug1_conc: concentration of drug 1
        % drug2_conc: concentration of drug 2
    % param : best-fit parameters from the surface fit
    % resnorm : residual norm from the surface fit
    % results : results from the single fits to set upper and lower bounds
    % for the fit
    
%OUTPUTS:
    % CI: matrix containing the 95% confidence interval for each parameter of
    % the surface. First line is the lower bound and second line is the
    % upper bound. Each column is a different parameter.
    % new_param: n x length(param) matrix. Each row is a set of parameters
    % estimated by the Monte Carlo simulation method. Each column is a
    % different parameter.

function [CI,new_param] = monteCarlo_fun(p,param,resnorm,results)
    
    % Number of times to perform the fit
    n = 200; 
    
    % Set drug concentrations
    drug1 = p.drug1_conc;
    drug2 = p.drug2_conc;

    % Best fit data
    E = surface_fun(drug1,drug2,param,p);
    
    % Initialize matrix to store new parameter values
    new_param = zeros(n,length(param));

    for k = 1:n
        % Generate best fit data + noise
        new_E = zeros(size(E,1),size(E,2));
    
        for i =1:size(E,1)
            for j = 1:size(E,2)
                new_E(i,j) = E(i,j) + normrnd(0,resnorm);
            end 
        end 
    
        %Perform fit
        new_fit_parameters = fit_surface_monteCarlo(p,new_E,results,param);

        %Store new parameters in matrix
        new_param(k,:) = new_fit_parameters;
    end 

    %Create matrix to store lower and upper bounds of the confidence
    %interval
    CI = zeros(2,length(param));

     %Create confidence intervals for the parameters
    for i = 1:length(param)
        S = std(new_param(:,i));

        CI(1,i) = param(i)-2*S; %lower bound
        CI(2,i) = param(i)+2*S; %upper bound
    end 
