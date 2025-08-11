
% Function to fit the surface to data, perform Monte Carlo simulations and
% make the figure. 
% Inputs:
    % p : 'struct' object containing 
        % map: own color for the surface display
        % scaled_conc: cytokines scaled concentration from Yale et al. (2021)
        % monotonic: 
            % if monotonic == 0: drugs 1 and 2 are non-monotonic
            % if monotonic == 1: drug 1 is monotonic and drug 2 is
            % non-monotonic
            % if monotonic == 2: drug 1 is non-monotonic and drug 2 is
            % monotonic
            % if monotonic ==3: drugs 1 and 2 are monotonic
       % days: indicates days of the experiment
       % cell: indicates cell type
       % drug1_name and drug2_name: names of the cytokines in the
       % interaction
       % E0, E1, E2, E3: initial guesses for the fit
       % idx1 and idx2: contains the indices of the cytokines in the cell
       % array with the names
       % c01 and c02: contains the converting values for the scaled
       % concentrations from Yale et al. (2021)
       % int_name: name of the interaction
       % drug1_conc and drug2_conc: concentrations of the cytokines in
       % ng/mL
       % drug1_coeff and drug2_coeff: polynomial coefficients of each
       % cytokine
       % int_coeff: polynomial coefficients for the interaction
    %results: results from the single fits to set upper and lower bounds
            %and initial guesses for the fit

%OUTPUTS:
    % fit_parameters: vector containing best-fit parameters for the interaction
    % resnorm: residual norm from the surface fit
    % resnorm_normalized : residual norm from the surface fit normalized by
    % the number of data points
    % parameters_MC : fit_parameters from the Monte Carlo simulations
    % CI: confidence intervals following the Monte Carlo simulations
    % output: output from the lsqnonlin MATLAB function for the surface
    % fitting
    % iterationData: shows the iteration progression from the fit

     

function [fit_parameters,resnorm_normalized,parameters_MC,CI,output,iterationData] = runCode_fun(p,results)

    %Evaluate polynomial
    [data,nanNbr] = polySurface_fun(p);

    N = size(data,1)*size(data,2); %Number of data points
    N = N-nanNbr; %Number data points that are not NaN
    
    %Perform fit
    [fit_parameters,resnorm,~,~,output,iterationData] = fit_surface(p,data,results.(p.days).(p.cell));
    resnorm_normalized = resnorm/N;
    disp('Fit done')

    %Perform Monte Carlo
    rms = sqrt(resnorm/N);
    [CI,parameters_MC] = monteCarlo_fun(p,fit_parameters,rms,results.(p.days).(p.cell));
    disp('Monte Carlo done')

    %Make figure
    [~,~] = figure_fun(p,data,fit_parameters);

end 