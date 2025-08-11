%This function performs the fit on the data.

%INPUTS
    % X is a vector of drug concentrations
    % Y is a vector of drug effect
    % p is a 'struct' object storing parameter guesses and lower and upper
    % bounds

%OUTPUTS
    % param: contains the fit parameters
    % fit_results1 and fit_results2: contain the fitting results as follows:
        % residual norm
        % residual vector
        % exit flag
        % output of the fitting algorithm
        % lambda of the fitting algorithm
        % jacobian of the fit

    function [param,fit_results1,fit_results2] = fit_fun(X,Y,p)

    param_guess1 = p.param_guess1;
    param_guess2 = p.param_guess2;

    lb1 = p.lb1;
    ub1 = p.ub1;
    lb2 = p.lb2;
    ub2 = p.ub2;

    %setting the optimisation routine specifics
        options = optimoptions(@lsqnonlin, ...
          'FunctionTolerance',1e-6,'StepTolerance',1e-6,...,
          'OptimalityTolerance', ...
             1e-6,'display','off','Algorithm','levenberg-marquardt'); 

    %Find index at which concentration is maximal
    idx = find(Y==max(Y));

    if (idx == length(Y) || idx ==1) %If the dose-response curve is monotonic
        X1 = X;
        Y1 = Y;

        [fit_parameters,resnorm1,residual1,exitflag1,output1,lambda1,jacobian1] = lsqnonlin(@residualsfunction1, param_guess1, lb1, ub1, options); % Invoking optimiser
        
        %Save results of the fit
        fit_results1.resnorm = resnorm1;
        fit_results1.residual = residual1;
        fit_results1.exitflag = exitflag1;
        fit_results1.output = output1;
        fit_results1.lambda = lambda1;
        fit_results1.jacobian = jacobian1;
        fit_results2 = [];
        
        %Define best-fit parameters in vector 'param'
        param = fit_parameters;
    
    else %When the dose-response is not monotonic (perform fit twice)
        
        X1 = X(1:idx);
        X2 = X(idx:end);

        Y1 = Y(1:idx);
        Y2 = Y(idx:end);

        [fit_parameters1,resnorm1,residual1,exitflag1,output1,lambda1,jacobian1] = lsqnonlin(@residualsfunction1, param_guess1,lb1, ub1, options);   % Invoking optimiser

        [fit_parameters2,resnorm2,residual2,exitflag2,output2,lambda2,jacobian2] = lsqnonlin(@residualsfunction2, param_guess2,lb2, ub2, options);   % Invoking optimiser

        fit_results1.resnorm = resnorm1;
        fit_results1.residual = residual1;
        fit_results1.exitflag = exitflag1;
        fit_results1.output = output1;
        fit_results1.lambda = lambda1;
        fit_results1.jacobian = jacobian1;
        
        fit_results2.resnorm = resnorm2;
        fit_results2.residual = residual2;
        fit_results2.exitflag = exitflag2;
        fit_results2.output = output2;
        fit_results2.lambda = lambda2;
        fit_results2.jacobian = jacobian2;

        param = [fit_parameters1, fit_parameters2];

    end 

%------------------------------------------------------------------------
function val = residualsfunction1(param)
    
    E01 = param(1);
    C1 = param(2);
    h1 = param(3);
    Emax1 = param(4);
    
    E = @(d) E01 + (Emax1-E01)./(1+10.^((C1-d)*h1));
    eff = E(X1);
   
    val = reshape(abs(Y1 - eff),[],1);
end

function val = residualsfunction2(param)
    
    E02 = param(1);
    C2 = param(2);
    h2 = param(3);
    Emax2 = param(4);
    
    E = @(d) E02 + (Emax2-E02)./(1+10.^((C2-d)*h2));
    eff = E(X2);
   
    val = reshape(abs(Y2 - eff),[],1);
end

%------------------------------------------------------------------------
end