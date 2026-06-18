 %This function performs the Monte Carlo simulations for the surfaces
 function [CI,new_param] = monteCarlo_fun(x1,x2,p,param,resnorm,n,ci1,ci2)
    
    % Set drug concentrations
    drug1 = x1;
    drug2 = x2;

    % Best fit data
    E = surface_fun(drug1,drug2,param,p);
    
    % Initialize storage
    new_param = NaN(n,length(param));
    exitflags = NaN(n,1);
    iters = NaN(n,1);


    for k = 1:n
        % Generate noisy synthetic data
        new_E = E + resnorm .* randn(size(E));

        start_param = param;

        % Perform fit
        [new_fit_parameters,~,~,out,exitflag] = ...
            fit_surface_monteCarlo(x1,x2,p,new_E,start_param,ci1,ci2);

        exitflags(k) = exitflag;
        iters(k) = out.iterations;

        if exitflag > 0 && all(isfinite(new_fit_parameters))
            new_param(k,:) = new_fit_parameters;
        end
    end 

    good = exitflags > 0 & all(isfinite(new_param),2);
    new_param = new_param(good,:);

    % Create CI
    CI = NaN(2,length(param));

    if isempty(new_param)
        warning('No successful Monte Carlo fits.');
        return
    end

    for i = 1:length(param)
        S = std(new_param(:,i), 'omitnan');
        CI(1,i) = param(i)-2*S;
        CI(2,i) = param(i)+2*S;
    end 
end

