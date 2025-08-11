%Function to generate noisy data based on the polynomial coefficients.

%INPUTS:
    % y : function handle for the poylnomial
    % conc : vector of scaled concentration
    % coeff : coefficients of the function handle
    % se : standard error of the coefficients of the function handle
%OUTPUTS: 
    %noise_data : matrix filled with the generated noise. Each row is a new
    %set of noisy data and each column is a concentration.

function [noise_data] = noise_fun(y,conc,coeff,se)
    n = 1000; %Number of times to generate new coefficients
    
    y_max = max(y(conc).^2); %Square because the data is the square root of the number of cells.
    
    new_coeff = zeros(n,length(coeff)); %Initialize matrix that will store new coefficients
    
    k=1;
    while (k<=n)
            for j = 1:length(coeff)
                test(j) = coeff(j) + normrnd(0,se(j));
            end 
            
            for j = 1:length(coeff)
                 if (test(j) < coeff(j) - 2*se(j) || test(j) > coeff(j) + 2*se(j))
                     test(j)=-1; %set to -1 if it doesn't fall within plus/minus 2 standard deviations
                 end 
            end 

            if any(test==-1)
                k=k; %if coefficient doesn't fall within the interval, stay at k
            else
                new_coeff(k,:) = test; 
                k = k+1; %if coefficient falls within the interval, iterate
            end 
    end 

    noise_data = zeros(n,length(conc)); %Initialize matrix with the noise data.

    for i = 1:n
         if length(coeff) == 2 %Case when the polynomial is of order one
            new_fun = @(x) new_coeff(i,1) + new_coeff(i,2)*x;
        elseif length(coeff) == 3 %Case when the polynomial is of order two
            new_fun = @(x) new_coeff(i,1) + new_coeff(i,2)*x + new_coeff(i,3)*x.^2;
        elseif length(coeff) == 4 %Case when the polynomial is of order three
            new_fun = @(x) new_coeff(i,1) + new_coeff(i,2)*x + new_coeff(i,3)*x.^2 + new_coeff(i,4)*x.^3;
        end 

        noise_data(i,:) = new_fun(conc).^2/y_max;
    end 

end 