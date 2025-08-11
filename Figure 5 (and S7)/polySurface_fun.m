% Function to generate points from the polynomial coefficients.

%INPUTS:
    %p 'struct' containing:
        % drug1_coeff : vector containing the polynomial coefficients of
        % cytokine 1
        % drug2_coeff: vector containing the polynomial coefficients of
        % cytokine 2

%OUTPUTS:
    % sol: matrix containing the normalized number of cells for the surface
    % nanNbr: number of NaN values for the surface (i.e. surface values
    % that were estimated to be negative with the polynomial)

function [sol,nanNbr] = polySurface_fun(p)
    
    intercept = p.drug1_coeff(1);
    
    syms x1 x2 
    %Create symbolic function
    if ((length(p.drug1_coeff) == 3) && (length(p.drug2_coeff) ==3)) %both quadratic
        y = intercept + p.drug1_coeff(2)*x1 + p.drug1_coeff(3)*x1.^2 + ...
            p.drug2_coeff(2)*x2 +p.drug2_coeff(3)*x2.^2 + p.int_coeff*x1.*x2;
        thirdOrder = 0; %no third order coefficient

    elseif ((length(p.drug1_coeff) == 4) && (length(p.drug2_coeff) == 3)) %first drug polynomial of order 3 and second drug quadratic
        y = intercept + p.drug1_coeff(2)*x1 + p.drug1_coeff(3)*x1.^2 + ...
            p.drug1_coeff(4)*x1.^3 + p.drug2_coeff(2)*x2 + p.drug2_coeff(3)*x2.^2 + ...
            p.int_coeff*x1.*x2;
        thirdOrder = 1; %third order coefficient for drug 1
        firstOrder = 0;

    elseif ((length(p.drug1_coeff) == 3) && (length(p.drug2_coeff) == 4)) %first drug quadratic and second drug polynomial of order 3
        y = intercept + p.drug1_coeff(2)*x1 + p.drug1_coeff(3)*x1.^2 + ...
            p.drug2_coeff(2)*x2 + p.drug2_coeff(3)*x2.^2 + p.drug2_coeff(4)*x2.^3 + ...
            p.int_coeff*x1.*x2;
        thirdOrder = 2; %third order coefficient for drug 2
        firstOrder = 0;

    elseif ((length(p.drug1_coeff) == 2) && (length(p.drug2_coeff)==3)) %drug 1 is monotonic and drug 2 quadratic
        y = intercept + p.drug1_coeff(2)*x1+p.drug2_coeff(2)*x2 + p.drug2_coeff(3)*x2.^2 + ...
            p.int_coeff*x1.*x2;
        firstOrder = 1; %first order coefficient for drug 1
        thirdOrder = 0;

    elseif ((length(p.drug1_coeff) == 3) && (length(p.drug2_coeff)==2)) %drug 1 is quadratic and drug 2 is monotonic
        y = intercept + p.drug1_coeff(2)*x1+p.drug1_coeff(3)*x1.^2 + p.drug2_coeff(2)*x2 + ...
            p.int_coeff*x1.*x2;
        firstOrder = 2; %first order coefficient for drug 2
        thirdOrder = 0;

    elseif ((length(p.drug1_coeff) == 2) && (length(p.drug2_coeff)==4)) %drug 1 is monotonic and second drug polynomial of order 3
        y = intercept + p.drug1_coeff(2)*x1+p.drug2_coeff(2)*x2 + p.drug2_coeff(3)*x2.^2 + ...
            p.drug2_coeff(4)*x2.^3+ p.int_coeff*x1.*x2;
        firstOrder = 1; %first order coefficient for drug 1
        thirdOrder = 2; %third order coefficient for drug 2

    elseif ((length(p.drug1_coeff) == 4) && (length(p.drug2_coeff)==2)) %drug 2 is monotonic and first drug polynomial of order 3
        y = intercept +p.drug1_coeff(2)*x1 + p.drug1_coeff(3)*x1.^2 + p.drug1_coeff(4)*x1.^3+ ...
            p.drug2_coeff(2)*x2+  p.int_coeff*x1.*x2;
        firstOrder = 2; %first order coefficient for drug 2
        thirdOrder = 1; %third order coefficient for drug 1

    elseif ((length(p.drug1_coeff) == 4) && (length(p.drug2_coeff)==4)) %both are of third order
        y = intercept +p.drug1_coeff(2)*x1 + p.drug1_coeff(3)*x1.^2 + p.drug1_coeff(4)*x1.^3+ ...
            p.drug2_coeff(2)*x2+ p.drug2_coeff(3)*x2.^2 + p.drug2_coeff(4)*x2.^3 + p.int_coeff*x1.*x2;
        thirdOrder = 3; %both drugs are of third order
    end     

    %Find partial derivatives of symbolic function
        fx1 = diff(y,x1);    
        fx2 = diff(y,x2); 
        
        if thirdOrder == 0 %when both are quadratic
            sol_fx1_zero = solve(fx1==0,x1);
            sol_fx2_zero = solve(fx2==0,x2);
    
            fx1_zero = matlabFunction(sol_fx1_zero);
            fx2_zero = matlabFunction(sol_fx2_zero);

        elseif ((thirdOrder == 1) && (firstOrder == 0)) %when 1st drug is of degree 3 and second drug is quadratic
            sol_fx1_zero = solve(fx1==0,x1,"Real",true,"ReturnConditions",true);
            sol_fx2_zero = solve(fx2==0,x2);

            fx1_zero = matlabFunction(sol_fx1_zero.x1(2));
            fx2_zero = matlabFunction(sol_fx2_zero);

        elseif ((thirdOrder == 2) && (firstOrder == 0)) %when 2nd drug is of degree 3 and first drug is quadratic
            sol_fx1_zero = solve(fx1==0,x1);
            sol_fx2_zero = solve(fx2==0,x2,"Real",true,"ReturnConditions",true);
    
            fx1_zero = matlabFunction(sol_fx1_zero);
            fx2_zero = matlabFunction(sol_fx2_zero.x2(2));

        elseif ((thirdOrder == 0) && (firstOrder == 1)) %first drug monotonic and second quadratic
            sol_fx2_zero = solve(fx2==0,x2);

            fx1_zero = [];
            fx2_zero = matlabFunction(sol_fx2_zero);

        elseif ((thirdOrder == 0) && (firstOrder == 2)) %first drug is quadratic and second is monotonic
            sol_fx1_zero = solve(fx1==0,x1);

            fx1_zero = matlabFunction(sol_fx1_zero);
            fx2_zero = [];

        elseif ((thirdOrder == 2) && (firstOrder == 1))
            sol_fx2_zero = solve(fx2==0,x2,"Real",true,"ReturnConditions",true);

            fx1_zero = [];
            fx2_zero = matlabFunction(sol_fx2_zero.x2(2));
            
        elseif ((thirdOrder == 1) && (firstOrder == 2))
            sol_fx1_zero = solve(fx1==0,x1,"Real",true,"ReturnConditions",true);

            fx1_zero = matlabFunction(sol_fx1_zero.x1(2));
            fx2_zero = [];

        elseif thirdOrder == 3
            sol_fx1_zero = solve(fx1==0,x1,"Real",true,"ReturnConditions",true);
            sol_fx2_zero = solve(fx2==0,x2,"Real",true,"ReturnConditions",true);
    
            fx1_zero = matlabFunction(sol_fx1_zero.x1(1));

            if p.monotonic == 2 
                fx2_zero = [];
            else
                fx2_zero = matlabFunction(sol_fx2_zero.x2(2));
            end 
        end 
    
    %Convert symbolic function to matlab function for the surface
    f = matlabFunction(y);

    [X,Y] = meshgrid(p.scaled_conc,p.scaled_conc);

    sol = f(X,Y);

    sol = sol./(max(max(sol)));

    sol(sol<0) = NaN;

    nanNbr = sum(sum(isnan(sol))); %Gives the sum of NaN values

end 