function [fit_parameters,resnorm,residual,out,exitflag] = ...
    fit_surface_monteCarlo(x1,x2,p,data,bestFit,ci1,ci2)

X = x1;
Y = x2;

p.alpha12 = 1;
p.alpha21 = 1;

E0_lb = 0;
E0_ub = 1.05;
E1_lb = 0;
E1_ub = 1.05;
E2_lb = 0;
E2_ub = 1.05;
E3_lb = 0;
E3_ub = 1.05;

warning('off','MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:singularMatrix')

warning('off','all')


options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective', ...
    'FunctionTolerance',1e-6,'StepTolerance',1e-6,'OptimalityTolerance', ...
    1e-6,'display','off');

alpha_lb = 10^(-2); 
alpha_ub = 10^2; 

if p.monotonic == 0
    % both drugs non-monotonic
    p.E0 = bestFit(1);
    p.E1 = bestFit(2);
    p.E2 = bestFit(3);
    p.E3 = bestFit(4);

    p.h11 = bestFit(5);   
    p.h12 = bestFit(6);   
    p.h21 = bestFit(7);   
    p.h22 = bestFit(8);   

    p.C11 = bestFit(9);
    p.C12 = bestFit(10);
    p.C21 = bestFit(11);
    p.C22 = bestFit(12);

    p.alpha12 = bestFit(13);
    p.alpha21 = bestFit(14);

    param_guess = [p.E0,p.E1,p.E2,p.E3,p.h11,p.h12,p.h21,p.h22,...
        p.C11,p.C12,p.C21,p.C22,p.alpha12,p.alpha21];

    C11_lb = ci1(1,3); 
    C11_ub = ci1(2,3);
    C12_lb = ci1(1,5); 
    C12_ub = ci1(2,5);

    C21_lb = ci2(1,3);
    C21_ub = ci2(2,3);
    C22_lb = ci2(1,5); 
    C22_ub = ci2(2,5);

    h11_lb = ci1(1,4);
    h11_ub = ci1(2,4);
    h12_lb = ci1(1,6);
    h12_ub = ci1(2,6);

    h21_lb = ci2(1,4);
    h21_ub = ci2(2,4);
    h22_lb = ci2(1,6);
    h22_ub = ci2(2,6);


    lb = [E0_lb E1_lb E2_lb E3_lb ...
        h11_lb h12_lb h21_lb h22_lb ...
        C11_lb C12_lb C21_lb C22_lb  ...
        alpha_lb alpha_lb];
    
    ub = [E0_ub E1_ub E2_ub E3_ub ...
        h11_ub h12_ub h21_ub h22_ub ...
        C11_ub C12_ub C21_ub C22_ub  ...
        alpha_ub alpha_ub];

    [fit_parameters,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@residualsfunction, param_guess, lb, ub, options);

    out = output;

elseif p.monotonic == 1
    % first drug monotonic, second drug non-monotonic
    p.E0 = bestFit(1);
    p.E1 = bestFit(2);
    p.E2 = bestFit(3);
    p.E3 = bestFit(4);

    p.h11 = -bestFit(5);  
    p.h21 = bestFit(6);   
    p.h22 = bestFit(7);   

    p.C11 = bestFit(8);
    p.C21 = bestFit(9);
    p.C22 = bestFit(10);

    p.alpha12 = bestFit(11);
    p.alpha21 = bestFit(12);

    param_guess = [p.E0,p.E1,p.E2,p.E3,p.h11,p.h21,p.h22,...
        p.C11,p.C21,p.C22,p.alpha12,p.alpha21];


    C11_lb = ci1(1,2); 
    C11_ub = ci1(2,2);

    C21_lb = ci2(1,3);
    C21_ub = ci2(2,3);
    C22_lb = ci2(1,5); 
    C22_ub = ci2(2,5);

    h11_lb = -ci1(2,3);
    h11_ub = -ci1(1,3);

    h21_lb = ci2(1,4);
    h21_ub = ci2(2,4);
    h22_lb = ci2(1,6);
    h22_ub = ci2(2,6);

    
    lb = [E0_lb E1_lb E2_lb E3_lb ...
        h11_lb h21_lb h22_lb ...
        C11_lb C21_lb C22_lb  ...
        alpha_lb alpha_lb];
    
    ub = [E0_ub E1_ub E2_ub E3_ub ...
        h11_ub h21_ub h22_ub ...
        C11_ub C21_ub C22_ub  ...
        alpha_ub alpha_ub];


    [fit_parameters,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@residualsfunction1, param_guess, lb, ub, options);

    out = output;

elseif p.monotonic == 2
    % first drug non-monotonic, second drug monotonic
    p.E0 = bestFit(1);
    p.E1 = bestFit(2);
    p.E2 = bestFit(3);
    p.E3 = bestFit(4);

    p.h11 = bestFit(5);   
    p.h12 = bestFit(6);   
    p.h21 = -bestFit(7);   

    p.C11 = bestFit(8);
    p.C12 = bestFit(9);
    p.C21 = bestFit(10);

    p.alpha12 = bestFit(11);
    p.alpha21 = bestFit(12);

    param_guess = [p.E0,p.E1,p.E2,p.E3,p.h11,p.h12,p.h21,...
        p.C11,p.C12,p.C21,p.alpha12,p.alpha21];


    C11_lb = ci1(1,3); 
    C11_ub = ci1(2,3);
    C12_lb = ci1(1,5); 
    C12_ub = ci1(2,5);

    C21_lb = ci2(1,2); 
    C21_ub = ci2(2,2);

    h11_lb = ci1(1,4);
    h11_ub = ci1(2,4);
    h12_lb = ci1(1,6);
    h12_ub = ci1(2,6);

    h21_lb = -ci2(2,3);
    h21_ub = -ci2(1,3);


    lb = [E0_lb E1_lb E2_lb E3_lb ...
        h11_lb h12_lb h21_lb ...
        C11_lb C12_lb C21_lb ...
        alpha_lb alpha_lb];
    
    ub = [E0_ub E1_ub E2_ub E3_ub ...
        h11_ub h12_ub h21_ub ...
        C11_ub C12_ub C21_ub  ...
        alpha_ub alpha_ub];

    [fit_parameters,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@residualsfunction2, param_guess, lb, ub, options);

    out = output;

elseif p.monotonic == 3
    % both drugs monotonic
    p.E0 = bestFit(1);
    p.E1 = bestFit(2);
    p.E2 = bestFit(3);
    p.E3 = bestFit(4);

    p.h11 = -bestFit(5);   
    p.h21 = -bestFit(6);   

    p.C11 = bestFit(7);
    p.C21 = bestFit(8);

    p.alpha12 = bestFit(9);
    p.alpha21 = bestFit(10);

    param_guess = [p.E0,p.E1,p.E2,p.E3,p.h11,p.h21,...
        p.C11,p.C21,p.alpha12,p.alpha21];

    C11_lb = ci1(1,2); 
    C11_ub = ci1(2,2);

    C21_lb = ci2(1,2); 
    C21_ub = ci2(2,2);

    h11_lb = -ci1(2,3);
    h11_ub = -ci1(1,3);

    h21_lb = -ci2(2,3);
    h21_ub = -ci2(1,3);

    lb = [E0_lb E1_lb E2_lb E3_lb ...
        h11_lb h21_lb ...
        C11_lb C21_lb ...
        alpha_lb alpha_lb];
    
    ub = [E0_ub E1_ub E2_ub E3_ub ...
        h11_ub h21_ub ...
        C11_ub C21_ub  ...
        alpha_ub alpha_ub];


    [fit_parameters,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@residualsfunction3, param_guess, lb, ub, options);

    out = output;
else
    error('Unknown p.monotonic value: %d', p.monotonic);
end

    
% -------------------------------------------------------------------------
 function val = residualsfunction(param)
    
    E0 = param(1);
    E1 = param(2);
    E2 = param(3);
    E3 = param(4);

    h11 = param(5);
    h12 = param(6);
    h21 = param(7);
    h22 = param(8);

    C11 = param(9);
    C12 = param(10);
    C21 = param(11);
    C22 = param(12);
    
    alpha12 = param(13);
    alpha21 = param(14);

    E = zeros(length(X),length(Y));

   for l = 1:length(X)
      for m = 1:length(Y)
        d1 = X(l); %drug concentration in ng/mL
        d2 = Y(m); %drug concentration in ng/mL

            r11 = 100;
            r12 = 100;

            r21 = 100;
            r22 = 100;

            rneg11 = r11*10^(C11*h11);
            rneg12 = r12*10^(C12*h12);

            rneg21 = r21*10^(C21*h21);
            rneg22 = r22*10^(C22*h22);
    
            M11 = -(r11*rneg12*d1^h11 + r12*rneg11*d1^h12 + r11*r12*d1^(h11+h12) + ...
                r21*rneg22*d2^h21 + r22*rneg21*d2^h22 + r21*r22*d2^(h21+h22));
    
            M12 = rneg11*rneg12;
    
            M13 = rneg21*rneg22;
    
            M14 = 0;
    
            M21 = r11*rneg12*d1^h11 + r12*rneg11*d1^h12 + r11*r12*d1^(h11+h12);
    
            M22 = -(M12 + r21*rneg22*(alpha12*d2)^h21 + r22*rneg21*(alpha12*d2)^h22 + r21*r22*(alpha12*d2)^(h21+h22));
    
            M23 = 0;
    
            M24 = rneg21*rneg22;
    
            M31 = r21*rneg22*d2^h21 + r22*rneg21*d2^h22 + r21*r22*d2^(h21+h22);
    
            M32 = 0;
    
            M33 =  -(M13 + r11*rneg12*(alpha21*d1)^(h11)+r12*rneg11*(alpha21*d1)^(h12)+ r11*r12*(alpha21*d1)^(h11+h12));
    
            M34 = rneg11*rneg12;
    
            M41 = 1;
    
            M42 = 1;
    
            M43 = 1;
    
            M44 = 1;
    
           M = [M11 M12 M13 M14;
                 M21 M22 M23 M24;
                 M31 M32 M33 M34;
                 M41 M42 M43 M44];

        test = M\([0 0 0 1]');

        if isnan(data(l,m))
            E(l,m) = NaN;
        else 
            E(l,m) = [E0 E1 E2 E3]*test;
        end 

      end 
  end 
    
    val = reshape(data-E,[],1);
    val(isnan(val))=[];

end 

function val = residualsfunction1(param)

    E0 = param(1);
    E1 = param(2);
    E2 = param(3);
    E3 = param(4);

    h11 = param(5);
    h21 = param(6);
    h22 = param(7);

    C11 = param(8);
    C21 = param(9);
    C22 = param(10);
    
    alpha12 = param(11);
    alpha21 = param(12);

    E = zeros(length(X),length(Y));

     for l = 1:length(X)
        for m = 1:length(Y)
            d1 = X(l); %drug concentration in ng/mL
            d2 = Y(m); %drug concentration in ng/mL
    
            r11 = 100;

            r21 = 100;
            r22 = 100;

            rneg11 = r11*10^(C11*h11);

            rneg21 = r21*10^(C21*h21);
            rneg22 = r22*10^(C22*h22);
            
    
            M11 = -(r11*d1^h11  + ...
                r21*rneg22*d2^h21 + r22*rneg21*d2^h22 + r21*r22*d2^(h21+h22));
    
            M12 = rneg11;
    
            M13 = rneg21*rneg22;
    
            M14 = 0;
    
            M21 = r11*d1^h11;
    
            M22 = -(M12 + r21*rneg22*(alpha12*d2)^h21 + r22*rneg21*(alpha12*d2)^h22 + r21*r22*(alpha12*d2)^(h21+h22));
    
            M23 = 0;
    
            M24 = rneg21*rneg22;
    
            M31 = r21*rneg22*d2^h21 + r22*rneg21*d2^h22 + r21*r22*d2^(h21+h22);
    
            M32 = 0;
    
            M33 =  -(M13 + r11*(alpha21*d1)^(h11));
    
            M34 = rneg11;
    
            M41 = 1;
    
            M42 = 1;
    
            M43 = 1;
    
            M44 = 1;
    
            M = [M11 M12 M13 M14;
                 M21 M22 M23 M24;
                 M31 M32 M33 M34;
                 M41 M42 M43 M44];

            test = M\([0 0 0 1]');

            if isnan(data(l,m))
                E(l,m) = NaN;
            else 
                E(l,m) = [E0 E1 E2 E3]*test;
            end 
        end 
     end 
     val = reshape(data-E,[],1);
    val(isnan(val))=[];

end

function val = residualsfunction2(param)
        E0 = param(1);
        E1 = param(2);
        E2 = param(3);
        E3 = param(4);
    
        h11 = param(5);
        h12 = param(6);
        h21 = param(7);
    
        C11 = param(8);
        C12 = param(9);
        C21 = param(10);
        
        alpha12 = param(11);
        alpha21 = param(12);
    
        E = zeros(length(X),length(Y));
    
         for l = 1:length(X)
            for m = 1:length(Y)
                d1 = X(l); %drug concentration in ng/mL
                d2 = Y(m); %drug concentration in ng/mL
        
                r11 = 100;
                r12 = 100;

                r21 = 100;

                rneg11 = r11*10^(C11*h11);
                rneg12 = r12*10^(C12*h12);

                rneg21 = r21*10^(C21*h21);
                
        
                M11 = -(r11*rneg12*d1^h11 + r12*rneg11*d1^h12 + r11*r12*d1^(h11+h12) + ...
                    r21*d2^h21);
        
                M12 = rneg11*rneg12;
        
                M13 = rneg21;
        
                M14 = 0;
        
                M21 = r11*rneg12*d1^h11 + r12*rneg11*d1^h12 + r11*r12*d1^(h11+h12);
        
                M22 = -(M12 + r21*(alpha12*d2)^h21);
        
                M23 = 0;
        
                M24 = rneg21;
        
                M31 = r21*d2^h21;
        
                M32 = 0;
        
                M33 =  -(M13 + r11*rneg12*(alpha21*d1)^(h11)+r12*rneg11*(alpha21*d1)^(h12)+ r11*r12*(alpha21*d1)^(h11+h12));
        
                M34 = rneg11*rneg12;
        
                M41 = 1;
        
                M42 = 1;
        
                M43 = 1;
        
                M44 = 1;
        
                M = [M11 M12 M13 M14;
                     M21 M22 M23 M24;
                     M31 M32 M33 M34;
                     M41 M42 M43 M44];
    
                test = M\([0 0 0 1]');

                if isnan(data(l,m))
                    E(l,m) = NaN;
                else 
                    E(l,m) = [E0 E1 E2 E3]*test;
                end 
            end 
         end

    val = reshape(data-E,[],1);
    val(isnan(val))=[];
end 

function val = residualsfunction3(param)
        E0 = param(1);
        E1 = param(2);
        E2 = param(3);
        E3 = param(4);
    
        h11 = param(5);
        h21 = param(6);

        C11 = param(7);
        C21 = param(8);
        
        alpha12 = param(9);
        alpha21 = param(10);
    
        E = zeros(length(X),length(Y));
    
         for l = 1:length(X)
            for m = 1:length(Y)
                d1 = X(l); %drug concentration in ng/mL
                d2 = Y(m); %drug concentration in ng/mL
        
                r11 = 100;
                r21 = 100;
        
                rneg11 = r11*10^(C11*h11);
                rneg21 = r21*10^(C21*h21);

                M11 = -(r11*d1^h11 + r21*d2^h21);
        
                M12 = rneg11;
        
                M13 = rneg21;
        
                M14 = 0;
        
                M21 = r11*d1^h11;
        
                M22 = -(M12 + r21*(alpha12*d2)^h21);
        
                M23 = 0;
        
                M24 = rneg21;
        
                M31 = r21*d2^h21;
        
                M32 = 0;
        
                M33 =  -(M13 + r11*(alpha21*d1)^(h11));
        
                M34 = rneg11;
        
                M41 = 1;
        
                M42 = 1;
        
                M43 = 1;
        
                M44 = 1;
        
                M = [M11 M12 M13 M14;
                     M21 M22 M23 M24;
                     M31 M32 M33 M34;
                     M41 M42 M43 M44];
    
                test = M\([0 0 0 1]');

                if isnan(data(l,m))
                    E(l,m) = NaN;
                else 
                    E(l,m) = [E0 E1 E2 E3]*test;
                end 
            end 
         end 
    val = reshape(data-E,[],1);
    val(isnan(val))=[];
end 
end 