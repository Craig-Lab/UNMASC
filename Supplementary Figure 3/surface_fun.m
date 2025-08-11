% Function to evaluate surface at data points of interest
% INPUTS: 
    % drug1 and drug2: vectors of cytokine concentrations (in ng/mL)
    % param: vector of parameters for the fit
    % p : 'struct' object containing 
    % monotonic: 
            % if monotonic == 0: drugs 1 and 2 are non-monotonic
            % if monotonic == 1: drug 1 is monotonic and drug 2 is
            % non-monotonic
            % if monotonic == 2: drug 1 is non-monotonic and drug 2 is
            % monotonic
            % if monotonic ==3: drugs 1 and 2 are monotonic
%OUTPUTS:
    % E: matrix containing the points resulting from the evaluation of the
    % surface at each concentration pair


function E = surface_fun(drug1,drug2,param,p)

warning('off','MATLAB:nearlySingularMatrix')

if p.monotonic == 0

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
    
        E = zeros(length(drug1),length(drug2));
    
         for l = 1:length(drug1)
            for m = 1:length(drug2)
                d1 = drug1(l); %drug concentration in ng/mL
                d2 = drug2(m); %drug concentration in ng/mL
        
                r11 = 100;
                r12 = 100;

                r21 = 100;
                r22 = 100;

                rneg11 = r11*10^(C11*h11);
                rneg12 = r12*10^(C12*h12);

                rneg21 = r21*10^(C21*h21);
                rneg22 = r22*10^(C22*h22);
        
                M11 = -(r11*rneg12*d1^(h11) + r12*rneg11*d1^(h12) + r11*r12*d1^(h11+h12) + ...
                    r21*rneg22*d2^(h21) + r22*rneg21*d2^(h22) + r21*r22*d2^(h21+h22));

                M12 = rneg11*rneg12;
        
                M13 = rneg21*rneg22;
        
                M14 = 0;
        
                M21 = r11*rneg12*d1^(h11) + r12*rneg11*d1^(h12) + r11*r12*d1^(h11+h12);
        
                M22 = -(M12 + r21*rneg22*(alpha12*d2)^(h21) + r22*rneg21*(alpha12*d2)^(h22) + r21*r22*(alpha12*d2)^(h21+h22));
        
                M23 = 0;
        
                M24 = rneg21*rneg22;
        
                M31 = r21*rneg22*d2^(h21) + r22*rneg21*d2^(h22) + r21*r22*d2^(h21+h22);
        
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
                
                 E(m,l) = [E0 E1 E2 E3]*test;
            end 
         end 

elseif p.monotonic == 1
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
    
        E = zeros(length(drug1),length(drug2));
    
         for l = 1:length(drug1)
            for m = 1:length(drug2)
                d1 = drug1(l); %drug concentration in ng/mL
                d2 = drug2(m); %drug concentration in ng/mL
        
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
                
                 E(m,l) = [E0 E1 E2 E3]*test;
            end 
         end 

elseif p.monotonic == 2

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
    
        E = zeros(length(drug1),length(drug2));
    
         for l = 1:length(drug1)
            for m = 1:length(drug2)
                d1 = drug1(l); %drug concentration in ng/mL
                d2 = drug2(m); %drug concentration in ng/mL
        
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
                
                 E(m,l) = [E0 E1 E2 E3]*test;
            end 
         end 

elseif p.monotonic == 3 %both drugs are monotonic (reduces to MuSyC)

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
    
        E = zeros(length(drug1),length(drug2));
    
         for l = 1:length(drug1)
            for m = 1:length(drug2)
                d1 = drug1(l); %drug concentration in ng/mL
                d2 = drug2(m); %drug concentration in ng/mL
        
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
                
                E(m,l) = [E0 E1 E2 E3]*test;
            end 
         end 

end 
end 