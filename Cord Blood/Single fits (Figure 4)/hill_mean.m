function mu = hill_mean(z, x, behavior)

    switch behavior
        case {'monotonic_increasing','monotonic_decreasing'}
            bottom   = z(1);
            logEC50  = z(2);
            logSlope = z(3);
            top      = z(4);

            slope = exp(logSlope);

            
           mu = bottom + (top - bottom) ./ (1 + 10.^((logEC50 - x) .* slope));
           
        case 'non-monotonic'
            E0       = z(1);
            E1       = z(2);
            C1       = z(3);
            logAbsH1 = z(4);
            C2       = z(5);
            logAbsH2 = z(6);

            h1 = -exp(logAbsH1);   % increasing Hill
            h2 =  exp(logAbsH2);   % decreasing Hill

            H1 = 1 ./ (1 + 10.^((x - C1) .* h1));
            H2 = 1 ./ (1 + 10.^((x - C2) .* h2));

            mu = E1 + (E0 - E1) .* H1 .* H2;

        otherwise
            error('Unknown behavior: %s', behavior);
    end
end
