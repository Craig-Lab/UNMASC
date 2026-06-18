
function lp = logposterior_hill(z, x, y, behavior)

    mu = hill_mean(z, x, behavior);

    if any(~isfinite(mu))
        lp = -Inf;
        return
    end

    % Log-normal prior parameters 
    mu_logSlope = log(1.5);
    sd_logSlope = 0.3;

    switch behavior
        case {'monotonic_increasing','monotonic_decreasing'}
            bottom   = z(1);
            logEC50  = z(2);
            logSlope = z(3);
            top      = z(4);
            logSigma = z(5);

            xmin = min(x);
            xmax = max(x);
            
            if logEC50 < xmin || logEC50 > xmax
                lp = -Inf;
                return
            end

            slope = exp(logSlope);
            sigma = exp(logSigma);

            if ~isfinite(slope) || ~isfinite(sigma) || sigma <= 0
                lp = -Inf;
                return
            end
            
            if slope > 5 
                lp = -Inf;
                return
            end

            if bottom < 0 || bottom > 1.1 || top < 0 || top > 1.1
                lp = -Inf;
                return
            end

            if (strcmp(behavior,"monotonic_increasing") && bottom > top)
                lp = -Inf;
                return
            elseif (strcmp(behavior,"monotonic_decreasing") && bottom < top)
                lp = -Inf;
                return
            end

            lower = min(bottom, top);
            upper = max(bottom, top);

            lp_top = -0.5*((upper - max(y))/0.001)^2;
            lp_bottom = -0.5*((lower - min(y))/(0.01))^2;

            resid = y - mu;
            ll = sum(-0.5*log(2*pi) - log(sigma) - 0.5*(resid./sigma).^2);
            lp_logEC50  = -0.5*((logEC50 - median(x))/1)^2;
            lp_logSigma = -0.5*(logSigma/2)^2;

            lp_slope = -0.5*((logSlope - mu_logSlope)/sd_logSlope)^2;

            lp = ll + lp_logEC50 + lp_slope + lp_bottom + lp_top + lp_logSigma;

        case 'non-monotonic'
            E0         = z(1);
            E1         = z(2);
            logEC501   = z(3);
            logAbsH1   = z(4);
            logEC502   = z(5);
            logAbsH2   = z(6);
            logSigma   = z(7);

            xmin = min(x);
            xmax = max(x);
            
            if logEC501 < xmin || logEC501 > xmax || ...
               logEC502 < xmin || logEC502 > xmax
                lp = -Inf;
                return
            end

            absH1 = exp(logAbsH1);
            absH2 = exp(logAbsH2);
            h1    = -absH1;
            h2    =  absH2;
            sigma = exp(logSigma);

            if ~isfinite(absH1) || ~isfinite(absH2) || ~isfinite(sigma) || sigma <= 0
                lp = -Inf;
                return
            end

            if absH1 > 5 || absH2 > 5
                lp = -Inf;
                return
            end

            if E0 < 0 || E0 > 1.1 || E1 < 0 || E1 > 1.1
                lp = -Inf;
                return
            end

            if logEC501 >= logEC502
                lp = -Inf;
                return
            end

            if E0 < E1
                lp = -Inf;
                return
            end

            resid = y - mu;
            ll = sum(-0.5*log(2*pi) - log(sigma) - 0.5*(resid./sigma).^2);

            sy = max(std(y), 1e-3);


            lp_E0 = -0.5*((E0 - 1)/0.01)^2;
            lp_E1       = -0.5*((E1 - min([y(1), y(end)]))/(2*sy + 1e-6))^2;
            lp_logEC501 = -0.5*((logEC501 - median(x))/1)^2;
            lp_logEC502 = -0.5*((logEC502 - median(x))/1)^2;
            lp_logSigma = -0.5*(logSigma/2)^2;

          
            lp_h1 = -0.5*((logAbsH1 - mu_logSlope)/sd_logSlope)^2;
            lp_h2 = -0.5*((logAbsH2 - mu_logSlope)/sd_logSlope)^2;

            lp = ll + lp_E0 + lp_E1 + lp_logEC501 + lp_logEC502 ...
                   + lp_h1 + lp_h2 + lp_logSigma;

        otherwise
            error('Unknown behavior: %s', behavior);
    end
end


