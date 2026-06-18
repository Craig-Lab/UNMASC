function [samples,ci95] = mcmc_hill(x, y, z0, behavior)

    logpdf = @(z) logposterior_hill(z, x, y, behavior);

    nsamples = 20000;
    burnin   = 5000;
    thin     = 1;

    raw = slicesample(z0(:)', nsamples, ...
        "logpdf", logpdf, ...
        "burnin", burnin, ...
        "thin", thin);

    switch behavior
        case {'monotonic_increasing','monotonic_decreasing'}
            bottom   = raw(:,1);
            logEC50  = raw(:,2);
            slope    = exp(raw(:,3));
            top      = raw(:,4);
            sigma    = exp(raw(:,5));

            samples = [bottom, logEC50, slope, top, sigma];

        case 'non-monotonic'
            E0       = raw(:,1);
            E1       = raw(:,2);
            logEC501 = raw(:,3);
            h1       = -exp(raw(:,4));   % enforce negative slope
            logEC502 = raw(:,5);
            h2       =  exp(raw(:,6));   % enforce positive slope
            sigma    = exp(raw(:,7));

            samples = [E0, E1, logEC501, h1, logEC502, h2, sigma];

        otherwise
            error('Unknown behavior: %s', behavior);
    end

    ci95 = prctile(samples, [2.5 97.5]);
end
