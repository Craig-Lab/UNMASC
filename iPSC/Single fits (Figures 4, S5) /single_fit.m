function [samples,ci95] = single_fit(x,y)

    y_check = y(1:end);

    [peakVal, peakIdx] = max(y_check);

    dropRight = peakVal - y_check(end);
    riseLeft  = peakVal - y_check(1);

    rangeY = max(y_check) - min(y_check);
    if rangeY < 1e-6
        rangeY = 1e-6;
    end

    tolPeak = 0.05;
    tolStep = 0.10 * rangeY;

    dy = diff(y_check);
    hasUpStep   = any(dy >  tolStep);
    hasDownStep = any(dy < -tolStep);

    if peakIdx > 1 && peakIdx < numel(y_check) && ...
       dropRight > tolPeak && riseLeft > tolPeak && ...
       hasUpStep && hasDownStep

        behavior = 'non-monotonic';

    else
        if y_check(end) >= y_check(1)
            behavior = 'monotonic_increasing';
        else
            behavior = 'monotonic_decreasing';
        end
    end

    if strcmp(behavior,"monotonic_increasing") || strcmp(behavior,"monotonic_decreasing")

        x = log10(x(1:end));
        y = y(1:end);

        sigma0 = log(max(std(y),1e-6));

        if strcmp(behavior,"monotonic_increasing")
            z0 = [min(y), mean([min(x), max(x)]), log(2), max(y), sigma0];
            [samples,ci95] = mcmc_hill(x,y,z0,'monotonic_increasing');
        elseif strcmp(behavior, "monotonic_decreasing")
            z0 = [max(y), mean([min(x), max(x)]), log(2), min(y), sigma0];
            [samples,ci95] = mcmc_hill(x,y,z0,'monotonic_decreasing');
        end

    else
        x = log10(x(1:end));
        y = y(1:end);

        E0_0        = max(y);                  
        E1_0        = min([y(1), y(end)]);    
        logEC501_0  = min(x) + 0.25*(max(x)-min(x));
        logAbsH1_0  = log(2);                  
        logEC502_0  = min(x) + 0.75*(max(x)-min(x));
        logAbsH2_0  = log(2);                  
        logSigma0   = log(max(std(y), 1e-6));

        z0 = [E0_0, E1_0, logEC501_0, logAbsH1_0, ...
              logEC502_0, logAbsH2_0, logSigma0];

        [samples, ci95] = mcmc_hill(x,y,z0,'non-monotonic');
    end
end


