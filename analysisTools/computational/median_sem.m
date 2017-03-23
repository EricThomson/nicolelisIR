function [med, se_median] = median_sem(data, numSamples)
%[med, se_median]=median_sem(data, bootstraps)
%
%Calculates median and standard error of the median using the bootstrap
%for the data vector data. Uses numSamples bootstrap samples.
%
%Based on method in 'An introduction to the bootstrap' by Efron and Tibshirani

try
    med = median(data);
    data = data(~isnan(data));
    if length(data) > 1
        med_samples = bootstrp(numSamples, @median, data);  %bootstrp(numSamples, @median, data');
        se_median = std(med_samples);
    else
        se_median = 0;
    end
catch
    'debug median sem'
end


