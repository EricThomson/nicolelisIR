%Parameters to set for neuronal analysis, used for all animals--called in loadingSavingPrep

%Following around line 166 in channelPsthPlotScript (line 331 of neurAnal6)
prestimPsthDuration = 0.75; %plot time before onset of microstim for psth (usually around a second)
poststimPsthDuration = prestimPsthDuration;  %time after offset of final microstim for psth (usually around a second)
artifactDuration = 0.085;  %How long you set to stimulate (plus some slush)--you usually stimulate for 75ms 
countDuration = 0.06; %was minimumInterArtifactTime, but variable and always greter than 0.6
binWidth = 0.01;
numSubstimMin = 6;  % for psth includes the numSubstimMin/2 first and last stim (so trials with <numSubstimMin stim are not included)
plotFull = 0;  %plot detailed results, including bits from plotMinimal
plotMinimal = 1; %plots psth, mean stimulus popvec, proportion activated, and receptive fields
printOn = 0;  %print results