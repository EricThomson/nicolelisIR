function plot_peakZ(mnPeakZ,medianPeakZ,peakZScores,numChannels,animalName,dateStr)
%plots peak Z scores
% mnPeakZ = mean of Z scores
% medianPeakZ = median of Z scores
% peakZScores = array of peak z scores
% numChannels = total number of channels
% animalName = name of animal
% dateStr = date of session
% Nicolelis lab 4/7/17

    figure;
    plot([0 numChannels+1],[mnPeakZ mnPeakZ], 'r');hold on
    plot([0 numChannels+1],[medianPeakZ medianPeakZ], 'b')
    legend({'Mean', 'Median'}, 'Location', 'NorthWest')
    scatter([1: numChannels], sort(peakZScores,'ascend'),'k','filled');hold on
    axis([0 numChannels+1 -1 max(peakZScores)+0.1*max(peakZScores)])
    plot([0 numChannels+1],[0 0], 'k')
    title([animalName ' ' strrep(dateStr, '_', '/') ' peak z score (from filtered RF)'])
 