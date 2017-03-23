function frequencyViewer(stimSequences, stimTimeSequences, eventTimes, eventNames, trialFirstClockTimes, timeWindow, colors)
%Note finished documenting:
%frequencyViewer(stimSequences, stimTimeSequences, eventTimes, eventNames, trialFirstClockTimes, timeWindow, colors)
%
%Check to see that binary encoding/decoding is working properly by plotting reference (clock) channel and
%one of the other channels with binary as well as number decoded superimposed.
%
%Plot reference clock events, channel binary encoded events, first reference event indicated in green.
%For each substimulus, show decoded stimulus (e.g., 010 is 2) superimposed on the stimulus pulses).

numTrials = length(stimSequences);
%Pick trials to plot
trialsBelowWindow = cellfun(@(x) sum(find(x < timeWindow(1))), stimTimeSequences);
if sum(trialsBelowWindow)
    firstTrial = find(trialsBelowWindow, 1, 'last')+1;
else
    firstTrial = 1;
end
trialsAboveWindow = cellfun(@(x) sum(find(x > timeWindow(2))), stimTimeSequences);
if sum(trialsAboveWindow)
    lastTrial = find(trialsAboveWindow, 1, 'first')-1;
else
    lastTrial = numTrials;
end

trialsToPlot = [firstTrial: lastTrial];

%Overlay stimulus calculated with event view to see if they match
eventViewer(eventTimes, eventNames, timeWindow, {'m','b'}); hold on
for trialNum = trialsToPlot
    tmpTimes = stimTimeSequences{trialNum}; %onset time if none
    tmpStim = stimSequences{trialNum};  %nan if none
    %scatter(tmpTimes, zeros(1, length(tmpTimes)), 15, 'filled', 'CData', rand(1,3));  %strip at bottom of same color
    if isnan(tmpStim)
        text(tmpTimes, 1, num2str('nada'));
    else
        for stimNum = 1: length(tmpStim)
            text(tmpTimes(stimNum), 1, num2str(tmpStim(stimNum)), 'FontSize', 10);
        end
    end
end

scatter(trialFirstClockTimes(trialsToPlot), ones(1,length(trialsToPlot)), 20, 'g', 'filled');
title(['Stim event pulses with rate overlaid. Channel ' eventNames{2}]);
set(gca, 'XLim', timeWindow);
shg

