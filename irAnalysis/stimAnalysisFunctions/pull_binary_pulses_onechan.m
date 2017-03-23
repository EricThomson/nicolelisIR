function [stimSequences, stimTimeSequences] = ...
    pull_binary_pulses_onechan(channelEventTimes, clockTimes, trialOnsetTimes, trialOffsetTimes, jitter)

% [stimSequences, stimTimeSequences] = ...
%     pull_binary_pulses_onechan(channelEventTimes, clockTimes, trialOnsetTimes, trialOffsetTimes, jitter)
%
%Pulls out the decimal representation of stim sequences on each trial (e.g., 011 --> 3). 
%
%Inputs:
%channelEventTimes:
%   Times of all events in the channel (numerical array)
%clockTimes:
%   Time of reference event (the stim clock) that gives a tick on every possible output.
%trialOnsetTimes/trialOffsetTimes:
%   Time of each trial onset/offset (time to start/stop looking for clocktimes): numerical array
%jitter:
%   Number of time window around a clock tick event to look for events in channel.
%
%Outputs
%stimSequences: 
%   1 x numTrials cell array, each contains numeric array of stimulus value 
%   for that trial. NaN if none
%stimTimeSequences: 
%   1xnumTrials cell array, contains numeric array of onset times for each 
%   substimulus. Trial onset time if none.

numTrials = length(trialOnsetTimes);

stimSequences = cell(1,numTrials); %each contains set of stimuli for that channel
stimTimeSequences = cell(1,numTrials);
for trialNum = 1 : numTrials
%     if trialNum == 6
%        'hi'
%     end
    startTrial = trialOnsetTimes(trialNum);
    endTrial = trialOffsetTimes(trialNum);
    clockPulseInds = find(clockTimes >= startTrial & clockTimes <= endTrial);
    numClockPulses = length(clockPulseInds);
    if rem(numClockPulses, 3) == 0
        numStim = numClockPulses/3;
    else
        error(['Clock pulses not multiple of 3. Trial num ' num2str(trialNum)])
    end
    clockPulseTimes = clockTimes(clockPulseInds);
    stimOnsetTimes = clockPulseTimes(1:3:end);
    if length(stimOnsetTimes) ~= numStim
        error(['Stim onset times not matching up on trial ' num2str(trialNum)]);
    end
    
    stimSequence = [];
    for stimNum = 1:numStim
        clockStimInds = 3*(stimNum-1)+1: 3*(stimNum-1)+3;
        clockTimesStim = clockPulseTimes(clockStimInds);
        binResp='000';
        for digitNum = 1:3
            indOn = find(channelEventTimes >= (clockTimesStim(digitNum) - jitter) & ...
                    channelEventTimes <= (clockTimesStim(digitNum) + jitter));
            if indOn
                binResp(digitNum)='1';
            end
        end %checking each digit
        stimSequence(end+1) = bin2dec(binResp);
    end %1:numStim
    if length(stimSequence) ~= length(stimOnsetTimes)
        error(['Length of stimsequence and onsettimes not matched for trial ' num2str(trialNum)]);
    end
    if isempty(stimSequence)
        stimSequences{trialNum} = NaN;   
        stimTimeSequences{trialNum} = startTrial;
    else
        stimSequences{trialNum} = stimSequence;  
        stimTimeSequences{trialNum} = stimOnsetTimes;
    end
end


