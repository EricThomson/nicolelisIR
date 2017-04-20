%channelPsthPlotScript.m
%
% Calculate and plot PSTH for a single channel, blanked during artifact period
% Usually you cycle through this for multiple channels (typically called in neuronalAnalysisN.m)
%
%Note for baseline counts and psth:
%    For graphical representation of the baseline bit see baselineTweaking.ai.
%    You are basically trying to go clearly into a time between trials, not right before the trial starts (i.e.,
%       not right before the artifact begins, partly because there are sometimes artifacts there, but also because
%       you are trying to get a kind of task-free baseline where they're sort of doing nothing.

      
%% Remove "spikes" (event times) within artifact windows
neur_suffix={'sig'};
%raw psth, including artifact
eventTimes = eval(channelName);

%Pull psth, with data blanked within stim artifact windows
artifactBuffer = 0; %time after stimulation offset to keep blank
numArtifacts = length(artifactOnTimes);  %this is just totalNumSubstim (see neuronalAnalysisN.m)
blankEventInds = [];
for artifactNum = 1 : numArtifacts
    artifactWindow = [artifactOnTimes(artifactNum) artifactOffTimes(artifactNum)];
    artifactInds = find(eventTimes >= artifactWindow(1) & eventTimes <= artifactWindow(2) + artifactBuffer);
    blankEventInds = [blankEventInds ; artifactInds]; %indices to exclude data b/c of stim artifact
end
eventTimesBlanked = eventTimes;
eventTimesBlanked(blankEventInds) = [];


%% Extract spike counts before/after stimulation for same neuron 
% Probably switch to just this after all it's just the previous, but better
%Pull spike count in countDuration time window after stim artifact offset (with offset buffer)
% Probably just set countDuration to 60ms so it is the same for all animals/sessions
%
%responseVectors{trialNum}: 1 x numSubstim array of spike counts
%substimCounts(substimNum)  1 x numTotalSubstim array of counts 
%responseTimes{trialNum}:   1x numSubstim array of onset of count window
%baselineCounts(trialNum):  spike count in window before first stimulus starts
%baselineTimes(trialNum)    onset of baseline period
%
%
%Recall you previously calculated:
% artifactWindows{trialNum} = numSubstim x 2 window of artifact times
% artifactOnTimes(substimNum)  1 x totalNumSubstim includes all onset times
% artifactOffTimes(substimNum) 1 x totalNumSubstim includes all offset times
% artifactDurations(substimNum) 1 x totalNum
responseVectors = cell(1,numTrials);  %each is 1 x numSubstim array of spike counts in countDuration window
responseTimes = cell(1,numTrials);  %1x numSubstim array of onset of count window
baselineCounts = zeros(1,numTrials); %spike count in window before first stimulus starts
baselineTimes = zeros(1,numTrials); %onset of baseline period

for trialNum = 1:numTrials
    %if rem(trialNum, 50) == 0
    %    display(['Count analysis trial ' num2str(trialNum) ' out of ' num2str(numTrials)]);
    %end
    artifactWindowsTmp = artifactWindows{trialNum}; %all artifacts in trial (numSubstim x 2)
    artifactOffsetsTmp = artifactWindowsTmp(:,2);
    artifactOnsetTmp = artifactWindowsTmp(1,1); %time of first artifact in trial
    numSubstim = numSubstimAll(trialNum);
    responseVecTmp = zeros(1, numSubstim);
    responseTimesTmp = zeros(1,numSubstim);
    baseWindow = [artifactOnsetTmp-countDuration-prestimPsthDuration/2 artifactOnsetTmp-prestimPsthDuration/2];
    %was >= and <= for some strange reason, frankly not sure (maybe b/c it didn't matter because we
    %are not dealing with contiguous time chunks because of artifact gaps)
    baseSpikeInds = find(eventTimes > baseWindow(1) & eventTimes <= baseWindow(2));
    baselineCounts(trialNum) = length(baseSpikeInds);
    baselineTimes(trialNum) = baseWindow(1); % + countDuration/5;
    for substimNum = 1: numSubstim
        artifactOffTime = artifactOffsetsTmp(substimNum);
        responseTimesTmp(substimNum) = artifactOffTime; %artifact off, start counting
        responseWindow = [artifactOffTime artifactOffTime+countDuration];
        responseSpikeInds = find(eventTimes > responseWindow(1) & eventTimes <= responseWindow(2));
        responseVecTmp(substimNum) = length(responseSpikeInds); 
    end
    responseVectors{trialNum} = responseVecTmp;  %1xnumSubstim spike counts for trial
    responseTimes{trialNum} = responseTimesTmp; %1xnumSubstim time of onset of spike counts
end

substimCounts = cell2mat(responseVectors);

%% plot count superimposed on blanked-psth
if plotIndividualChannels  
    timeSubwindow = [100 500];  %
    figure;
    [psth, bin_edges]=session_psth_extract(eventTimesBlanked, timeSubwindow, binWidth, 1);hold on

    maxPSTHTmp = max(psth);
    onArtifacts = artifactOnTimes(find(artifactOnTimes > timeSubwindow(1) & ...
                                       artifactOnTimes < timeSubwindow(2)));
    offArtifacts = artifactOffTimes(find(artifactOffTimes > timeSubwindow(1) & ...
                                       artifactOffTimes < timeSubwindow(2)));
    maxVal = max(psth);                            
    line([onArtifacts'; onArtifacts'], [zeros(1,length(onArtifacts)); ...
        maxVal*ones(1,length(onArtifacts))], 'Color', 'g', 'LineWidth', 0.5, 'LineStyle', '-'); 
    line([offArtifacts'; offArtifacts'], [zeros(1,length(offArtifacts)); ...
        maxVal*ones(1,length(offArtifacts))], 'Color', 'm', 'LineWidth', 0.5, 'LineStyle', '-'); 
    title('Artifact removed. Artifact onset: green, artifact offset: magenta');
    
    windowTrialInds = find(trialOnsetTimes > timeSubwindow(1) & trialOnsetTimes < timeSubwindow(2));
    numWindowTrials = length(windowTrialInds);
    for trialNum = 1 : numWindowTrials
        trialIndex = windowTrialInds(trialNum);
        baselineTime = baselineTimes(trialIndex);
        baselineCount = baselineCounts(trialIndex);
        text(baselineTime, maxVal-.2, num2str(baselineCount), 'FontSize',8);

        tmpResponses = responseVectors{trialIndex};
        tmpTimes = responseTimes{trialIndex};
        numSubstim = length(tmpTimes);
        for substim = 1:numSubstim
            text(tmpTimes(substim), maxVal-0.2, num2str(tmpResponses(substim)), 'FontSize',8);
        end
    end

end



%% Analyze first and last substimNum/2 substimuli. First, pull event times, and mean stimulus population vectors.
%Pick out trials where at least 6 substim.
if plotIndividualChannels
    hist(numSubstimAll, [0: max(numSubstimAll)]);axis tight
    title('Number of substimuli')
end

%Set up 3d array for easy averaging of stimuli
trialsWithMinimalSubstimInds = find(numSubstimAll >= numSubstimMin);
numTrialsWithMinimalSubstim = length(trialsWithMinimalSubstimInds);
nStimFreqs = zeros(numTrialsWithMinimalSubstim, 4, numSubstimMin); %trial x channel x substimNum (1,2, ... end-1, end)

%Set up relevant artifact onset and offset times (includes baseline) (numTrials x 7)
trialSubsetEventTimes = zeros(numTrialsWithMinimalSubstim, numSubstimMin+1); %base, first, second, penult, last

for trialNum = 1: numTrialsWithMinimalSubstim
    trialInd = trialsWithMinimalSubstimInds(trialNum);
    %frequency vector
    stimFreqsTmp = stimvecSequencesFreq{trialInd};    
    %Event times
    artifactWindowsTmp = artifactWindows{trialInd};
    %#? why the -prestimPsthDuration and other stuff?
    
    %baseline time
    trialSubsetEventTimes(trialNum, 1) = baselineTimes(trialInd) - prestimPsthDuration - 0.1 - rand*0.1;
%     For graphical representation of the baseline bit see baselineTweaking.ai.
%     You are basically trying to go clearly into a time between trials, not right before the trial starts (i.e.,
%       not right before the artifact begins, partly because there are sometimes artifacts there, but also because
%       you are trying to get a kind of task-free baseline where they're sort of doing nothing.
    %Recall you got baselineTimes above:
    %   baselineTimes(trialNum) = baseWindow(1);
    %   baseWindow = [artifactOnsetTmp-countDuration-prestimPsthDuration/2 artifactOnsetTmp-prestimPsthDuration/2]; 
    %It is used to calculate baseline PSTH in the cell after this one. Here is the code it will use:
    %   tmpEventTimes = trialSubsetEventTimes(trialNum,:);
    %   baseTmp = histc(eventTimes, [tmpEventTimes(1): binWidth: tmpEventTimes(1) + prestimPsthDuration])';

    
    for substimEpoch = 1: numSubstimMin
        if substimEpoch <= numSubstimMin/2
            %Plug in the first n/2 frequencies
            trialSubsetEventTimes(trialNum, substimEpoch+1) = artifactWindowsTmp(substimEpoch,2);
            nStimFreqs(trialNum, :, substimEpoch) = stimFreqsTmp(substimEpoch,:);
        else %plug in the last n/2 frequencies
            indOffset = numSubstimMin - substimEpoch;
            trialSubsetEventTimes(trialNum, substimEpoch+1) = artifactWindowsTmp(end-indOffset,2);
            nStimFreqs(trialNum, :, substimEpoch) = stimFreqsTmp(end - indOffset,:);
        end
    end
end
%Get mean and sem 
detectorCoords = [1 -1 -1 1;1 1 -1 -1];  %
%popvecTmp(subStim,:) =  detectorCoords * stimvecFreq(subStim, :)'; 
      
mnNStimFreqs = zeros(numSubstimMin, 4);
mnNPopvecs = zeros(numSubstimMin, 2);
for substimEpoch = 1: numSubstimMin
    mnNStimFreqs(substimEpoch,:) = mean(nStimFreqs(:, :, substimEpoch));  %mean frequency
    mnNPopvecs(substimEpoch,:) = detectorCoords * mnNStimFreqs(substimEpoch,:)';  %stimulus population vector
end


%% Get trial psths (prestPsth, numSubstimMin psths, and poststim psth) 
%allPsthTrials contains cell array with psth numTrialsMinSubstim x (numSubstimMin+1)
allPsthTrials = cell(1, numSubstimMin+1);  %7 (baseline and up to response to substimN)
allPsthTrials{1} = zeros(numTrialsWithMinimalSubstim, length([0:binWidth: prestimPsthDuration]));  %before first substim ('baseline')
allPsthTrials{end} = zeros(numTrialsWithMinimalSubstim, length([0:binWidth: poststimPsthDuration])); %after last substimulus
for psthEpoch = 1: numSubstimMin-1
    allPsthTrials{psthEpoch+1} = zeros(numTrialsWithMinimalSubstim, length([0:binWidth: 0.06]));
end
for trialNum = 1: numTrialsWithMinimalSubstim
    tmpEventTimes = trialSubsetEventTimes(trialNum,:);
    %baseline
    baseTmp = histc(eventTimes, [tmpEventTimes(1):binWidth: tmpEventTimes(1) + prestimPsthDuration])';
    if ~isempty(baseTmp)
        allPsthTrials{1}(trialNum,:) = baseTmp;
    end
    %final, after last substim is off
    lastTmp = histc(eventTimes, [tmpEventTimes(end): binWidth: tmpEventTimes(end) + poststimPsthDuration])'; 
    if ~isempty(lastTmp)
        allPsthTrials{end}(trialNum,:) = lastTmp;
    end
    %middle guys (after each substim)
    for substimEpoch = 1: numSubstimMin - 1
        tmpHist = histc(eventTimes, [tmpEventTimes(substimEpoch+1): binWidth: tmpEventTimes(substimEpoch+1) + 0.06])';
        if ~isempty(tmpHist)
            allPsthTrials{substimEpoch+1}(trialNum,:) = tmpHist;
        end
    end     
end

%% Concatenate mean psth and full z-score psth
%Taking mean of the full trial psths surrounding substim etc into a full psth for display and such
numZeroBinsBetween = 20;  %huh?
mnBaseline = mean(mean(allPsthTrials{1}(:,1:end-1))); %end-1 b/c last count is final edge, which = 0
stdBaseline = std(mean(allPsthTrials{1}(:,1:end-1)));
fullMeanPsth = [];
fullMeanPsthZ = [];
for substimEpoch = 1: numSubstimMin + 1
    if substimEpoch <= numSubstimMin
        tmpPsth = [mean(allPsthTrials{substimEpoch}(:,1:end-1)) zeros(1, 7)]; %end-1 b/c last count is final edge = 0
        tmpPsthZ =  [(mean(allPsthTrials{substimEpoch}(:,1:end-1))-mnBaseline)/stdBaseline zeros(1, 7)];
    else
        tmpPsth = [mean(allPsthTrials{substimEpoch}(:,1:end-1))];
        tmpPsthZ =  [(mean(allPsthTrials{substimEpoch}(:,1:end-1))-mnBaseline)/stdBaseline];
    end
    
    if substimEpoch == numSubstimMin/2 + 1
        tmpPsth = [tmpPsth zeros(1,numZeroBinsBetween)];
    end
    
    fullMeanPsth = [fullMeanPsth tmpPsth];
    fullMeanPsthZ = [fullMeanPsthZ tmpPsthZ];    

end



