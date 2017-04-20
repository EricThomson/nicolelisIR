%Analysis stimuli (microstimulation patterns) for task
% Note there is a wrinkle as the integration task has two ir indicators per trial, as
% two IR lights turn on each trial (distractor and target)! You use the IR indicators to pull
% the number of trials, so herein use protocol name to individuate analysis.
%**This is probably not the greatest strategy.
%
%To do:
%   Try running through with multiple plot options set.
%   Fix plotting bounds to show max/min rather than numbers

%% Load plex data (regarding stim, at least) and session_info from computer
%datFile is from loadingSavingPrep
cd(datFolder)

eval(['load ' animalName '_' dateStr]) %data from Matlab session info
sessDat = session_info; 
clear session_info;
numTrialsSess = sessDat.data.numTrials;
pcSess = sessDat.data.percent_correct;

%analysis of data
load stimEvents  %previously saved from plexon, saved events 11 (if applicable), and 12-16


%Get stim clock, stim events, ir onset data
%note events 7-10 are visible light on
stimClock = Event012;
chanA = Event013;
chanB = Event014;
chanC = Event015;
chanD = Event016;
try
    ir1On = Event003;
catch
    error('Need event003')
end
ir2On = Event004;
ir3On = Event005;
ir4On = Event006;

onIR = sort([ir1On' ir2On' ir3On' ir4On']);
%integration task has two ir indicators per trial
if isequal(sessDat.protocol, 'irVisualIntegrationTask')
    onIR= onIR(1:2:end);
end
    

%Compare N from sessdat, and N from plexon events
numTrialsEvents = length(onIR);
if numTrialsEvents ~= numTrialsSess
    warning('Num trials from irOn events is not same as from session');
    numTrials = numTrialsSess;
else
    display('N is the same using irOn events and sessDat...Bullseye!');
    numTrials = numTrialsSess;
end



%View raw events in event viewer
allEvents ={onIR, stimClock, chanA, chanB, chanC, chanD};
stimEvents = {chanA, chanB, chanC, chanD};

chanVals = {'A','B','C','D'};
eventColors ={'m', 'b','g','r', 'k', 'c'};
eventNames = {'onIR', 'stimClock', 'Chan A', 'Chan B', 'Chan C', 'Chan D'};
timeWindow = [Start Stop];

if plotFull
    figure;
    eventViewer(allEvents, eventNames, timeWindow, eventColors);hold on
    title('Stim, clock, and irOn events')
end

%% Inspect time between trials 
interTrialTimes = diff(onIR);
interTrialMin = min(interTrialTimes);
interTrialMax = max(interTrialTimes);
[mnITI, semITI] = mean_sem(interTrialTimes');
display(['Min intertrial interval: ' num2str(interTrialMin) ' . mean/sem: ' num2str(mnITI) ' ' num2str(semITI)]);

%Plot time between irOn events
if plotFull
    figure
    scatter(onIR(2:end)/60,interTrialTimes);hold on
    axis([-1 max(onIR)/60+2 -.05 max(interTrialTimes)+1]);
    xlabel('Trial onset time (s)'); title('Time(trial_n)-Time(trial_{n-1}) ');
    ylabel('\Delta t')
    grid on;shg
end

%% Get time of first and last stimClock on each trial (trialFirstClockTimes, trialLastClockTimes)
%   Note on some trials there may be none, for instance if the button was locked on
%   so it was an immediate error. For these, use trial onset time, onset time + .05.
trialOnsetTimes = onIR(1:numTrials);

trialFirstClockTimes = zeros(1,numTrials);
trialLastClockTimes = zeros(1,numTrials);
for trialNum = 1:numTrials
    trialOnsetTime = trialOnsetTimes(trialNum);
    if trialNum < numTrials %if not last trial, then use next trial onset as upper bound
        nextTrialTime = trialOnsetTimes(trialNum + 1);
        trialClockInds = find(stimClock > trialOnsetTime & stimClock < nextTrialTime);
        if ~isempty(trialClockInds) %clock events during trial
            trialFirstClockTimes(trialNum) = stimClock(trialClockInds(1));
            trialLastClockTimes(trialNum) = stimClock(trialClockInds(end));
        else %case of no clock events (Jasmin 11/15 around 14 minutes in)
            trialFirstClockTimes(trialNum) = trialOnsetTime;
            trialLastClockTimes(trialNum) = trialOnsetTime + 0.05;  %this is arbitrary
        end
    else %for final trial, relevant clock tics are simply those greater than trial onset
        trialClockInds = find(stimClock > trialOnsetTime);
        if ~isempty(trialClockInds)
            trialFirstClockTimes(trialNum) = stimClock(trialClockInds(1));
            trialLastClockTimes(trialNum) = stimClock(end);
        else
            trialFirstClockTimes(trialNum) = trialOnsetTime;
            trialLastClockTimes(trialNum) = trialOnsetTime + 0.05;
        end
    end
end 

%Check to see if trialFirstClockTimes and trialLastClockTimes seem reasonable
if plotFull
    figure;
    eventViewer({stimClock}, {'Clock'}, timeWindow, {'b'}); hold on
    scatter(trialFirstClockTimes, 0.5*ones(1,numTrials), 20,'filled','g');
    scatter(trialLastClockTimes,  0.5*ones(1,numTrials), 20,'filled','m'); 
    title('Green is clock onset, magenta offset, for a given trial')
end
    

%% percent correct analysis, and compare with sessDat
%integration task has two ir indicators per trial
if isequal(sessDat.protocol, 'irVisualIntegrationTask')
    %sessDat.data.stimuli_visir is correct
    corrIncorrVals = sessDat.data.stimuli_visir(1:numTrials) == sessDat.data.behavior(1:numTrials); 
else
    irOnlyVals = sessDat.data.ir_only_list(1:numTrials);   %0: ir+vis, 1: ir only
    corrIncorrVals = sessDat.data.stimuli(1:numTrials) == sessDat.data.behavior(1:numTrials); 
end

corrTrialInds = find(corrIncorrVals);
incorrTrialInds = find(~corrIncorrVals);
numTrialsCorrect = length(corrTrialInds);
numTrialsIncorrect = length(incorrTrialInds);

%Calculate PC based on plexon events
pcOverall = mean(corrIncorrVals);  %check this against sessDat
corrOnsetTimes = trialOnsetTimes(corrIncorrVals);
incorrOnsetTimes = trialOnsetTimes(~corrIncorrVals);
corrOffsetTimes = trialLastClockTimes(corrIncorrVals);
incorrOffsetTimes = trialLastClockTimes(~corrIncorrVals);

%compare pc calculated from session info
%For integrated, only want visual+ir (not vis only, or ir only)
sessBeh = sessDat.data.behavior(1:numTrials);
if isequal(sessDat.protocol, 'irVisualIntegrationTask')
    sessStim = sessDat.data.stimuli_visir(1:numTrials);
else
    sessStim = sessDat.data.stimuli(1:numTrials);
end
sessCorrIncorr = sessBeh == sessStim;
sessPCCalc = mean(sessCorrIncorr);
corrIncorrCompare = sum(sessCorrIncorr - corrIncorrVals);
if corrIncorrCompare
    error('Different pc calculated from sessdat and plexon data')
else
    display(['Bullseye! PCs match up from sessdat and plexon: ' num2str(pcOverall)])
end

%% Decode binary for all channels, pulling stim frequency indices, values, and clock times
%    for trial with no stim, stimSequences{trial} is NaN, and stimTime{} is just irOn onset of trial.
%
%stimSequencesClock{stimChan}{trial} contains 1 x sequenceLength array of substim frequencies for each trial 
%substimOnsetsClock{stimChan}{trial} contains 1 x sequenceLength array of substim onset times,
%   as calculated using clock events, for each  trial

jitter=0.002; 
stimSequencesClock=cell(1,4);
substimOnsetsClock =  cell(1,numTrials);
for chanInd = 1 : 4
    display(['Pulling stim from channel ' num2str(chanInd)]);
    tmpChannelDat = stimEvents{chanInd};
    [stimSequencesClock{chanInd}, substimOnsetsClock{chanInd}] = pull_binary_pulses_onechan(tmpChannelDat, stimClock, ...
        trialOnsetTimes, trialLastClockTimes , jitter);
end
allSubstimOnsetsClock = cell2mat(substimOnsetsClock{1}');

%% Plot stimulus values superimposed over events to check
if plotFull
    figure;
    chanVal = 'A';  %
    timeWindow = [200 350];  %for weird stuff in jasmin: maxFreq*2 870
    chanInd = find(ismember(chanVals, chanVal));
    tmpChannelDat = stimEvents{chanInd};
    eventTimes = {stimClock, tmpChannelDat};
    eventNames = {'Clock', ['Channel ' chanVal]};
    colors = {'m','b'};
    frequencyViewer(stimSequencesClock{chanInd}, substimOnsetsClock{chanInd}, ...
                    eventTimes, eventNames, trialFirstClockTimes, timeWindow, colors)
end

%% Create stimulation vector/population vector sequences for all channels/trials
%   stimvecSequencesInd{trialNum} : numSubstim x 4 (stim vector for each substim)
%   stimvecSequencesFreq{trialNum}: similar, but with actual frequencies delievered
%   popvecSequences{trialNum}: similar but contains numSubstim x 2 population vectors
%       note: popvec = sum(magnitude_i*direction_i)
%   allStimvecsInd: numStimOverall x 4 (stim vector combined over all trials into one matrix)
%   allStimvecsFreq: similar, but actual stimulus frequencies
%   allPopvecs : similar, but all population vectors (sum(magnitude_i*direction_i)
%   numSubstimAll: 1xnumTrials number of substimuli in each stim sequence
%   numSubstimTotal: scalar of all
%stimSpaceExplore.m explores the set of possible stimuli


%How universal is the following? need to improve
if isequal(sessDat.protocol, 'irVisualIntegrationTask')
    VOLT_DIVISIONS = [eps   1.11  6.46     8.86     9.84    10.35   10.76];
    FREQ_DIVISIONS = round(linspace(maxFreq,0,length(VOLT_DIVISIONS)));
    stimFreqs = sort(FREQ_DIVISIONS);
else
    stimFreqs = [0 10 100 150 250 350 425]; %cold: 0/1 warm: 100-250/ hot: 350+
end
    
stimInds = [0 1 2 3 4 5 6];
numFreqs = length(stimFreqs);
detectorCoords = [1 -1 -1 1;1 1 -1 -1];  %2x4, channel i has column i (A B C D)
stimvecSequencesInd = cell(1,numTrials);
stimvecSequencesFreq = cell(1,numTrials);
popvecSequences = cell(1, numTrials);
allStimvecsInd =[]; 
allStimvecsFreq = [];
allPopvecs = [];
numSubstimAll = zeros(1,numTrials); 
for trialNum = 1 : numTrials
    if isnan(stimSequencesClock{1}{trialNum})
        numSubstim = 1;
        stimvecInd = [0 0 0 0];
        stimvecFreq = [0 0 0 0];  
        popvecTmp = [0 0];
    else
        numSubstim = length(stimSequencesClock{1}{trialNum}); %pick just one instance to get length
        stimvecInd = zeros(numSubstim, 4);
        stimvecFreq = zeros(numSubstim,4);
        popvecTmp = zeros(numSubstim, 2);
        for subStim = 1 : numSubstim
            stimvecInd(subStim,:)= [stimSequencesClock{1}{trialNum}(subStim) stimSequencesClock{2}{trialNum}(subStim) ...
                        stimSequencesClock{3}{trialNum}(subStim) stimSequencesClock{4}{trialNum}(subStim)];
            for stimInd = 1:numFreqs
                indsVal = find(stimvecInd(subStim,:) == stimInds(stimInd));
                stimvecFreq(subStim,indsVal) = stimFreqs(stimInd);
            end   
            popvecTmp(subStim,:) =  detectorCoords * stimvecFreq(subStim, :)'; 
        end
    end
    
    numSubstimAll(trialNum) = numSubstim;
    stimvecSequencesInd{trialNum} = stimvecInd;
    stimvecSequencesFreq{trialNum} = stimvecFreq;
    popvecSequences{trialNum} = popvecTmp;
    allStimvecsInd = [allStimvecsInd; stimvecInd];
    allStimvecsFreq = [allStimvecsFreq; stimvecFreq];
    allPopvecs = [allPopvecs; popvecTmp];

        
        
end

numSubstimTotal = sum(numSubstimAll);
totalNumSubstim = sum(numSubstimAll);

%% Calculate entropy of stimulus frequency vectors
uniqueStimvecs = unique(allStimvecsFreq , 'rows');
numUniqueStimvecs = size(uniqueStimvecs,1);
stimvecProbabilities = zeros(1,numUniqueStimvecs);
stimvecCount = zeros(1,numUniqueStimvecs);
for stimvecNum = 1: numUniqueStimvecs
   stimvecCount(stimvecNum) = sum(ismember(allStimvecsFreq , uniqueStimvecs(stimvecNum,:), 'rows'));
   stimvecProbabilities(stimvecNum) = stimvecCount(stimvecNum)/totalNumSubstim;
end
positiveProbInds = find(stimvecProbabilities);
stimvecEntropy = sum(-stimvecProbabilities(positiveProbInds) .* log2(stimvecProbabilities(positiveProbInds))); %6.5


%% Calculate entropy of pop vecs
uniquePopvecs = unique(allPopvecs, 'rows');
numUniquePopvecs = size(uniquePopvecs,1);
popvecProbabilities = zeros(1,numUniquePopvecs);
popvecCount = zeros(1,numUniquePopvecs);
totalNumSubstim = sum(numSubstimAll);
for popvecNum = 1: numUniquePopvecs;
   popvecCount(popvecNum) = sum(ismember(allPopvecs, uniquePopvecs(popvecNum,:), 'rows'));
   popvecProbabilities(popvecNum) = popvecCount(popvecNum)/totalNumSubstim;
end
positiveProbInds = find(popvecProbabilities);
stimEntropy = sum(-popvecProbabilities(positiveProbInds) .* log2(popvecProbabilities(positiveProbInds)));


%% Pull frequent (>=30) substim for classification in neuronal analysis 
minCount = 30;
if plotFull
    stem(sort(popvecCount,'descend'))% recall this is unique popvecs only
    set(gca,'YTick', minCount);grid on
end
highProbInds = find(popvecCount >= minCount);
numHigh = length(highProbInds);
highProbPopvecs = uniquePopvecs(highProbInds,:);
%NO following is wrong highProbStimFreqs = allStimvecsFreq(highProbInds,:);

%Just do zero versus all the others, get a nice map of discrimination.
classifyZeroInds =  find(ismember(allPopvecs, [0 0],'rows'));
numZero = length(classifyZeroInds);
probZero = numZero/totalNumSubstim;

highZeroInd = find(ismember(highProbPopvecs, [0 0], 'rows'));
classifyPopvecs = highProbPopvecs;
classifyPopvecs(highZeroInd,:) =[];
numNonzeroHigh = numHigh - 1;
classifyPopvecInds = cell(1, numNonzeroHigh); %
for vecNum = 1: numNonzeroHigh
    popvecTmp = classifyPopvecs(vecNum,:);
    classifyPopvecInds{vecNum} = find(ismember(allPopvecs, popvecTmp, 'rows'));
end
allPopvecs(classifyPopvecInds{1},:);
popvecProbabilities(highProbInds);

%To save for classification with neuronal data:
%   classifyZeroInds, classifyPopvecInds, classifyPopvecs

%% Plot high probability ones
%Don't use anymore as rarely are they this high.
% if plotFull
%     figure;
%     scatter(highProbPopvecs(:,1), highProbPopvecs(:,2), 'k', 'filled');hold on
%     plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0],'k','LineWidth', 1); 
%     axis equal;title('high probability population vectors');grid on;gridfix([0.7 0.7 0.7])
% end

%% Pull and plot top X
[sortedPopvecCount, sortedPopvecInds] = sort(popvecCount, 'descend');
numTopPopvecs = 8;
topPopvecCounts = sortedPopvecCount(1:numTopPopvecs);
topPopvecs = uniquePopvecs(sortedPopvecInds(1:numTopPopvecs),:);
topPopvecInds = sortedPopvecInds(1:numTopPopvecs);
if plotFull
    figure
    scatter(topPopvecs(:,1), topPopvecs(:,2), 15, 'k', 'filled');hold on
    plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0], 'k','LineWidth', 1); 
    grid on; axis equal; 
    axis([-maxFreq*2 maxFreq*2 -maxFreq*2 maxFreq*2])
    set(gca,'XTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'XTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    set(gca,'YTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'YTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    gridfix([0.7 0.7 0.7])
    title('highest frequency popvecs')

    for topVecNum = 1:numTopPopvecs
        text(topPopvecs(topVecNum,1)+25, topPopvecs(topVecNum,2), num2str(topPopvecCounts(topVecNum)));
    end
end

%get four clusters, one including zero, and use that to find the furthest points?

%% Plot all possible popvecs, and observed population vectors
%All possible stimulus combinations/population vectors
allStimCombos=combvec(stimFreqs, stimFreqs, stimFreqs, stimFreqs)';
numCombos = size(allStimCombos,1);
%Population vector P:
%P = sum(magnitude_i*direction_i)
popvecSpace = zeros(numCombos, 2);
for comboNum = 1: numCombos
    tmpCombo = allStimCombos(comboNum,:);
    popvecSpace(comboNum,:) = detectorCoords*tmpCombo';   
end

if plotFull
    figure;
    plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0],'k','LineWidth', 1); hold on
    scatter(popvecSpace(:,1), popvecSpace(:,2), 25, 'filled', 'CData', [.6 .6 .6])
    scatter(allPopvecs(:,1), allPopvecs(:,2),25,'r', 'filled');
    grid on;
    title([animalName ' ' strrep(dateStr, '_', '/') ' . Entropy: ' num2str(stimEntropy) ' .N = ' num2str(numTrials)])
    xlabel('Mediolateral')
    ylabel('Anteroposterior');
    axis tight; axis equal
    axis([-maxFreq*2 maxFreq*2 -maxFreq*2 maxFreq*2])
    set(gca,'XTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'XTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    set(gca,'YTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'YTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    gridfix([.7 .7 .7])
    shg
    if printOn
       print
    end
end

%% Plot 3d frequency histogram of population vectors
%hist3 takes m x 2 matrix
if plotFull
    figure;
    frequencyBinCenters = [-maxFreq*2: maxFreq/4: maxFreq*2];
    zeroBinInd = find(frequencyBinCenters == 0);
%     binNumZero = find(frequencyBinCenters
    N =  hist3(allPopvecs, {frequencyBinCenters, frequencyBinCenters} );
    histVals = sort(N(:), 'descend');
    secondLargestHist = histVals(2);
    hist3(allPopvecs, {frequencyBinCenters, frequencyBinCenters} ); hold on
    xlabel('Mediolateral')
    ylabel('Anteroposterior');
    plot3([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0],10*[1 1 1 1 1],'k','LineWidth', 1); hold on
    plot3([-maxFreq*2 maxFreq*2],[0 0], [10 10],'r')
    plot3([0 0],[-maxFreq*2 maxFreq*2], [10 10],'r')
    set(gca,'ZLim',[0 secondLargestHist])
    set(gca,'XTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'XTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    set(gca,'YTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'YTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    title('Population vector frequency histogram')
    set(gcf,'Renderer','OpenGL')   
    %if printOn
    %    print
    %end
end


%% Plot distribution of stimsequence length per trial(basically this is reaction time)
%break up by corr/incorr
if plotFull
    figure;
    maxSubstim = max(numSubstimAll);
    hnum = hist(numSubstimAll, [0:maxSubstim]);
    bar([0:maxSubstim], hnum,1); axis tight;shg
    [mnNumSubstim, semNumSubstim] = mean_sem(numSubstimAll'); %
    xlabel('Num substim');ylabel('Count');
    title(['Num substimuli per trial. Mean ' num2str(mnNumSubstim) ', SEM: ' num2str(semNumSubstim)])
end


%% Calculate histograms of stimulus frequency versus chan num
allstimvec=allStimvecsInd(:);
yCents = [0:6];
ha=hist(allstimvec, [0:6]);
h1=hist(allStimvecsInd(:,1),[0:6]); h1 = h1/sum(h1);
h2=hist(allStimvecsInd(:,2),[0:6]); h2 = h2/sum(h2);
h3=hist(allStimvecsInd(:,3),[0:6]); h3 = h3/sum(h3);
h4=hist(allStimvecsInd(:,4),[0:6]); h4 = h4/sum(h4);

hmax = max([h1 h2 h3 h4]);

if plotFull
    %Get said distribution for four different channels
    figure
    subplot(4,1,1)
    h=bar(yCents, h1, 1); hold on
    set(h,'FaceColor','k', 'EdgeColor','k')     
    axis([-.5 6.5 0 hmax+hmax*.1])
    set(gca,'XTick', [0:6], 'YTick',[0:.1: ceil(hmax*10)/10] );
    title(chanVals{1})

    subplot(4,1,2)
    h=bar(yCents, h2, 1); hold on
    set(h,'FaceColor','k', 'EdgeColor','k')     
    set(gca,'XTick', [0:6], 'YTick',[0:.1: ceil(hmax*10)/10] );
    axis([-.5 6.5 0 hmax+hmax*.1])
    title(chanVals{2})

    subplot(4,1,3)
    h=bar(yCents, h3, 1); hold on
    set(h,'FaceColor','k', 'EdgeColor','k')     
    set(gca,'XTick', [0:6], 'YTick',[0:.1: ceil(hmax*10)/10] );
    axis([-.5 6.5 0 hmax+hmax*.1])
    title(chanVals{3})

    subplot(4,1,4)
    h=bar(yCents, h4, 1); hold on
    set(h,'FaceColor','k', 'EdgeColor','k')     
    set(gca,'XTick', [0:6], 'YTick',[0:.1: ceil(hmax*10)/10] );
    axis([-.5 6.5 0 hmax+hmax*.1])
    title(chanVals{4})

    suplabel('Histogram of ICMS frequencies', 'x');
end

%%
save workingAnalysis1
%load workingAnalysis1




%% Get all rows where first and second are 350 or so
%Haven't done this in a while
% Plot different stimulus types (recall stimFreqs = [0 10 100 150 250 350 425])
%#? not sure what the following is
%Plot frequency histogram of all possibly stimuli (6^4 = 1296 possible stimuli = 10.3 bits if equiprobable)
%Plot frequncy of top fo
%Number of guys with [N N 0 0]
% freqVal1 = 350;
% freqVal2 = 350;
% freqInds1 = find(ismember(allStimvecsFreq ,[freqVal1 freqVal1 0 0],'rows'));  %pinky 11/14: 21 (66 for 425)
% freqInds2 = find(ismember(allStimvecsFreq ,[freqVal1 freqVal1 10 0],'rows'));  
% freqInds3 = find(ismember(allStimvecsFreq ,[freqVal1 freqVal1 0 10],'rows'));  
% freqInds4 = find(ismember(allStimvecsFreq ,[freqVal1 freqVal1 10 10],'rows'));  
% % freqInds5 = find(ismember(allStimvecsFreq ,[freqVal2 freqVal2 0 0],'rows'));  %pinky 11/14: 21 (66 for 425)
% % freqInds6 = find(ismember(allStimvecsFreq ,[freqVal2 freqVal2 10 0],'rows'));  
% % freqInds7 = find(ismember(allStimvecsFreq ,[freqVal2 freqVal2 0 10],'rows'));  
% % freqInds8 = find(ismember(allStimvecsFreq ,[freqVal2 freqVal2 10 10],'rows'));  
% 
% allNum = [length(freqInds1), length(freqInds2), length(freqInds3), length(freqInds4)];, ...
%     %length(freqInds5), length(freqInds6), length(freqInds7), length(freqInds8)];
% sum(allNum);


%%%%%***Analysis of individual population vector trajectories****%%%%%%%%
trialNum = 1;  %need to have this separate so it doesn't reset every time
%% Plot selection of population vector trajectories
%Note doesn't plot final substim
if plotPopvecCycles
    figure;
    keepGoing = 1;
    while keepGoing   
        popVec =  popvecSequences{trialNum};
        numSubstim = size(popVec,1);
        colorVals = linspace(.95,0,numSubstim);
        for subStimNum = 1: numSubstim
            colSub=colorVals(subStimNum)*ones(1,3);
            scatter(popVec(subStimNum, 1), popVec(subStimNum, 2), 'filled', 'CData', colSub);hold on
            if subStimNum < numSubstim
                plot([popVec(subStimNum, 1), popVec(subStimNum+1, 1)], ...
                    [popVec(subStimNum, 2), popVec(subStimNum+1, 2)],'Color',colSub);
            end
            %text(popVec(1, subStimNum)+.025, popVec(2, subStimNum), num2str(subStimNum));
        end
        %set(gca,'XTick',0, 'YTick',0)
        plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0],'k','LineWidth', 1); %the square
        grid on
        axis equal;shg

        title(['Trial ' num2str(trialNum) '/' num2str(numTrials) ': corr/incorr ' num2str(corrIncorrVals(trialNum))] )
        trialNum = trialNum + 1;
        
        gridfix([.9 .9 .9])
        %button = questdlg('Continue?','checking plots','yes','no', 'yes');
        button = MFquestdlg([.4 .35], 'Continue?', 'Continue?', 'yes', 'no', 'yes');
        if isequal(button, 'no')
            keepGoing = 0;
        else
            close;
        end
    end %while
end

%% Plot a few popvec trajectories: this is ugly: for figures will want to show individuals
% %   individual trials
% %Not working: seems to be overwriting...
% plotTrials=[18 28 29 33];
% numToPlot = length(plotTrials);
% for trialNum = 1: numToPlot
%     trialNum 
%     trialToPlot = plotTrials(trialNum);
%     popVec=popvecSequences{trialToPlot};
%     numSubstim = size(popVec,1);
%     colorVals = linspace(.95,0,numSubstim);
%     for subStimNum = 1: numSubstim
%         colSub=colorVals(subStimNum)*ones(1,3);
%         scatter(popVec(subStimNum, 1), popVec(subStimNum, 2), 'filled', 'CData', colSub);hold on
%         if subStimNum < numSubstim  %do this to get color different for each line
%             plot([popVec(subStimNum, 1), popVec(subStimNum+1, 1)], ...
%                 [popVec(subStimNum, 2), popVec(subStimNum+1, 2)],'Color',colSub);
%         end
%         %text(popVec(1, subStimNum)+.025, popVec(2, subStimNum), num2str(subStimNum));
%     end
% end     
% 
%     plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0],'k','LineWidth', 1); %the square
%     grid on
%     axis equal;shg;gridfix([.9 .9 .9])
% shg    
% %print -depsc2 cmSamples


%% Plot all population vector trajectories
%Note if it won't print:
%   set(gcf,'Renderer','OpenGL')
%All: green correct/red incorrect
if plotPopvecs | plotFull
    % plot all
    plot_stimulus_population(popvecSequences, animalName, maxFreq, dateStr); 
    
%     % plot correct
%     popvecSequencesSubset ={};
%     for trialInd = 1:numTrialsCorrect
%         popvecSequencesSubset{end+1}=popvecSequences{corrTrialInds(trialInd)};
%     end
%     plot_stimulus_population(popvecSequencesCorrect,animalName, maxFreq, dateStr);
%     
%     % plot incorrect
%     popvecSequencesSubset ={};
%     for trialInd = 1:numTrialsIncorrect
%         popvecSequencesSubset{end+1}=popvecSequences{incorrTrialInds(trialInd)};
%     end
%     plot_stimulus_population(popvecSequencesSubset,animalName, maxFreq, dateStr);
%     
    if printOn
        print
    end
end

%% Following is if you are printing for paper
%{ 
    cd C:\Users\Eric\Dropbox\IR6\S1V1\s1_v1_paper\Figures\s1v1Compare
    print -depsc2 popVecV1
    cd(datFolder) 
%}
    
%set(gcf,'Renderer','OpenGL')    %if printer won't work


%% save all data, and subset of data for subsequent neuronal analysis
saveString = ['save stimResults_' animalName '_' dateStr];  %stimResults_animalName_date (e.g., stimResults_pinky_11_11_14)
eval(saveString);


%For neuronal analysis
% corrTrialInds             1 x numCorr array with indices of correct trials
% incorrTrialInds           1 x numIncorrect array with indices of incorrect trials
% corrIncorrVals            1 x numTrials with 0 for incorrect, 1 for correct       
% trialOnsetTimes           1 x numTrials time at which IR light turned on   
% allSubstimOnsetsClock     totalNumSubstim x 1 (clock-based onset time of each substim)
% substimOnsetsClock        1 x numTrials cell array, each contains 1 x numSubstim array of onset times
% trialLastClockTimes       1 x numTrials final clock tic for each trial
% trialFirstClockTimes      1 x numtrials first clock tic for each trial
% stimvecSequencesInd       1 x numTrials cell array, each contains numSubstim x 4 vector of stim index
% stimvecSequencesFreq      1 x numTrials cell array, each contains numSubstim x 4 vector of stim frequency
% allStimvecsInd            totalNumSubstim x 4: vector of stim indices
% allStimvecsFreq           totalNumSubstim x 4: vector of stim frequencies
% popvecSequences           1 x numTrials cell array: each is numSubstim x 2 population vector
% allPopvecs                totalNumSubstim x 2: all population vectors
% numSubstimAll             1 x numTrials array of number of substim on each trial
% stimFreqs                 1 x 7 stimulus frequencies used for microstimulation
% Start, Stop
cd(datFolder)    
saveNeuroString = ['save stimAnalysisNeuro_' animalName '_' dateStr ...
            ' corrTrialInds incorrTrialInds corrIncorrVals allSubstimOnsetsClock ' ...
            ' substimOnsetsClock trialOnsetTimes trialLastClockTimes trialFirstClockTimes ' ...
            ' stimvecSequencesInd stimvecSequencesFreq allStimvecsInd allStimvecsFreq ' ...
            ' popvecSequences allPopvecs numSubstimAll Stop Start stimEntropy' ...
            ' classifyZeroInds classifyPopvecInds classifyPopvecs totalNumSubstim ' ...
            ' topPopvecCounts numTopPopvecs topPopvecs stimFreqs'];
eval(saveNeuroString);

%display('Done evaluating stimulus!')



