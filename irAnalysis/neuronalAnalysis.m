%Initial analysis of neuronal data: psth after stim, mean response versus
%   stimulus frequency, and other basic statistics, plot contour/surf plots
%   get "receptive field" center.
%
%With plotFull on, plots everything: [fill this in]
%
%With plotMinimal on, plots:
%   psth
%   resp versus number of stim sites
%   contour plot
%   discrimination results
%
%With printOn, will print:
%
%
%Pull relevant data from stimulus analysis:
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
%
% To do or think about
%     Get offset of artifacts at end of each trial, because they are usually cut short (i.e., the 
%       stimulation doesn't last the full 75 ms, or whatever).
%REsponseInds depends on stim duration and period! 
% responseInds = [83:88; 96:101; 109:114; 122:127; 135:140; 148:153];

%% Load results of stimulus analysis, and load/plot artifact events
close all
cd(datFolder)
stimAnalysisFilename = ['stimAnalysisNeuro_' animalName '_' dateStr];  
load(stimAnalysisFilename)

%Pull artifact event data from artifactChannels
cd(datFolderNeuro)
load neuronalEvents
numArtifactChannels = length(artifactChannels);
artifactEventNames = cell(1, numArtifactChannels);
artifactEvents = cell(1, numArtifactChannels);
for artifactChannelInd = 1 : numArtifactChannels
    artifactChannel = artifactChannels(artifactChannelInd);
    if artifactChannel < 10
        chanName = ['sig00' num2str(artifactChannel) 'a'];
    elseif artifactChannel < 100
        chanName = ['sig0' num2str(artifactChannel) 'a'];
    end
    artifactEventNames{artifactChannelInd} = chanName;  
    artifactEvents{artifactChannelInd} = eval(chanName)';
end

%%%
%Plot artifact events, clock events, and stim on/off times
artifactsStimOn = artifactEvents;
artifactsStimOn{numArtifactChannels+1} = allSubstimOnsetsClock;
artifactsStimOn{numArtifactChannels + 2} = trialOnsetTimes;
artifactsStimOn{numArtifactChannels + 3} = trialLastClockTimes;
eventColorsArtifact ={'k','c', 'm', 'b','g','r', 'c', 'g', 'b', 'r', 'k', 'c', 'm', 'k', 'b'};
if ~(length(eventColorsArtifact) >= numArtifactChannels + 3)
   error('Need to match up color array to number of artifact channels');
end
eventNamesArtifact = artifactEventNames;
eventNamesArtifact{numArtifactChannels+1} = 'substim Onset Clock';
eventNamesArtifact{numArtifactChannels+2} = 'Trial Onset Times';
eventNamesArtifact{numArtifactChannels+3} = 'Trial Last Clock Times';

timeWindowFull = [Start Stop];
timeSubwindow = [1000 1050]; %[173 176]; %[1028.55 1028.9];
if plotFull
    figure;
    eventViewer(artifactsStimOn, eventNamesArtifact, timeSubwindow, eventColorsArtifact); hold on
end

%% Pull super-fine-grained artifact psth for each substim 
%(basically makes a very fine-grained psth--this is not standard psth)
%
% allArtifactPsths{trialNum}  is a 1 x numSubStim cell array, each substim elt contains 
%   a numArtifactChannels x numBins psth.
% sumArtifactPsths{trialNum} is a 1 x numSubstim cell array, eahc contains 1 x numBins psth
%   the sum of the numArtifactChannel psths in allArtifactPsths
% allArtifactBinEdges{trialNum} is 1xnumSubStim cell array that gives bin edges for the psths 
%   contained in allArtifactPsths
%
%Note: they can have different bin numbers because it makes psth from (putative) onset of artifact until
%  onset of next substim (or 0.2 seconds, if it is the last substim in a trial), and this is a variable # b/c matlab

allArtifactTimes = cell2mat(artifactEvents);
%Pull stimulus artifact onset time and offset time after each pulse
allArtifactPsths = cell(1,numTrials); %each trial is numSubstim x 2 x numArtifactChannels
allArtifactBinEdges = cell(1, numTrials);
sumArtifactPsths = cell(1,numTrials);
artifactBinWidth = .0015;   %to count as onset, majority must have artifact within jitter ms of each other
artifactPsthDuration = 0.2; %looks in this time window if last substim in trial

for trialNum = 1: numTrials
    if rem(trialNum, 50) == 0
        disp(['Making artifact psth trial ' num2str(trialNum) ' out of ' num2str(numTrials)]);
    end
    trialOnsetTime = trialOnsetTimes(trialNum);
    substimOnsetTimesTmp = substimOnsetsClock{1}{trialNum};
    numSubstim = length(substimOnsetTimesTmp);
    artifactSubstimPsths = cell(1, numSubstim);   
    artifactSubstimPsthSums = cell(1,numSubstim);
    allSubstimBinEdges = cell(1, numSubstim);
    
        for substimNum = 1: numSubstim
            tmpSubstimFreqs = stimvecSequencesFreq{trialNum}(substimNum,:);
            %Pull ir light onset time, and time of next ir onset if there is one
            substimOnTimeTmp = substimOnsetTimesTmp(substimNum);
            
            if substimNum < numSubstim %if it isn't the last substimulus in the sequence
                nextSubstimOnsetTime = substimOnsetTimesTmp(substimNum + 1); %except when numSubstim
            else  %for the last substimulus in the sequence
                nextSubstimOnsetTime =trialLastClockTimes(trialNum) + artifactPsthDuration;
            end
            
            artifactsAbove = allArtifactTimes(find(allArtifactTimes >= substimOnTimeTmp));
            firstArtifactTime = min(artifactsAbove);
            substimBinEdges = [firstArtifactTime - .0001: artifactBinWidth: nextSubstimOnsetTime];
            allSubstimBinEdges{substimNum} = substimBinEdges;
            numBins = length(substimBinEdges);
            artifactSubstimPsth = zeros(numArtifactChannels, numBins);
            
            for artifactChannel = 1: numArtifactChannels
                tmpChanArtifactInds = ...
                    find(artifactEvents{artifactChannel} > substimOnTimeTmp & ...
                        artifactEvents{artifactChannel} < nextSubstimOnsetTime);
                if ~isempty(tmpChanArtifactInds)
                    tmpChanArtifactTimes = artifactEvents{artifactChannel}(tmpChanArtifactInds);
                    tmpHist = histc(tmpChanArtifactTimes, substimBinEdges);
                    artifactSubstimPsth(artifactChannel,:) = tmpHist >= 1;
                end
            end %for artifact channel
            
            artifactSubstimPsths{substimNum} = artifactSubstimPsth;
            artifactSubstimPsthSums{substimNum} = sum(artifactSubstimPsth);
        end %for substimNum
        
        allArtifactPsths{trialNum} = artifactSubstimPsths;
        allArtifactBinEdges{trialNum} = allSubstimBinEdges;
        sumArtifactPsths{trialNum} = artifactSubstimPsthSums;
end %for trialNum
disp('Done pulling artifact psth');

%%Following is for looking at individual cases you suspect are messed up
%Pretty much is for debugging purposes
% % %% Plot trial artifact psth (60/7: [1028.55 1028.9]), (8/3: [173 176])
% if plotFull
%     trial =  50; 
%     substim = 4;
%     dat = sumArtifactPsths{trial}{substim};
%     filtVariance = 1;
%     numels = 3;
%     gauss_smoothed = filtfilt(gausswin(3, 1),sum(gausswin(numels, filtVariance)),dat);
%     plot(gauss_smoothed,'r');hold on; 
%     bar(dat);axis tight;
%     title(['Artifact psth, trial ' num2str(trial) ', substim ' num2str(substim)])
% end

%% Pull artfifact onset/offset times for each substim (use Use allArtifactPsths and artifactPsthBinEdges)
%Recall what we just made
% sumArtifactPsths{trialNum} is a 1 x numSubstim cell array, eahc contains 1 x numBins psth
% allArtifactBinEdges{trialNum} is 1xnumSubStim cell array that gives bin edges for the psths 
%
%What we are now making:
% artifactWindows{trialNum} = numSubstim x 2 window of artifact times
% artifactOnTimes(substimNum)  1 x totalNumSubstim includes all onset times
% artifactOffTimes(substimNum) 1 x totalNumSubstim includes all offset times
% artifactDurations(substimNum) 1 x totalNumSubstim duration of artifact (though this is sort of moot)
%Trial 60, substim 7 is a bitch   (pinky 11/5/14)
%   substimOnsetsClock{60}(7)
artifactWindows = cell(1,numTrials);
artifactMajority = (1/2)*numArtifactChannels; %majority must have artifact for it to count as onset.
artifactOnTimes = [];
artifactOffTimes =[];
artifactDurations =[];

for trialNum = 1 : numTrials
    numSubstim = numSubstimAll(trialNum);
    substimOnsetTimesTmp = substimOnsetsClock{1}{trialNum};
    trialWindows = [];
    if numSubstim == 0
        artifactOnTime = substimOnsetTimesTmp +.005; 
        artifactOffTime = artifactOnTime + artifactDuration;  
        artifactOnTimes = [artifactOnTimes; artifactOnTime];
        artifactOffTimes = [artifactOffTimes; artifactOffTime];
        artifactDurations = [artifactDurations; artifactOffTime - artifactOnTime];
        trialWindows = [trialWindows; [artifactOnTime artifactOffTime]];
    else
        for substimNum = 1: numSubstim
            tmpSubstimFreqs = stimvecSequencesFreq{trialNum}(substimNum,:);
            %First check if there was no stimulation, in which case set manually
            %Otherwise 
            if isequal(tmpSubstimFreqs, [0 0 0 0])
                artifactOnTime = substimOnsetTimesTmp(substimNum)+.005; 
            else
                tmpHist = sumArtifactPsths{trialNum}{substimNum};
                tmpHistBinEdges = allArtifactBinEdges{trialNum}{substimNum};
                histThreshCrossBin = find(tmpHist >= artifactMajority,1,'first');   %what bin, in the hist, does it cross?
                if isempty(histThreshCrossBin)
                    %disp('No majority')
                    artifactOnTime = substimOnsetTimesTmp(substimNum)+.005; 
                else   
                    artifactOnTime = tmpHistBinEdges(histThreshCrossBin);
                end           
            end %check substim onset artifact time
            artifactOffTime = artifactOnTime + artifactDuration;  
            artifactOnTimes = [artifactOnTimes; artifactOnTime];
            artifactOffTimes = [artifactOffTimes; artifactOffTime];
            artifactDurations = [artifactDurations; artifactOffTime - artifactOnTime];
            trialWindows = [trialWindows; [artifactOnTime artifactOffTime]];
        end %substim
    end
        
    artifactWindows{trialNum} = trialWindows;  
end


%% Look at basic properties of artifact substim windows (max duration, time between artifacts).

%Pull max artifact duration
maxArtifactDuration = max(artifactDurations);
%Calculate and plot time between artifacts just to get a sense
%overall (this is not very meaningful)
betweenArtifactAll = artifactOnTimes(2:end) - artifactOffTimes(1:end-1) ;
[mnDurAll, semDurAll] = mean_sem(betweenArtifactAll);
% plot(betweenArtifactAll);axis tight;shg;hold on
% set(gca,'YLim', [-1 max(betweenArtifactAll)])
% badSubstimDurs = find(betweenArtifactAll < .06) %which means between offset 755 and onset 756
% title('overall inter-artifact times')

%within each stimsequence
betweenArtifactTimes =[];
for trialNum = 1:numTrials
    artifactWindowsTmp = artifactWindows{trialNum};
    if ~isempty(artifactWindowsTmp)
        betweenTimes =  artifactWindowsTmp(2:end,1) - artifactWindowsTmp(1:end-1,2);
    else
        trialNum
        betweenTimes = NaN
    end
    betweenArtifactTimes = [betweenArtifactTimes ; betweenTimes];

end
[mnDur, semDur] = mean_sem(betweenArtifactTimes); %.14/ 9e-4
disp(['Max artifact duration is ' num2str(maxArtifactDuration)])
minimumInterArtifactTime = min(betweenArtifactTimes);
disp(['Min inter-artifact time in stim sequences: ' num2str(minimumInterArtifactTime)])
disp(['Mean/sem of inter-artifact times: ' num2str(mnDur) ' / ' num2str(semDur)]) 

if plotFull
    figure;
    subplot(2,1,1)
    plot(betweenArtifactTimes); hold on
    title('Inter-arfifact times (within sequence)')
    axis tight;
    subplot(2,1,2)
    hist(betweenArtifactTimes);hold on %
    title('Inter-artifact times (within sequence)')
end


     
%% Eyeball check: Plot artifact windows on top of data: general reasonableness check
if plotFull
    figure
    timeSubwindow = [100 150]; %[1027 1030];  %pinky 11/5 1028.56 it's an ugly one...
    disp(['Plotting artifact in window ' num2str(timeSubwindow) '. This can take a few seconds.'])
    numEvents = length(eventNamesArtifact);
    onArtifacts = artifactOnTimes(find(artifactOnTimes > timeSubwindow(1) & ...
                                       artifactOnTimes < timeSubwindow(2)));
    offArtifacts = artifactOffTimes(find(artifactOffTimes > timeSubwindow(1) & ...
                                       artifactOffTimes < timeSubwindow(2)));
    eventViewer(artifactsStimOn, eventNamesArtifact, timeSubwindow, eventColorsArtifact); hold on
    line([onArtifacts'; onArtifacts'], [zeros(1,length(onArtifacts)); ...
        numEvents*ones(1,length(onArtifacts))], 'Color', 'g', 'LineWidth', 0.5, 'LineStyle', ':'); 
    line([offArtifacts'; offArtifacts'], [zeros(1,length(offArtifacts)); ...
        numEvents*ones(1,length(offArtifacts))], 'Color', 'm', 'LineWidth', 0.5, 'LineStyle', ':'); 
    title('Artifact bounds: on (green), off (magenta)')
end

cd(datFolderNeuro)
save(['neurAnalysisWorking1'])
%load(['neurAnalysisWorking1']);


%% Pull channel names
workspace_variables=who;
neur_suffix={'sig'};
[allChanNames, foo]=get_names(neur_suffix, workspace_variables);  
chanNamesMultiunit = {};
for channelNumber = 1: length(allChanNames)
    if allChanNames{channelNumber}(end) == 'i'
        chanNamesMultiunit{end+1} = allChanNames{channelNumber};
    end
end
numChannels = length(chanNamesMultiunit);
channelNumbers = pullChanNumbersMultiunit(chanNamesMultiunit);


%% Calculate full-trial psths and calculate mean/std from baseline period
%This is where channelPsthPlotScript is called, which is itself a nontrivial script
%It does the psth for a single channel

%indices for the six response epochs (note will need to change if you change numSubstimMin)
%#? Is this locked into bin width too? Shouldn't this be automated?
%This is going to be fubar for ford, where we have different duration stim!

% How did I even get these response Inds? 
%   Go into channelPsthPlotScript for nitty-gritty of psth bits until now you just have 
%   things in units of windows and times, not bin indices. You shoudl probably be calculating
%   "respondInds" inside the psth guts
responseInds = [83:88; 96:101; 109:114; 122:127; 135:140; 148:153];  %these should be calculated in channelPsthPlotScript not set!
%Calculate based on countDuration, artifact offsets, maybe other bits? Well, you represent stimuli and artifact as identical, which
% is why you did this, no?

baselineMeans = zeros(1, numChannels);
allBaselineCounts = zeros(numTrials, numChannels);
allSubstimCounts = zeros(totalNumSubstim, numChannels);
mnPsthZAll = zeros(numChannels, numSubstimMin);
plotIndividualChannels = 0;  %this is part of channelPsthPlotScript don't want to plot every psth detail
for channelNumber = 1 : numChannels
    channelName = chanNamesMultiunit{channelNumber};
    disp(['Calculating psth ' channelName])
    
    channelPsthPlotScript;  %script gets fullMeanPsth, fullMeanPsthZ
    for responseNum = 1: numSubstimMin
        mnPsthZAll(channelNumber, responseNum) = mean(fullMeanPsthZ(responseInds(responseNum,:)));
    end
    baselineMeans(channelNumber) = mean(baselineCounts);
    allBaselineCounts(:, channelNumber) = baselineCounts';
    allSubstimCounts(:,channelNumber) = substimCounts;
end

if plotMinimal | plotFull
    suplabel([animalName ' ' strrep(dateStr, '_', '/') ' psth']);shg
    if printOn
        print
    end
    disp('Done plotting psths');
end


%% for creating figures
%{
    cd C:\Users\Eric\Dropbox\IR6\S1V1\s1_v1_paper\Figures\s1v1Compare
    print -depsc2 peterPSTH
    cd(datFolder)
%}

%% Plot mean 3x3 population vector trajectory (this only needs doing once each session!)
if plotMinimal | plotFull
    figure;
    colorVals = linspace(.75,0,numSubstimMin);
    popVec = mnNPopvecs;
    for substimNum = 1: numSubstimMin
        %display(substimNum)
        %colSub=[colorVals(substimNum) colorVals(substimNum) colorVals(substimNum)]; 
        colSub=[1 colorVals(substimNum) colorVals(substimNum)]; %for paper to superimpose
        scatter3(popVec(substimNum,1), popVec(substimNum,2), 10, 45, 'filled', 'CData', colSub);hold on
        if substimNum < numSubstimMin
           if substimNum ~= numSubstimMin/2
               plot3([popVec(substimNum,1), popVec(substimNum+1,1)], ...
                   [popVec(substimNum,2), popVec(substimNum+1,2)], [10 10], 'Color',colSub,'LineWidth',2);
           else  %dashed line connected 3..4
               plot3([popVec(substimNum,1), popVec(substimNum+1,1)], ...
                  [popVec(substimNum,2), popVec(substimNum+1,2)], [10 10], ':','Color',colSub,'LineWidth',1);
           end
        end
    end
    
    plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0],'k','LineWidth', 1); %the square
    title(['Mean trajectory start/end ' animalName ' ' strrep(dateStr, '_', '/')])
    grid on; axis equal;
    axis([-maxFreq*2 maxFreq*2 -maxFreq*2 maxFreq*2])
    set(gca,'XTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'XTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    set(gca,'YTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
        'YTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
    gridfix([.85 .85 .85])  
    shg;
    if printOn
        set(gcf,'Renderer','OpenGL')  
        print
        pause(0.1)
    end
end


%{
%to save



%}

%% Calculate one-neuron receptive fields
%For old filtering properties (width=30, sigma = 10), see neuronalAnalysisFilterCheck.m
%Recall old stimFreqs before we messed with it: [0 10 100 150 250 350 425];

%Base RF filter/spacing properties on maximum frequencies used...
limVal = (maxFreq*2)+25;
rfSpacing = round((maxFreq*2)/170);  %Small spacing--5 for old
filterSigma = round((maxFreq*2)/85); %Not huge--10 for old;
filterWidth = round((maxFreq*2)/28); %Pretty big--30 for old 
imageFilter=fspecial('gaussian', filterWidth, filterSigma);
%figure;surf(imageFilter);

rfCenters = zeros(numChannels,2);
rfDiameters = zeros(numChannels,1);
peakZScores = zeros(numChannels, 1);
if plotMinimal
    figure;
end
for channelNum = 1:numChannels
    %channelNum = 11
    channelName = chanNamesMultiunit{channelNum};
    disp(['RF calculation ' channelName]);
    channelBaseMean = baselineMeans(channelNum);
    baselineCounts = allBaselineCounts(:, channelNum);
    substimCounts = allSubstimCounts(:, channelNum);

    plot_on = 0;
    channelData = [allPopvecs allSubstimCounts(:, channelNum)];
    baselineCounts = allBaselineCounts(:, channelNum);
    [mnBaselineCount, stdBaselineCount] = mean_std(baselineCounts);

    %Interpolative calculation of RF over support set
    [rfVals, rfFunction, mnResponsesUnique, xGrid, yGrid] = rfCalculate(channelData, rfSpacing, limVal, plot_on);
    %
    xVals = xGrid(1,:)';
    yVals = yGrid(:,1);
    
    %binWidthHists = (maxFreq*2)/4;
    %binCenters  = [-maxFreq*2: binWidthHists:maxFreq*2];  %85
    %histActualStim = hist3(totalPopvecs, {binCenters, binCenters} ); 
    %hist3(totalPopvecs, {binCenters, binCenters} ); 
   
    %{
    h = surf(xGrid, yGrid, rfVals);hold on
    set(h, 'edgecolor','none') %remove lines
    shg
    %}
    rfValsZ = (rfVals - mnBaselineCount)/stdBaselineCount;
    %%Filter it
    plot_on = 0; %plot_on = 1;
    rfValsFiltered = imFilterNan(rfValsZ, imageFilter, plot_on);
    %{
    figure; contourf(xVals, yVals, rfValsZ); hold on; xlabel('x'); ylabel('y'); 
    plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0], 'k','LineWidth', 1); axis tight; axis equal; 
    colorbar
    %}
    
    %{
    figure; contourf(xVals, yVals, rfValsFiltered); hold on;
    xlabel('x'); ylabel('y'); 
    plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0], 'k','LineWidth', 1); axis tight;axis equal;
    %}
    
     peakZScores(channelNum) = max(rfValsFiltered(:));
    
    % Surf plot of response magnitude versus location
    %if plotFull
    %    surf(xVals, yVals, rfValsFiltered); 
    %    hold on
    %    %For mean you'd use: h=surf(binCenters, binCenters, mnFiltered');hold on;
    %    shading interp;
    %    xlabel('Mediolateral')
    %    ylabel('Rostrocaudal'); 
    %    %plot3([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0], [0 0 0 0 0], 'w','LineWidth', 1); %the square
    %    title(['Median ' animalName ' ' strrep(dateStr, '_', '/') ' all. Neuron ' channelName]);shg
    %    colorbar
    %end

    % Get location of max, and value at max
    rfMaxVal = max(rfValsFiltered(:));
    [yMaxInd, xMaxInd] = find(rfValsFiltered == rfMaxVal);
    xMaxVal = xVals(xMaxInd);
    yMaxVal = yVals(yMaxInd);
    maxCoordinate = [xMaxVal yMaxVal];
    rfCenters(channelNum,:) = maxCoordinate;
    
    % Get surface area at 0.75 max.
    rfHeightContour = 0.75;
    roundDigits = 4;
    if rfMaxVal <= 0  %negative values throws things off
        rfValsFiltered = rfValsFiltered + max(abs(rfValsFiltered(:)));
        rfMaxVal = max(rfValsFiltered(:));
    end
    rfMinVal = min(rfValsFiltered(:));
    rfValsFiltered = rfValsFiltered - rfMinVal;
    
    contourHeights = round([0 .1 .25 .5 .75 .95]*rfMaxVal*10^roundDigits)/10^roundDigits;  
    C = contourc(xVals, yVals, rfValsFiltered, contourHeights);
    [ A ] = contourArea( C );
    Arounded = round(A(1,:)*10^roundDigits)/10^roundDigits;
    Aind = find(Arounded == round(10^roundDigits*(rfMaxVal*rfHeightContour))/10^roundDigits);
    rfHeightAreas = A(2,Aind);
    rfHeightArea = max(rfHeightAreas); %in case there are multiple peaks
    rfDiameter = 2*sqrt(rfHeightArea/pi);
    
    if plotMinimal | plotFull
        subplot(4,4,channelNum);
        contour(xVals, yVals,  rfValsFiltered, contourHeights);hold on
        plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0], 'k','LineWidth', 1); %the square
        %title([animalName ' ' strrep(dateStr, '_', '/') '. Neuron ' channelName '. Diam=' num2str(rfDiameter)])
        scatter(xMaxVal, yMaxVal, 10, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerFaceColor', 'r');
        set(gca,'XTick', [], 'YTick', []);
        axis equal;
        %colorbar
        %plot_circle([xMaxVal, yMaxVal], rfDiameter/2, 'r', 1) ;shg
        axis tight;
        title(channelName);shg
    end
    rfDiameters(channelNum) = rfDiameter;

end
if plotMinimal | plotFull
    suplabel(['Filter size/sigma: ' num2str(filterWidth) ' ' num2str(filterSigma) ' ' animalName ' ' strrep(dateStr, '_', '/') ' contour']);shg
    if printOn
       print
    end
    disp('Done plotting RFs');
end

%% Peak z scores
if plotFull
    [mnPeakZ, semPeakZ] = mean_sem(peakZScores);
    [medianPeakZ, seMedianPeakZ] = median_sem(peakZScores, 500);
    figure;
    plot([0 numChannels+1],[mnPeakZ mnPeakZ], 'r');hold on
    plot([0 numChannels+1],[medianPeakZ medianPeakZ], 'b')
    legend({'Mean', 'Median'}, 'Location', 'NorthWest')
    scatter([1: numChannels], sort(peakZScores,'ascend'),'k','filled');hold on
    axis([0 numChannels+1 -1 max(peakZScores)+0.1*max(peakZScores)])
    plot([0 numChannels+1],[0 0], 'k')
    title([animalName ' ' strrep(dateStr, '_', '/') ' peak z score (from filtered RF)'])
    if printOn
        print
    end
end


%% Calculate mean/sem spike counts for activation of different numbers of IR channels
%uses stimFreqs

%% Recall original analysis that assumed stimFreqs = [0 10 100 150 250 350 425];
% New analysis instead uses indices (1&2 cold, 3-6 warm (active), 7 hot (excluded))
responseByNum_indexBased = nan(numChannels, 5);
mnRespByNum_indexBased = nan(1, 5);
calcWarmHotColdInds_indexBased2  %pull indices of substim of different types

for channelNum = 1 : numChannels
    channelName = chanNamesMultiunit{channelNum};
    disp(['Calculating response to #channels activated for channel ' channelName])
    channelBaseMean = baselineMeans(channelNum);
    baselineCounts = allBaselineCounts(:, channelNum);
    [mnBaseline, stdBaselineCount] = mean_std(baselineCounts);
    substimCounts = allSubstimCounts(:, channelNum);

    % Two channels (front, back, left right)
    %irFront is: 1/2 warm 3/4 cold
    frontSpikeCounts_indexBased = substimCounts(irFrontInds); 
    [medianFront_indexBased, seMedianFront_indexBased] = median_sem(frontSpikeCounts_indexBased , 500);
    [mnFront_indexBased, semFront_indexBased] = mean_sem(frontSpikeCounts_indexBased );
    maxFront_indexBased = max(frontSpikeCounts_indexBased);
    mnPopvecFront_indexBased = nanmean(allPopvecs(irFrontInds,:));
    %Back: 1/2 cold 3/4 warm
    backSpikeCounts_indexBased = substimCounts(irBackInds);
    [medianBack_indexBased, seMedianBack_indexBased] = median_sem(backSpikeCounts_indexBased , 500);
    [mnBack_indexBased, semBack_indexBased] = mean_sem(backSpikeCounts_indexBased );
    maxBack_indexBased = max(backSpikeCounts_indexBased);
    mnPopvecBack_indexBased = nanmean(allPopvecs(irBackInds,:));
    %Left: 2/3 warm 1/4 cold
    leftSpikeCounts_indexBased = substimCounts(irLeftInds);
    [medianLeft_indexBased, seMedianLeft_indexBased] = median_sem(leftSpikeCounts_indexBased , 500);
    [mnLeft_indexBased, semLeft_indexBased] = mean_sem(leftSpikeCounts_indexBased );
    maxLeft_indexBased = max(leftSpikeCounts_indexBased);
    mnPopvecLeft_indexBased = nanmean(allPopvecs(irLeftInds,:));
    %Right: 1/4 warm 2/3 cold
    rightSpikeCounts_indexBased = substimCounts(irRightInds);
    [medianRight_indexBased, seMedianRight_indexBased] = median_sem(rightSpikeCounts_indexBased , 500);
    [mnRight_indexBased, semRight_indexBased] = mean_sem(rightSpikeCounts_indexBased );
    maxRight_indexBased = max(rightSpikeCounts_indexBased);
    mnPopvecRight_indexBased = nanmean(allPopvecs(irRightInds,:));
    %Combine all two on
    mnTwoOn_indexBased = (nanmean([mnRight_indexBased mnLeft_indexBased mnBack_indexBased mnFront_indexBased]) - mnBaseline)/stdBaselineCount;

    % Three channels (front and left etc)
    %Front Left: IR channels 1+2+3
    [medianFrontLeft_indexBased, seMedianFrontLeft_indexBased] = median_sem(substimCounts(irFrontLeftInds), 500);
    [mnFrontLeft_indexBased, semFrontLeft_indexBased] = mean_sem(substimCounts(irFrontLeftInds));
    %Front right: IR channels 1+2+4
    [medianFrontRight_indexBased, seMedianFrontRight_indexBased] = median_sem(substimCounts(irFrontRightInds), 500);
    [mnFrontRight_indexBased, semFrontRight_indexBased] = mean_sem(substimCounts(irFrontRightInds));
    %Back Left: Channels 2 3 4
    [medianBackLeft_indexBased, seMedianBackLeft_indexBased] = median_sem(substimCounts(irBackLeftInds), 500);
    [mnBackLeft_indexBased, semBackLeft_indexBased] = mean_sem(substimCounts(irBackLeftInds));
    %Back Right: Channels 3 4 1
    [medianBackRight_indexBased, seMedianBackRight_indexBased] = median_sem(substimCounts(irBackRightInds), 500);
    [mnBackRight_indexBased, semBackRight_indexBased] = mean_sem(substimCounts(irBackRightInds));
     mnThreeOn_indexBased = (nanmean([mnFrontLeft_indexBased mnFrontRight_indexBased ...
                                  mnBackLeft_indexBased mnBackRight_indexBased]) - mnBaseline)/stdBaselineCount;

    % Four channels: Channels 1-4
    [mnAll_indexBased, semAll_indexBased] = mean_sem(substimCounts(irAllInds));
    [medianAll_indexBased, seMedianAll_indexBased] = median_sem(substimCounts(irAllInds), 500);
    mnPopvecAll_indexBased = nanmean(allPopvecs(irAllInds,:));
    mnAllOn_indexBased = (mnAll_indexBased - mnBaseline)/stdBaselineCount;

    % Individual channels
    %Channel 1 only
    [mnOne_indexBased, semOne_indexBased] = mean_sem(substimCounts(ir1OnlyInds));
    [medianOne_indexBased, seMedianOne_indexBased] = median_sem(substimCounts(ir1OnlyInds), 500);
    mnPopvecOne_indexBased = nanmean(allPopvecs(ir1OnlyInds,:));
    %Channel 2 only
    [mnTwo_indexBased, semTwo_indexBased] = mean_sem(substimCounts(ir2OnlyInds));
    [medianTwo_indexBased, seMedianTwo_indexBased] = median_sem(substimCounts(ir2OnlyInds), 500);
    mnPopvecTwo_indexBased = nanmean(allPopvecs(ir2OnlyInds,:));
    %Channel 3 only (wow)
    [mnThree_indexBased, semThree_indexBased] = mean_sem(substimCounts(ir3OnlyInds));
    [medianThree_indexBased, seMedianThree_indexBased] = median_sem(substimCounts(ir3OnlyInds), 500);
    mnPopvecThree_indexBased = nanmean(allPopvecs(ir3OnlyInds,:));
    %Channel 4 only
    [mnFour_indexBased, semFour_indexBased] = mean_sem(substimCounts(ir4OnlyInds));
    [medianFour_indexBased, seMedianFour_indexBased] = median_sem(substimCounts(ir4OnlyInds), 500);
    mnPopvecFour_indexBased = nanmean(allPopvecs(ir4OnlyInds,:));
    %Combine all ones
    mnOneOn_indexBased = (nanmean([mnOne_indexBased mnTwo_indexBased ...
                                mnThree_indexBased mnFour_indexBased]) - mnBaseline)/stdBaselineCount;

    % Zero channels (All cold)
    [mnOff_indexBased, semOff_indexBased] = mean_sem(substimCounts(irOffInds));
    [medianOff_indexBased, seMedianOff_indexBased] = median_sem(substimCounts(irOffInds), 500);
     mnAllOff_indexBased = (mnOff_indexBased - mnBaseline)/stdBaselineCount;

    responseByNum_indexBased(channelNum,:) = [mnAllOff_indexBased mnOneOn_indexBased ...
                                              mnTwoOn_indexBased mnThreeOn_indexBased mnAllOn_indexBased];
end %channel num

[mnRespByNum_indexBased, semRespByNum] = mean_sem(responseByNum_indexBased);

if plotMinimal | plotFull
    figure;
    subplot(2,1,1)
    bar([0:4], proportionChannelsActivated, 'k');
    ylabel('Proportion of substim');
    title([animalName ' ' strrep(dateStr, '_', '/') ' #active']);

    subplot(2,1,2)
    bar([0 1 2 3 4], mnRespByNum_indexBased, 'k');hold on;
    plot_errorbars([0 1 2 3 4], mnRespByNum_indexBased, semRespByNum, 0, 'k', 1);
    xlabel('Number of stimulators active (index based)')
    ylabel('Z Score')

    if printOn
        print
    end
end  %if plotMinimal


%% Some error-checking, if desired
% error_check = 1;
% if error_check
%     if isequal(mnRespByNum_indexBased, mnRespByNum)
%         disp('Bingo! Index-based and value based match!')
%     end
%     scatter([0:4],mnRespByNum_indexBased - mnRespByNum, 'filled');
%     set(gca,'XTick',[0:4]);
%     axis([-.25 4.25 -1 1]);
%     grid on;
% end %if error_check


%% Save that shit
cd(datFolderNeuro)
save neuronalAnalysis
cd(datFolder)
saveNeuroString = ['save neurAnalysis_' animalName '_' dateStr ...
            ' mnRespByNum_indexBased rfCenters rfDiameters numSubstimAll '  ...
            ' peakZScores responseByNum_indexBased baselineMeans channelNumbers ' ...
            ' mnPsthZAll allPopvecs stimFreqs proportionChannelsActivated maxFreq'];
eval(saveNeuroString);


%% Discrimination analysis is in neuronalDiscriminationIr.m run separately in loadingSavingPre;p
%

%%
disp(['Done with neuronal analysis ' animalName ' ' dateStr])
