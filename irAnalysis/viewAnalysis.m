% Replots data after initial analysis
% Run after stimAnalysis and neuronalAnalysis
% Nicolelis Lab 4/7/17
%% Load analysis results and initialize variables
close all
cd(datFolder)

stimAnalysisFilename = ['stimAnalysisNeuro_' animalName '_' dateStr];  
load(stimAnalysisFilename);

loadNeuroString  = ['load neurAnalysis_' animalName '_' dateStr];
eval(loadNeuroString);

numTrials = numTrialsAll(sessNum);
numChannels = length(chanNamesMultiunit);
%% plot stimulus population
plot_stimulus_population(popvecSequences, animalName, maxFreq, dateStr);
%% plot psth of neural recordings
figure;
for channelNumber = 1 : numChannels
    plot_psth(channelNumber, chanNamesMultiunit,  mnPsthAll{channelNumber} );
end
suplabel([animalName ' ' strrep(dateStr, '_', '/') ' psth']);shg
%% Plot mean 3x3 population vector trajectory (this only needs doing once each session!)
plot_3x3_popvec(animalName, maxFreq, mnNPopvecs, numSubstimMin, dateStr);
%% Calculate one-neuron receptive fields
filterSigma = round((maxFreq*2)/85); %Not huge--10 for old;
filterWidth = round((maxFreq*2)/28); %Pretty big--30 for old 
figure;
for channelNum = 1:numChannels
    channelName = chanNamesMultiunit{channelNum};

    xVals = rfData{1}{channelNum};
    yVals = rfData{2}{channelNum};
    xMaxVal = rfData{3}{channelNum};
    yMaxVal = rfData{4}{channelNum};
    contourHeights = rfData{5}{channelNum};
    rfValsFiltered = rfData{6}{channelNum};
    
    plot_single_neuron_rf (channelNum, channelName, maxFreq,contourHeights,xMaxVal,xVals,yMaxVal,yVals,rfValsFiltered)

end
suplabel(['Filter size/sigma: ' num2str(filterWidth) ' ' num2str(filterSigma) ' ' animalName ' ' strrep(dateStr, '_', '/') ' contour']);shg
%% Peak z scores
mnPeakZ = zScoreData{1};
medianPeakZ = zScoreData{3};
peakZScores = zScoreData{5};
plot_peakZ(mnPeakZ,medianPeakZ,peakZScores,numChannels,animalName,dateStr)
%% Plot number of substim
plot_number_substim(proportionChannelsActivated, responseByNum_indexBased,animalName,dateStr)
disp('Finished re-plotting analysis');
%% Optional printing code
if printOn
    set(gcf,'Renderer','OpenGL')  
    print
    pause(0.1)
end


