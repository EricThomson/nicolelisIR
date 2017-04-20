    
%%
%Central hub for loading and analyzing all stimuli and neuronal data for cohort of animals.
%
%%%%%%%%%%%%%%%%%%%%%%
%% NEURONAL ANALYSIS 
%%
%NEURONAL preparation (common to all animals)
clear; close all; common_neuronal_prep

%% Animal-specific neuronal prep
%% Amber
cd('C:\Users\Nicolab\Desktop\irData\amber');
amber_neuronal_prep

%% Run one neuronal (after common_neuronal prep and then animalName_neuronal_prep)
sessNum = 3;
if sessNum <= numSessions  
    dateStr = dateStrings{sessNum}; 
    numTrials = numTrialsAll(sessNum); 
    maxFreq = maxFreqAll(sessNum);
    artifactChannels = artifactChannelAll{sessNum};
    datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] ;
    datFolderNeuro = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr];
    neuronalAnalysis
else
    display('Out of range')
end

%% Run all neuronal (after neuronal prep)
for sessNum = 1 : numSessions
    fprintf('\n\n')
    display(['Analyzing session ' num2str(sessNum) '/' num2str(numSessions) ' for ' animalName]);
    
    if sessNum <= numSessions  
        dateStr = dateStrings{sessNum}; 
        numTrials = numTrialsAll(sessNum); 
        maxFreq = maxFreqAll(sessNum);
        artifactChannels = artifactChannelAll{sessNum};
        datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] ;
        datFolderNeuro = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr];
        neuronalAnalysis
    else
        display('Out of range')
    end
end % sessNum
beep;
fprintf('\n****************\n')
display(['Done analyzing data for ' animalName '!!']);
%close all %?

%FOR FIGURES;
%{
set(gcf,'Renderer','OpenGL')  
set(gcf,'Renderer', 'painters')
cd 'C:\Users\Eric\Dropbox\IR6\S1V1\s1_v1_paper\Figures\v1Recording'
print -depsc2 brianRFs
print -depsc2 brianPsths
export_fig brianPopvecTrajectories.pdf -pdf -transparent 

%}
%% Re-plot prior analysis
sessNum = 3;
if sessNum <= numSessions  
    dateStr = dateStrings{sessNum}; 
    numTrials = numTrialsAll(sessNum); 
    maxFreq = maxFreqAll(sessNum);
    artifactChannels = artifactChannelAll{sessNum};
    datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] ;
    datFolderNeuro = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr];
    viewAnalysis
else
    display('Out of range')
end

%% After running all neuronal animals, pool data:
%Go through all sessions for all animals and pull data such as rf diameters, rf centers, max z score dist, etc

%% First animal -- Amber-- this is different as you initialize and do first save of fullNeuronalData_thalamus
%%Neuronal preparation (common to all animals)
clear;close all; common_neuronal_prep;
amber_neuronal_prep;
combinedDataFolder = 'C:\Users\Nicolab\Desktop\irData\dataCombined';  
%%Recall individual variabels saved in neuronalAnalysisN
% mnRespByNum_indexBased x
% rfCenters x
% rfDiameters  x
% numSubstimAll x
% peakZScores x
% responseByNum_indexBased  x
% baselineMeans x
% channelNumbers x
% mnPsthZAll x
% allPopvecs x
% stimFreqs x
% proportionChannelsActivated 
% maxFreq

%Initialize variables
animalNamesAll = {};
sessionDatesAll = {};
mnRespByNumAll = [];
rfCentersAll = {};
rfDiametersAll = {};
numSubstimAllSessions = {};  %was AllTrials
peakZScoresAll = {};
respByNumAll = {};
baselineMeansAll = {};
channelNumbersAll = {};
mnPsthZAllSessions = {};
popvecsAll = {};
stimFreqsAll = []; %1x7 stimulus frequency array
proportionChannelsActiveAll = [];
maxFreqAll = [];

%First session is special needs to fill initial value of the variables.
sessNum = 1;
dateStr = dateStrings{sessNum}; 
datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] ;
cd(datFolder)
loadNeuroString  = ['load neurAnalysis_' animalName '_' dateStr];
eval(loadNeuroString);
disp(['Pooling data from ' animalName ' session ' num2str(sessNum) '/' num2str(numSessions)])
animalNamesAll = {animalName};
sessionDatesAll = {dateStr};
mnRespByNumAll = mnRespByNum_indexBased;
rfCentersAll = {rfCenters};
rfDiametersAll = {rfDiameters};
numSubstimAllSessions = {numSubstimAll};
peakZScoresAll = {peakZScores};
respByNumAll = {responseByNum_indexBased};
baselineMeansAll = {baselineMeans};
channelNumbersAll = {channelNumbers};
mnPsthZAllSessions = {mnPsthZAll};
popvecsAll = {allPopvecs};
stimFreqsAll = stimFreqs;
proportionChannelsActiveAll = proportionChannelsActivated;
maxFreqAll = maxFreq;

% Once first is done then can append data to existing arrays and structures.
for sessNum = 2: numSessions
    disp(['Pooling data from ' animalName ' session ' num2str(sessNum) '/' num2str(numSessions)])
    dateStr = dateStrings{sessNum}; 
    datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] ;
    cd(datFolder)
    loadNeuroString  = ['load neurAnalysis_' animalName '_' dateStr];
    eval(loadNeuroString);

    animalNamesAll{end+1} = animalName;
    sessionDatesAll{end+1} = dateStr;
    mnRespByNumAll = [mnRespByNumAll; mnRespByNum_indexBased]; %1x5
    rfCentersAll{end+1} = rfCenters;  %frCenters is numChannels x 2
    rfDiametersAll{end+1} = rfDiameters;  %rfDiameters is numChannels x 1
    numSubstimAllSessions{end+1} = numSubstimAll; %numSubstimAll is 1xnumTrials
    peakZScoresAll{end+1} = peakZScores;  %peakZScores is numChannels x 1
    respByNumAll{end+1} = responseByNum_indexBased; %numChan x 5
    baselineMeansAll{end+1} = baselineMeans; %baselineMeans is 1xnumChannels
    channelNumbersAll{end+1} = channelNumbers;  %NEW
    mnPsthZAllSessions{end+1} = mnPsthZAll;  %numChannels x 6 (6= num bins for psth)
    popvecsAll{end+1} = allPopvecs; %NEW
    stimFreqsAll = [stimFreqsAll; stimFreqs]; %NEW
    proportionChannelsActiveAll = [proportionChannelsActiveAll; proportionChannelsActivated]; %NEW
    maxFreqAll = [maxFreqAll; maxFreq]; %NEW
end
fprintf(['Done pooling data from ' animalName '\n'])

%Save that shit!
cd(combinedDataFolder)
try
    save fullNeuronalData_thalamus animalNamesAll sessionDatesAll mnRespByNumAll rfCentersAll rfDiametersAll ...
                               numSubstimAllSessions peakZScoresAll respByNumAll baselineMeansAll channelNumbersAll ...
                               mnPsthZAllSessions popvecsAll stimFreqsAll proportionChannelsActiveAll maxFreqAll
    fprintf(['\nBullseye! Data saved for ' animalName '! \n'])                       

catch exception 
   beep
   fprintf(['\n\nYa'' blew it! Something went wrong with saving data for ' animalName '! \n']);
   error(exception.message);
end
                       

%% Second, need to do other animals after initial animal is done (peter, quagmire, cleveland)
clear;close all
animals_to_run = {'quagmire'};  %would usually have multiple animals
combinedDataFolder = 'C:\Users\Nicolab\Desktop\irData\dataCombined';  
numberOfAnimals = length(animals_to_run);

for animalNum = 1: numberOfAnimals
    common_neuronal_prep  %includes clear and close all
    animal_specific_string = [animals_to_run{animalNum} '_neuronal_prep'];  %e.g., peter_neuronal_prep
    eval(animal_specific_string);

    cd(combinedDataFolder)
    load fullNeuronalData_s1v1
    for sessNum = 1: numSessions
        %disp(['Pooling data from ' animalName ' session ' num2str(sessNum) '/' num2str(numSessions)])
        dateStr = dateStrings{sessNum}; 
        datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] ;
        cd(datFolder)
        loadNeuroString  = ['load neurAnalysis_' animalName '_' dateStr];
        eval(loadNeuroString);

        animalNamesAll{end+1} = animalName;
        sessionDatesAll{end+1} = dateStr;
        mnRespByNumAll = [mnRespByNumAll; mnRespByNum_indexBased]; %1x5
        rfCentersAll{end+1} = rfCenters;  %frCenters is numChannels x 2
        rfDiametersAll{end+1} = rfDiameters;  %rfDiameters is numChannels x 1
        numSubstimAllSessions{end+1} = numSubstimAll; %numSubstimAll is 1xnumTrials
        peakZScoresAll{end+1} = peakZScores;  %peakZScores is numChannels x 1
        respByNumAll{end+1} = responseByNum_indexBased; %numChan x 5
        baselineMeansAll{end+1} = baselineMeans; %baselineMeans is 1xnumChannels
        channelNumbersAll{end+1} = channelNumbers;  %NEW
        mnPsthZAllSessions{end+1} = mnPsthZAll;  %numChannels x 6 (6= num bins for psth)
        popvecsAll{end+1} = allPopvecs; %NEW
        stimFreqsAll = [stimFreqsAll; stimFreqs]; %NEW
        proportionChannelsActiveAll = [proportionChannelsActiveAll; proportionChannelsActivated]; %NEW
        maxFreqAll = [maxFreqAll; maxFreq]; %NEW
    end
    %fprintf(['Done pooling data from ' animalName '\n'])
    %Save that shit!
    cd(combinedDataFolder)
    try
        save fullNeuronalData_s1v1 animalNamesAll sessionDatesAll mnRespByNumAll rfCentersAll rfDiametersAll ...
                                   numSubstimAllSessions peakZScoresAll respByNumAll baselineMeansAll channelNumbersAll ...
                                   mnPsthZAllSessions popvecsAll stimFreqsAll proportionChannelsActiveAll maxFreqAll
        fprintf(['\nBullseye! Data saved for ' animalName '! \n'])                       
    catch exception 
       beep
       fprintf(['\n\nYa'' blew it! Something went wrong with saving data for ' animalName '! \n']);
       error(exception.message);
    end

end %animalNumber
%%


