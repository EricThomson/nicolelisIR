    
%%
%Central hub for loading and analyzing all stimuli and neuronal data for cohort of animals.
%
%To do:
%1. Note one place you need to automate in neuronalAnalysis6: 
%   responseInds = [83:88; 96:101; 109:114; 122:127; 135:140; 148:153];
% This is part of channelPsthPlotScript. For someone like Ford where stim
% duration is variable, this will be fubar.
% 
% This looks at indices to quantify responses, but when you change duration etc (Ford) this will change
% 5b. Similar for countDuration = 0.06; This will vary depending on stimu duration\
% But obviously also 
% 
% How did I even get these response Inds? 
%   Go into channelPsthPlotScript for nitty-gritty of psth bits until now you just have 
%   things in units of windows and times, not bin indices. You shoudl probably be calculating
%
%2. RF width and centers will need to be normalized to stimulation max, no?
% This can be done once you have integrated/collated all the data.
%
%3. Response versus stim magnitude is going to depend on max frequency, so be
%   careful this will not be comparable across animals or even sessions if you use
%   variable maxFreq! We haven't worked this out yet.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%BEHAVIOR/STIM ANALYSIS
%General settings
close all
clear
plotFull =  1;  %Plots microdetails including tick marks for timing of events
plotPopvecs = 0;  %plot stimulus population vectors (all trajectories)
plotPopvecCycles = 0;  %plot individual population vector cycles to look through
printOn = 0; %popvec actual/possible (grey/red), and trajectories (depending on plot settings)

%% Brian (V1 implant IR only right now)
animalName = 'brian';
numSessions = 2;
dateStrings = {'02_29_16', '03_02_16'};
numTrialsAll =[231 250];
maxFreqAll = [425 425];  %

%% Quagmire (V1 implant, IR task only)
animalName = 'quagmire';
numSessions = 2;
dateStrings = {'03_02_16', '03_14_16'};
numTrialsAll =[208, 80];
maxFreqAll = [425 425];  %


%% Analyze one session (mainly to make sure things are working)
sessNum = 1;
if sessNum <= numSessions
    dateStr = dateStrings{sessNum}; 
    numTrials = numTrialsAll(sessNum); 
    maxFreq = maxFreqAll(sessNum);
    datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] 
    stimulusAnalysis
else
    warning('stimulusAnalysis run -- Check session number')
end
display('done running behavior for one')

%% Run all stim/behavior analysis for an animal
for sessNum = 1: numSessions
    fprintf('\n\n')
    display(['Analyzing session ' num2str(sessNum) '/' num2str(numSessions) ' for ' animalName]);
    if sessNum <= numSessions
        dateStr = dateStrings{sessNum}; 
        numTrials = numTrialsAll(sessNum); 
        maxFreq = maxFreqAll(sessNum);
        datFolder = ['C:\Users\Nicolab\Desktop\irData\' animalName '\' dateStr] ;
        stimulusAnalysis
    else
        disp('Check session number')
    end
    disp('**************');
    pause(1);
    close all;
end
beep;
fprintf('\n****************\n')
display(['Done analyzing data for ' animalName '!!']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%
%% NEURONAL ANALYSIS 
%%
%NEURONAL preparation (common to all animals)
clear; close all; common_neuronal_prep

%% Animal-specific neuronal prep
%% Brian (V1 implant IR only right now)
cd('C:\Users\Nicolab\Desktop\irData\brian');
brian_neuronal_prep

%% Quagmire
cd('C:\Users\Nicolab\Desktop\irData\quagmire')
quagmire_neuronal_prep



%% Run one neuronal (after common_neuronal prep and then animalName_neuronal_prep)
sessNum = 1;
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

%% After running all neuronal animals, pool data:
%Go through all sessions for all animals and pull data such as rf diameters, rf centers, max z score dist, etc

%% First animal -- Brian-- this is different as you initialize and do first save of fullNeuronalData_s1v1
%%Neuronal preparation (common to all animals)
clear;close all; common_neuronal_prep;
brian_neuronal_prep;
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


