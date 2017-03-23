%A. Analysis of data from awake animal
% 1. Loading and saving data: loadingSavingPrep 
%       This is the central hub for controlling loading, saving, and analysis of all data.
%       Use it to loop through stimulus analysis and neuronal analysis for each animal.
%
% 2.Stimulus analysis (stimulusAnalysis.m)
% 	-In Nex save all events in stimEvents.mat (irData/animal/date)      
% 	-Save session_info from behavioral session as animal_name_date (e.g., pinky_11_14_14)
% 	-Analyze  (use stimulusAnalysis.m) and it will save stimResults_name_date.mat  
%       It will also save stimAnalysisNeuro_name_date.mat (minimal subset of 
%         variables just for  later neuronal analysis) -- datFolder
%
%3. Neuronal (neuronalAnalysis.m)
%   Awake during task:
%   First, in OLS sort as many channels as you think you need to get good artifact vote. 
%       Usually four is enough (actually, use odd #), but if there is a lot of false positives, then more, and might want odd #. 
%       Be sure to click to turn on clock events in OLS (event 012).
%   In Nex save data with stim artifact sorted as units in neuronalEvents.mat (in Desktop)
%       C:\Users\Eric\Desktop\stimRecordData\name\date (or neuronalEventsSingle.mat for single units)
%   Use neuronalAnalysisN to analyze the data in awake animal
%       For discrimination analysis: neuronalDisciminationIr,
%       run in loadingSavingPrep2
%
%   Awake air puff:
%       airPuffAwakeAnalysis2 or single
%
%   Combined neuronal analysis: neuronalAllAnalysis (receptive fields, discrimination, etc).
%                               neuronalAllAnalysis_s1v1


%B. Other stuff
%   1. To analyze reliability of timing of events sent to plexon, see timeTesting.m (right now on desktop in lab)
%Previously you had used:
%   Center of mass analysis: neuronalAnalysisCM.m and neuronalAnalysisCycleCM (normalized) [before you settled on pop vec]
%   Actual stim values (neuronalAnalysisMags)
%
%
%Pull "sham" stim for control animals (pinky/brain)
%   Pulled data pool for control stimuli : stim values, and trial periods (analysis\poolPull.m)
%   learningParameters.mat and learnedParameters.mat
%
%Analysis of anesthetized basis set:  anesthetizedWorking