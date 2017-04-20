%neuronalAllAnalysis_s1v1
clear;close all
combinedDataFolder = 'C:\Users\Eric\Dropbox\IR6\S1V1\s1_v1_paper\data\dataCombined';  
cd(combinedDataFolder)
load fullNeuronalData_s1v1
numAnimals = length(unique(animalNamesAll));
%Loads:
%animalNamesAll         %1xnumSessions cell Array
%sessionDatesAll        %1xnumSessions cell Array
%mnRespByNumAll         %numSessions x 5 array response by number active
%respByNumAll           %1xnumSessions cell array contains numChannels x 5 
%rfCentersAll           %1xnumSessions cell Array
%rfDiametersAll         %1xnumSessions cell Array
%numSubstimAllSessions  %1xnumSessions cell Array (each is 1 x numTrials)
%peakZScoresAll         %1xnumSessions cell Array
%baselineMeansAll       %1xnumSessions cell Array
%channelNumbersAll      %1xnumSessions cell Array
%mnPsthZAllSessions     %1xnumSessions cell Array
%popvecsAll             %1xnumSessions cell Array
%stimFreqsAll                   %numSessions x 7 array of frequencies stimulated
%proportionChannelsActiveAll    %numSessions x 5 (proportion IR channels activated warm)
%maxFreqAll                     %numSessions x 1 maximum frequency of stimulation

numSessions = length(respByNumAll);
respByNumCombined = cell2mat(respByNumAll');
numChannelsTotal = size(respByNumCombined,1);
%all rfs
rfCentersCombined = cell2mat(rfCentersAll');
[mnRfAll, semRfAll] = mean_sem(rfCentersCombined);


%% 
%PSTH response z mean
psthZCombined = cell2mat(mnPsthZAllSessions');
[mnZPsth, semZPsth] = mean_sem(psthZCombined);
[p, table, stats] = anova1(psthZCombined);  %p = 1.5 e -15
multcompare(stats)
%% plot psth response
bar(mnZPsth, 'k');hold on
plot_errorbars([1:6], mnZPsth, semZPsth, 0, 'k', 1);
axis([0.5 6.5 0 57])
set(gca,'YTick',[0 10 20 30 40 50])
xlabel('Stimulus number')
ylabel('Response Magnitude (z)')
shg
%print -depsc2 magVsSubstim


%% mean/sem number of substim over all trials
numSubstimCombined = cell2mat(numSubstimAllSessions)';
[mnSubstimNum, semSubstimNum] = mean_sem(numSubstimCombined);  %7.4 +/- 0.05
%
%% plot histogram of all substim numbers 
histSubstimNum = hist(numSubstimCombined, [0:max(numSubstimCombined)])
bar([0:max(numSubstimCombined)], histSubstimNum/sum(histSubstimNum), 'k');hold on;
axis([0 35 0 0.34]);
xlabel('Number of substimuli')
ylabel('Frequency')
%print -depsc2 numSubstimHistogram


%********************
%% RF and stimvec analysis

%% Mask for histogram: there will be zeros outside diamond: these will strongly bias it
stimFreqs = stimFreqsAll(1,:); %[0 10 100 150 250 350 425]; 
maxFreq = maxFreqAll(1);
detectorCoords = [1 -1 -1 1;-1 -1 1 1];  %2x4, channel i has column i 
allStimCombos = combvec(stimFreqs, stimFreqs, stimFreqs, stimFreqs)';
popVecspace = (detectorCoords * allStimCombos')';


% limVal = (maxFreq*2)+25;
% rfSpacing = round((maxFreq*2)/170);  %Small spacing--5 for old
% filterSigma = round((maxFreq*2)/85); %Not huge--10 for old;
% filterWidth = round((maxFreq*2)/28); %Pretty big--30 for old 
% imageFilter=fspecial('gaussian', filterWidth, filterSigma);


binWidthHists = 85;  %old
binCenters  = [-850: binWidthHists:850];  %old

histPopvecSpace = hist3(popVecspace, {binCenters, binCenters} ); 
histMaskInds = find(histPopvecSpace == 0);
histPopvecSpace(histMaskInds) = NaN;


%% All population vectors actually presented on top of all possible
popvecsCombined = cell2mat(popvecsAll');
uniquePopvecs = unique(popvecsCombined, 'rows');
%%
%Plot all actual popvecs versus all possible
%possible
scatter(popVecspace(:,1),popVecspace(:,2), 40, 'filled', 'CData',[.7 .7 .7]);   hold on
%actual
scatter(uniquePopvecs(:, 1), uniquePopvecs(:, 2), 30, 'r','filled');hold on
grid on; axis equal;
plot([-850 0 850 0 -850], [0 850 0 -850 0],'k','LineWidth', 1); %the square   
axis([-850 850 -850 850]);
set(gca,'XTick',[-850:170:850], 'YTick', [-850: 170: 850])
%
%print -depsc2 everyPopvec

    
%% Population vector statistics
%Plot all popvecs that actually appeared
scatter(uniquePopvecs(:, 1), uniquePopvecs(:, 2), 30,'filled', 'k');hold on
plot([-850 0 850 0 -850], [0 850 0 -850 0],'k','LineWidth', 1); %the square
grid on; axis equal;
axis([-860 860 -860 860]);
set(gca,'XTick',[-850:170:850], 'YTick', [-850: 170: 850]);
xlabel('Mediolateral')
ylabel('Anteroposterior')
title('All population vectors that actually occurred');shg
gridfix([.7 .7 .7])
shg
%print -depsc2 actualPopvecsAll


%% 3d histogram of actual popvecs
% compare actual versus centers
%CHANGE METHOD: DO INTERPOLATION AND SMOOTHING AND COMPARE "OVERLAP"

binWidthHists = 85;
binCenters  = [-850: binWidthHists:850];  %85
[histActualStim, histActualEdges] = hist3(popvecsCombined, {binCenters, binCenters} ); 
hist3(popvecsCombined, {binCenters, binCenters} ); xlabel('x');ylabel('y');
histActualStim(histMaskInds) = NaN;  %histActualStimDiag'
histActualStim(6,6) = NaN;
%%

histRFCenters = hist3(rfCentersCombined, 'Edges', {binCenters, binCenters} ); 
hist3(rfCentersCombined, {binCenters, binCenters} ); xlabel('x');ylabel('y');
histRFCenters(histMaskInds) = NaN ;
histRFCenters(6,6) = NaN;
xlabel('mediolateral');ylabel('anteroposterior')
%%
% %Significant correlation? OLD way with raw data and no interpolation
% rfHistMasked = histRFCenters(:);
% rfHistMasked(isnan(rfHistMasked)) =[];
% popvecHistMasked = histActualStim(:);
% popvecHistMasked(isnan(popvecHistMasked)) =[]
% rfActualCorr  = corr2(popvecHistMasked, rfHistMasked);  %-0.03
% [slope,corr_co,p_value, y_intercept]=regression_fit(popvecHistMasked(:), rfHistMasked(:),1); ; %p=0.67


%% Plot antero and mediolat frequency histograms (popvecsCombined)


%% anteroposterior actual popvecs
figure;
binCenters  = [-850:85:850];  %85
pvY = popvecsCombined(:,2);
nY = hist(pvY, binCenters);
freqY = nY/sum(nY);
bar(binCenters, freqY, 1,'k');hold on
axis([-850 850 0 0.25])
xlabel('Anteroposterior position')
ylabel('Frequency of Stimulus')
title('Actual popvecs')
set(gca, 'YTick',[0 0.05 .1 .15 0.2 0.25])
set(gca, 'XTick', [-850:170:850])
shg;
%print -depsc2 yHistogramPopvecs


%compare mean of total popvecs and mean of rf centers


%-18, 215


%143/203
%
%[h, plocations] = ttest2(popvecsCombined(:,2), rfCentersCombined(:,2));
%[~, pMeanLocations] = manova1([popvecsCombined; rfCentersCombined], ...
%    [zeros(size(popvecsCombined,1),1); ones(size(rfCentersCombined,1),1)]); %0.1154
%print -depsc2 rfCentersAllPlot

%% mediolateral actual popvecs
figure;
binCenters  = [-850:85:850];  %85
pvX = popvecsCombined(:,1);
nX = hist(pvX, binCenters);
freqX = nX/sum(nX);
bar(binCenters, freqX, 1,'k');hold on
axis([-850 850 0 0.32])
xlabel('Mediolateral position')
ylabel('Frequency of Stimulus')
title('Actual popvecs')
% plot([-850 -850], [0 0.25],'k:')
% plot([850 850], [0 0.25],'k:')
set(gca, 'YTick',[0 0.05 .1 .15 0.2 .25 0.3])
set(gca, 'XTick', [-850:170:850])
shg;
%print -depsc2 xHistogramPopvecs







%% RFs 
%% Plot all rf centers
%rfCentersCombined = cell2mat(rfCentersAll');
%Plot all rf centers
scatter(rfCentersCombined(:,1), rfCentersCombined(:,2), 30, 'filled', 'k'); hold on; %'CData', [.5 .5 .5]); hold on;
plot([-850 0 850 0 -850], [0 850 0 -850 0],'k','LineWidth', 1); %the square
grid on; axis equal;
axis([-850 850 -850 850]);
set(gca,'XTick',[-850:170:850], 'YTick', [-850: 170: 850]);
title('RF centers')
gridfix([.7 .7 .7])
shg
%scatterhist?
%scatterhist(rfCentersCombined(:,1), rfCentersCombined(:,2))

%print -depsc2 rfCentersFull


%% Plot histogram of RF centers along anteroposterior
%for early versus late learners
%Either subplot, or stack bars
binCenters  = [-850:85:850];  %85
rfY = rfCentersCombined(:,2);
nY = hist(rfY, binCenters);
freqY = nY/sum(nY);
bar(binCenters, freqY,  1,'k');hold on
axis([-850 850 0 0.2])
xlabel('Anteroposterior position')
ylabel('Frequency')
title('RF centers')
set(gca, 'YTick',[0 0.05 .1 .15])
set(gca, 'XTick', [-850:170:850])
shg;
%print -depsc2 rfCenterYHistogram

%% Mediolateral RF center frequency
binCenters  = [-850:85:850];  %85
rfX = rfCentersCombined(:,1);
nX = hist(rfX, binCenters);
freqX = nX/sum(nX);
bar(binCenters, freqX, 1, 'k');hold on
axis([-850 850 0 0.14])
xlabel('Mediolateral Position')
ylabel('Frequency')
title('RF centers')
set(gca, 'YTick',[0 0.05 .12])
set(gca, 'XTick', [-850:170:850])
shg;
%print -depsc2 rfCenterXHistogram



%% Back to 2d correlation, but on interpolated density function

%% Create filtered version of each distribution
%popvecsCombined vs rfCombined
plotOn = 1;
plotOff = 0;

%% Get mean and 95% confidence interval of each
numBoots = 10000;
alphaLevel = 0.05;
% [ellipseValues,eigenVectors,average,majorAxis,minorAxis] = gaussianConfidenceBounds(data,confLevel)
[mnPopvec, semPopvec] = mean_sem(popvecsCombined)
[confidenceEllipse, ~, ~, ~, ~] = ...
    bootstrap2dConfidence(popvecsCombined, numBoots, alphaLevel, plotOn);

[mnRf, semRf] = mean_sem(rfCentersCombined)
[confidenceEllipseRf, ~, ~, ~, ~] = ...
    bootstrap2dConfidence(rfCentersCombined, numBoots, alphaLevel, plotOn)

%%****************************************
%% Next step: Look at correlations among distributions
densitySpacing = 10;  %there are 360 xVals/yVals
limVal = 860;
filterWidth = 30; %there are 360 xVals/yVals, about 100 space between jumps
filterSigma = 10;
imageFilter=fspecial('gaussian', filterWidth, filterSigma);
binCenters = [-limVal + densitySpacing/2 :densitySpacing: limVal - densitySpacing/2];



%% first get all four density functions and their maksed versions
%Setting up boundary bits
densitySpacing = 86;  %there are 360 xVals/yVals
binCenters = [-limVal + densitySpacing/2 :densitySpacing: limVal - densitySpacing/2];

%% Mask for histogram: there will be zeros outside diamond: these would bias it, so turn them to NaNs
allStimCombos=combvec(stimFreqs, stimFreqs, stimFreqs, stimFreqs)';
popVecspace = (detectorCoords * allStimCombos')';
histPopvecSpace = hist3(popVecspace, {binCenters, binCenters} ); 
histMaskInds = find(histPopvecSpace == 0);
histPopvecSpace(histMaskInds) = NaN;


%%
filterWidth = 7; %there are 360 xVals/yVals, about 100 space between jumps
filterSigma = 1.5;
imageFilter=fspecial('gaussian', filterWidth, filterSigma);
%figure;surf(imageFilter)

%%
% Popvecs 
[histCountsPopvecs, histBins] = hist3(popvecsCombined,{binCenters, binCenters});
histCountsPopvecs = histCountsPopvecs';
maxValPop = max(histCountsPopvecs(:));
[xIndMax, yIndMax] = find(histCountsPopvecs == maxValPop)
histCountsPopvecs(xIndMax, yIndMax) = NaN;

%Mask out the out of range bits
histCountsPopvecs(histMaskInds) = NaN;
%Filter it
popvecsFiltered = imFilterNan(histCountsPopvecs, imageFilter, plotOn);
suptitle('V1 popvecs')

%% Rfs 
[histCountsRfs, histBins] = hist3(rfCentersCombined,{binCenters, binCenters});
histCountsRfs = histCountsRfs';
%Mask out the out of range bits
histCountsRfs(histMaskInds) = NaN;
histCountsRfs(xIndMax, yIndMax) = NaN;
%Filter it
rfsFiltered = imFilterNan(histCountsRfs, imageFilter, plotOn);
subtitle('V1 rfs')

%% 2d  popvecs linearize and filter out the NaNs
popvecsLinear = popvecsFiltered(:);
popvecsLinear(isnan(popvecsLinear)) =[];


%% 2d corr popvecs and RF centers
rfsLinear = rfsFiltered(:);
rfsLinear(isnan(rfsLinear)) =[];

%% Do the correlation
%Note for S1 data this was: %p=7e-12, rho=.44
[slope, corr_co, p_value, y_intercept]=regression_fit(popvecsLinear, rfsLinear, plotOn); %p=0.7; corr_co (rho) = -.03 
title('V1 RF and STIM correlations')
xlabel('pop vec');ylabel('rfs')
%p = 2.2 e-12 rho = 0.04? 
%p= 0.7, rho = 0.03???


save s1v1Analysis_Working1

%{

clear;close all
combinedDataFolder = 'C:\Users\Eric\Dropbox\IR6\S1V1\s1_v1_paper\data\dataCombined';  
cd(combinedDataFolder)
load s1v1Analysis_Working1

%}


%% Below stuff is extra, not used



%%
%%
%%
%%%%%FOLLOWING NOT CONVERTED TO V1 YET!!! IT IS IMPORTED FROM S1
%% NEW PLOT Plot rfs and stim together
%%Population vectors
subplot(1,2,1)
contourf(binCenters, binCenters,popvecsFiltered );  set(gca,'YDir', 'normal');hold on
axis equal;
axis tight
title('population vectors V1');
axis([-860 860 -860 860]);
set(gca,'XTick',[-850:425:850], 'YTick', [-850: 425: 850]);
%%RFs
subplot(1,2,2)
contourf(binCenters, binCenters, rfsFiltered); set(gca,'YDir', 'normal');hold on
axis equal;
axis tight
axis([-860 860 -860 860]);
set(gca,'XTick',[-850:425:850], 'YTick', [-850: 425: 850]);
title('rf centers V1')
suptitle('V1 (p = 0.7, rho = -0.03)')

%{
print -depsc2 rf_popvec_compare_v1
%}


%% Rose plot for RF angles--do this with 2016b using polarhistogram which lets you enter the histogram edges...
[thetaRFs, rhoRFs] = cart2pol(rfCentersCombined(:,1), rfCentersCombined(:,2));
figure;rose(thetaRFs, 16) %, roseBinCenters)
roseCardinalCentersDeg = [-180:45:180]
roseCardinalCentersRad = deg2rad(roseCardinalCentersDeg);
n = hist(thetaRFs, roseCardinalCentersRad)
bar(roseCardinalCentersDeg, n);shg
%print -depsc2 rfRose_v1


%% RF diameters...
rfDiametersCombined = cell2mat(rfDiametersAll');
[mnRfDiametersV1, semRfDiametersV1] = mean_sem(rfDiametersCombined); %279.6+/-9.06

%% get s1 diamters
load rfDiametersS1
[mnRfDiametersS1, semRfDiametersS1] = mean_sem(rfDiametersS1);

[hRfDiam, pRfDiam] = ttest2(rfDiametersCombined, rfDiametersS1);  %p = 0.42
%% Plot
bar([0 1], [mnRfDiametersV1 mnRfDiametersS1],'k');hold on;
axis([-.5 1.5 0 320])
plot_errorbars([0 1], [mnRfDiametersV1 mnRfDiametersS1],[semRfDiametersV1 semRfDiametersS1],1, 'k', 2);
set(gca,'XTickLabel',{'S1','V1'})
ylabel('RF Diameter');
title('RF diameter in S1 vs V1 (p = 0.42)')

print -depsc2 rfDiameterCompare

%% 
save s1v1Analysis_Working1


%% Interpolate--not sure what this is
% interpSpacing = 10;  %there are 360 xVals/yVals
% % densitySpacing = 86;  %there are 360 xVals/yVals
% % binCenters = [-limVal + densitySpacing/2 :densitySpacing: limVal - densitySpacing/2];
% limVal = 860;
% binCentersInterp = [-limVal + interpSpacing/2 : interpSpacing: limVal - interpSpacing/2];
% [xGridInterp, yGridInterp] = meshgrid(binCentersInterp, binCentersInterp);
% interpDensity = TriScatteredInterp(binCentersInterp', binCentersInterp', histCountsPopvecs);
% %scatter3(uniqueXy(:,1), uniqueXy(:,2), freqHist);
% densityVals = densityFunction(xGrid,yGrid);

%% FIlter it
filterWidth = 5; %there are 360 xVals/yVals, about 100 space between jumps
filterSigma = 1;
imageFilter=fspecial('gaussian', filterWidth, filterSigma);surf(imageFilter)
densityPopvecsFiltered = imFilterNan(histCountsPopvecs, imageFilter, plotOn);


%%
subplot(1,2,1);hist3(popvecsCombined,{binCenters, binCenters});
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode', 'auto');
xlabel('x');ylabel('y')
subplot(1,2,2); contourf(binCenters, binCenters, histCountsPopvecs);shg
surf(binCenters, binCenters, histCountsPopvecs)
%%
%Do this, then interpolate...then filter...then get correlation

%Popvecs  using old way that sort of sucked
[densityPopvecs, ~, ~, ~, ~] = density2dCalculate(popvecsCombined, densitySpacing, limVal, plotOn);
densityPopvecsFiltered = imFilterNan(densityPopvecs, imageFilter, plotOn);
normPopvecNanMask = isnan(densityPopvecsFiltered);


%RF 
[densityRfs, densityFunction, countHist, xGrid, yGrid] = ...
    density2dCalculate(rfCentersCombined, densitySpacing, limVal, plotOn);
densityRfsFiltered = imFilterNan(densityRfs, imageFilter, plotOn);
[histCounts, histBins] = hist3(rfCentersCombined,{[-860:86:860], [-860:86:860]});


%%  population vectors
% popvecIndsZero = find(ismember(popvecsCombined, [0 0], 'rows'));
% totalPopvecsNozero = popvecsCombined;
% totalPopvecsNozero(popvecIndsZero,:) = [];
plotOn = 1;
plotOff = 0;
[densityPopvecs, densityFunction, countHist, xGrid, yGrid] = density2dCalculate(popvecsCombined, densitySpacing, limVal, plotOn);
densityPopvecsFiltered = imFilterNan(densityPopvecs, imageFilter, plotOn);


%% Build uber-mask (or the NaNs to pull all NaNs) and then get 2d Corr.
normNanMask = isnan(densityPopvecsFiltered);
fullNanMask = normNanMask; size(fullNanMask)

%% plot masks
subplot(1,3,1); contourf(normNanMask);axis equal; axis tight;grid on;
subplot(1,3,3); contourf(fullNanMask);axis equal;axis tight;shg;grid on

%%
masked = densityPopvecsFiltered;
masked(fullNanMask) = NaN;

%% plot
contourf(masked);;axis equal; axis tight;grid on;


%% 2d correlation for population vectors
%pull out non-nan guys:
%rfActualCorr  = corr2(popvecHistMasked, rfHistMasked);  %-0.05
%[slope,corr_co,p_value, y_intercept]=regression_fit(popvecHistMasked(:), rfHistMasked(:),1); ; %p=0.52
normPopvecs = masked(:);
normPopvecs(isnan(normPopvecs)) = [];
[slope, corr_coDiag, p_valueDiag, y_intercept]=regression_fit(diagPopvecs, normPopvecs, plotOn); ; %p=0 (hyp that; rho = 0.79
[r,p]=corrcoef(diagPopvecs, normPopvecs) 
%Each p-value is the probability of getting a correlation as large as the observed value by 
%random chance, when the true correlation is zero. If the correlation is
%significant, then the p value will be low.
corr_co = r(1,2);
p_corr = p(1,2);

%%FOLLOWING USEFUL?
%% compare to correlation coefficnet with scrambled bits
randTestR = zeros(1,1000);
for runNum =  1: 1000
    [rscram, pscram] = corrcoef(diagPopvecs, normPopvecs(randperm(length(normPopvecs))));
    randTestR(runNum) = rscram(1,2);
end
%%


%%%2dcorr  for rfs:  rfCentersCombined and rfCentersDiagonal
%% RF centers  and popvec distribution...
%already have  popvec distribution: masked
figure;scatter(rfCentersCombined(:,1), rfCentersCombined(:,2))

densitySpacing = 10;
limVal = 900;
[densityRfs, densityFunction, countHist, xGrid, yGrid] = ...
    density2dCalculate(rfCentersCombined, densitySpacing, limVal, plotOn);
densityRfsFiltered = imFilterNan(densityRfs, imageFilter, plotOn);
hist3(rfCentersCombined, {binCenters, binCenters} ); 

gridx1 = 0:.05:1;
gridx2 = 5:.1:10;
X = [0+.5*rand(20,1) 5+2.5*rand(20,1);
    .75+.25*rand(10,1) 8.75+1.25*rand(10,1)];
ksdensity2d(X,gridx1,gridx2);

densitySpacing = 10;
bandOne = 75;
bandWidth = [bandOne bandOne];
binCenters = [-limVal + densitySpacing/2 :densitySpacing: limVal - densitySpacing/2];
densityRfal = ksdensity2d(rfCentersCombined,binCenters,binCenters, bandWidth);
densityRfalMasked = densityRfal';
densityRfalMasked(normNanMask) = NaN;
contourf(binCenters, binCenters,densityRfalMasked);

densityPopvecs = ksdensity2d(popvecsCombined(1:1000,:),binCenters,binCenters, bandWidth);


%% TO DO
%Need to use your algorithm for large ones, and can use this for smaller
%--Do for norm popvec/RF center
%--Do for diag popvec/RF center
%
%%
subplot(1,2,1); contourf(binCenters, binCenters,densityPopvecsFiltered);axis equal; axis tight; grid on
subplot(1,2,2); contourf(binCenters, binCenters,densityRfalMasked);axis equal; axis tight; grid on
%%

contourf(densityRfal');
surf(densityPopvecsFiltered


%% Plot
[xGrid, yGrid] = meshgrid(binCenters, binCenters);
h=surf(xGrid, yGrid, f');
%set(gca,'YDir','reverse', 'XDir', 'reverse')
set(h, 'edgecolor','none') %remove lines
xlabel('x')
ylabel('y');shg
%%



%%
% KSDENSITY2D Compute kernel density estimate in 2D.
% F = KSDENSITY2D(X,GRIDX,GRIDX2,BW) computes a nonparametric estimate of
% the probability density function of the sample in the N-by-2 matrix X.
% F is the vector of density values evaluated at the points in the grid
% defined by the vectors GRIDX1 and GRIDX2. The estimate is based on a
%  kernel function, using a window parameter (bandwidth) that is a
% function of the number of points in X.



[densityRfDiag, densityFunctionDiag, countHist, xGrid, yGrid] = density2dCalculate(popvecsFullDiag, densitySpacing, limVal, plotOn);
densityPopvecsFilteredDiag = imFilterNan(densityPopvecsDiag, imageFilter, plotOn);



%% Compare the two histograms...
%Popvec
numal = size(popvecsCombined,1);
numDiag = size(popvecsFullDiag, 1);
combinedPopvecs = [popvecsFullDiag; popvecsCombined];
combinedPopvecsGroups = [ones(size(popvecsFullDiag,1),1); 2*ones(size(popvecsCombined,1),1)];
[d, pvals] = manova1(combinedPopvecs, combinedPopvecsGroups)



%% anteroposterior diagonal
popvecYDiag = popvecsFullDiag(:,2);
nYDiag = hist(popvecYDiag, binCenters);
freqYDiag = nYDiag/sum(nYDiag);
bar(binCenters, freqYDiag,  1,'k');hold on
axis([-850 850  0 0.2])
xlabel('Anteroposterior position')
ylabel('Frequency')
title('Population vector actual/diagonal')
% plot([-850 -850], [0 0.27],'k:')
% plot([850 850], [0 0.27],'k:')
set(gca, 'YTick',[0 .1 .2])
set(gca, 'XTick', [-850:170:850])
shg;
%print -depsc2 popvecYHistogramDiag

%% mediolateral diagonal
popvecXDiag = popvecsFullDiag(:,1);
nXDiag = hist(popvecXDiag, binCenters);
freqXDiag = nXDiag/sum(nXDiag);
bar(binCenters, freqXDiag,  1,'k');hold on
axis([-850 850 0 0.3])
xlabel('Mediolateral position')
ylabel('Frequency')
title('Popvec actual/diagonal')
% plot([-850 -850], [0 0.35],'k:')
% plot([850 850], [0 0.35],'k:')
set(gca, 'YTick',[0 .1 .2 0.3])
set(gca, 'XTick', [-850:170:850])
shg;
%print -depsc2 popvecXHistogramDiag



%% Get quadrant of each recording cahnnel
load channelNumbersAll
electrodeLocationMap = {[19, 20, 23, 24], [9, 10, 13, 14], [27, 28, 31, 32], [1, 2, 5, 6]};
uniqueChannelNumbers = unique(channelNumbersAll);
numUniqueChannels = length(uniqueChannelNumbers);
channelLocationsAll = zeros(numChannelsTotal, 1);
for chanNumberInd = 1: numUniqueChannels
    channelNumber = uniqueChannelNumbers(chanNumberInd);
    channelNumberInds = find(channelNumbersAll == channelNumber);
    channelLocation = find(cellfun(@(x) sum(find(x == channelNumber)), electrodeLocationMap));
    channelLocationsAll(channelNumberInds) = channelLocation;
end
meanResponse = grpstats(rfCentersCombined,channelLocationsAll)
[d, p, stats] = manova1(rfCentersCombined, channelLocationsAll) ;
manovacluster(stats)

%% Plot four RF centers
mnOverallRF = mean(rfCentersCombined);
scatter(mnOverallRF(1), mnOverallRF(2),'k','filled');hold on;
scatter(means(1,1), means(1,2),'r', 'filled');hold on
scatter(means(2,1), means(2,2),'g', 'filled')
scatter(means(3,1), means(3,2),'b', 'filled')
scatter(means(4,1), means(4,2),'c', 'filled')
legend({'All', 'AL (1)', 'AR (2)', 'PL (3)', 'PR (4)'})
plot([-850 0 850 0 -850], [0 850 0 -850 0],'k','LineWidth', 1); %the square
axis([-860 860 -860 860]);
set(gca,'XTick',[-850:170:850], 'YTick', [-850: 170: 850]);
title('Mean RF center versus recording quadrant')
grid on; axis equal;
gridfix([.7 .7 .7])

%% PLot all
scatter(rfCentersCombined(find(channelLocationsAll == 1), 1), rfCentersCombined(find(channelLocationsAll == 1), 2), 'r');hold on
scatter(rfCentersCombined(find(channelLocationsAll == 2), 1), rfCentersCombined(find(channelLocationsAll == 2), 2), 'g');hold on
scatter(rfCentersCombined(find(channelLocationsAll == 3), 1), rfCentersCombined(find(channelLocationsAll == 3), 2), 'b');hold on
scatter(rfCentersCombined(find(channelLocationsAll == 4), 1), rfCentersCombined(find(channelLocationsAll == 4), 2), 'c');hold on


%% Build a classifer to predict electrode location, given RF center
%Horrible...
plotOn = 1;
classifier_type = 'knn';
[pc, ~, ~, ~, ~, rocComboStruct, mi] = ...
  classify_kfold(rfCentersCombined, channelLocationsAll, 10, plotOn, classifier_type)
auc = rocComboStruct.netAuc;




%% plot all rf centers, those for good and bad sessions, and their mean
rfCentersCombined = cell2mat(rfCentersAll');
%Plot all rf centers
scatter(rfCentersCombined(:,1), rfCentersCombined(:,2), 30, 'filled', 'CData', [.6 .6 .6]); hold on;
rfCentersUnder40 = cell2mat(rfCentersAll(below40SessInds)');
rfCentersAbove90 = cell2mat(rfCentersAll(above90SessInds)');
numUnder40 = size(rfCentersUnder40,1);
numAbove90 = size(rfCentersAbove90, 1);

scatter(rfCentersUnder40(:,1), rfCentersUnder40(:,2), 30, 'filled', 'CData', [.4 .4 1]);
scatter(rfCentersAbove90(:,1), rfCentersAbove90(:,2), 30, 'filled', 'CData', [.4 1 .4]);

[mnUnder40, semUnder40] = mean_sem(rfCentersUnder40);
[mnAbove90, semAbove90] = mean_sem(rfCentersAbove90);
scatter(mnUnder40(1), mnUnder40(2), 50, 'b', 'filled');hold on
%Errorbars under 40
plot([mnUnder40(1)-semUnder40(1) mnUnder40(1)+semUnder40(1)],[mnUnder40(2) mnUnder40(2)],'b', 'LineWidth', 2)
plot([mnUnder40(1) mnUnder40(1)],[mnUnder40(2)-semUnder40(2) mnUnder40(2)+semUnder40(2)],'b', 'LineWidth', 2)

scatter(mnAbove90(1), mnAbove90(2), 50, 'filled', 'g');
plot([mnAbove90(1)-semAbove90(1) mnAbove90(1)+semAbove90(1)],[mnAbove90(2) mnAbove90(2)],'g', 'LineWidth', 2)
plot([mnAbove90(1) mnAbove90(1)],[mnAbove90(2)-semAbove90(2) mnAbove90(2)+semAbove90(2)],'g', 'LineWidth', 2)

plot([-850 0 850 0 -850], [0 850 0 -850 0],'k','LineWidth', 1); %the square
grid on; axis equal;
axis([-860 860 -860 860]);
set(gca,'XTick',[-850:170:850], 'YTick', [-850: 170: 850]);
%also plot all popvecs that actually occurred in all sessions, in very
%light color
shg
[~, pMeanRfCenter] = manova1([rfCentersUnder40; rfCentersAbove90], [zeros(numUnder40,1); ones(numAbove90,1)]); %0.1154
%print -depsc2 rfCentersAllPlot


% covariance matrix
% covUnder40 = cov(rfCentersUnder40);
% covAbove90 = cov(rfCentersAbove90);



%% plot histogram of centers as distance from edge or something....likely get 
%higher as you go toward edge, I would guess


%% Plot histogram of centers along anteroposterior
%Either subplot, or stack bars
binCenters  = [-850:50:850];  %85
rfYUnder40 = rfCentersUnder40(:,2);
rfYAbove90 = rfCentersAbove90(:,2);

nUnder40 = hist(rfYUnder40, binCenters);
freqUnder40 = nUnder40/sum(nUnder40);
nAbove90 = hist(rfYAbove90, binCenters);
freqAbove90 = nAbove90/sum(nAbove90);
bothHists = [freqUnder40; freqAbove90]

% draw them
bar(binCenters, bothHists', 1, 'grouped');hold on
ax = get(gca);
cat = ax.Children;
%set the first bar chart style
set(cat(2),'FaceColor',[0 0 1],'BarWidth',1);
%set the second bar chart style
set(cat(1),'FaceColor',[0 1 0],'BarWidth',1);
axis([-925 950 0 0.45])
legend({'<40','>90'}, 'Location', 'NorthWest')
xlabel('y position')
ylabel('Frequency of peaks')
mnYUnder = mean(rfYUnder40);
mnYAbove = mean(rfYAbove90);
scatter(mnYUnder,0, 'filled', 'b');
scatter(mnYAbove,0, 'filled', 'g');
shg;
%print -depsc2 yCenterHistograms


%% Peak z score for those 40% and under, and those 90% and over
peakZScoresOverall = cell2mat(peakZScoresAll');

%plot all z scores
hist(peakZScoresOverall, [min(peakZScoresOverall)-0.1: 0.5: max(peakZScoresOverall)+0.1])
numAbove = length(find(peakZScoresOverall > 3))
percentAbove = numAbove/length(peakZScoresOverall)  %67 percent

peakZScoresUnder40 = cell2mat(peakZScoresAll(below40SessInds)');
peakZScoresAbove90 = cell2mat(peakZScoresAll(above90SessInds)');
[mnZBelow40, semZBelow40] = mean_sem(peakZScoresUnder40);
[mnZAbove90, semZAbove90] = mean_sem(peakZScoresAbove90);
[~, pZScores] = ttest2(peakZScoresUnder40, peakZScoresAbove90); %p = 0.91

%% Bar plot
bar([0 1], [mnZBelow40 mnZAbove90],'k');hold on;
plot_errorbars([0 1], [mnZBelow40 mnZAbove90], [semZBelow40 semZAbove90],0,'k',1);
axis([-0.5 1.5 0 30])
set(gca,'XTickLabel',{'<40', '>90'}, 'YTick', [0:10:40])

xlabel('Session Performance')
ylabel('Peak Z Score');
shg

print -depsc2 zBarCompare

%% Receptive field diameter
diametersOverall = cell2mat(rfDiametersAll');
[mnDiameter, semDiameter] = mean_sem(diametersOverall);  %533 +/- 8 (283 channels)
diamUnder40 = cell2mat(rfDiametersAll(below40SessInds)');
diamAbove90 = cell2mat(rfDiametersAll(above90SessInds)');
[mnDiamBelow40, semDiamBelow40] = mean_sem(diamUnder40);
[mnDiamAbove90, semDiamAbove90] = mean_sem(diamAbove90);
[~, pDiam] = ttest2(diamUnder40, diamAbove90); %0.009

%% Bar graph
bar([0 1], [mnDiamBelow40 mnDiamAbove90],'k');hold on;
plot_errorbars([0 1], [mnDiamBelow40 mnDiamAbove90], [semDiamBelow40 semDiamAbove90],0,'k',1);
axis([-0.5 1.5 0 300])
set(gca,'XTickLabel',{'<40', '>90'})

xlabel('Session Performance')
ylabel('RF Diameter');
shg

print -depsc2 rfDiamCompare

%% Draw all popvecs, and circle withmn diamter in middle
scatter(popVecspace(:,1),popVecspace(:,2), 50, 'filled', 'CData',[.7 .7 .7]);   hold on
plot([-850 0 850 0 -850], [0 850 0 -850 0],'k','LineWidth', 1); %the square
plot_circle([0 0], mnDiameter/2, 'k', 2)
grid on; axis equal;
axis([-850 850 -850 850]);
set(gca,'XTick',[-850:170:850], 'YTick', [-850: 170: 850]);shg



    
%% plot all possible with actual from one session
%Which session?
scatter(popVecspace(:,1), popVecspace(:,2),20, 'CData', [.7 .7 .7]);hold on;
%%Plot all possible, with actual on top
detectorCoords = [1 -1 -1 1;-1 -1 1 1];  %2x4, channel i has column i 
stimFreqs = [0 10 100 150 250 350 425]; 
allStimCombos=combvec(stimFreqs, stimFreqs, stimFreqs, stimFreqs)';
%Population vector P:
%P = sum(magnitude_i*direction_i)
popVecspace = (detectorCoords * allStimCombos')';
%Plot all population vectors
scatter(popVecspace(:,1),popVecspace(:,2), 20, 'filled', 'CData',[.7 .7 .7]);   hold on
    for trialNum = 1: numTrials
        if rem(trialNum,40) == 0
            disp(['Plotting trial ' num2str(trialNum) '/' num2str(numTrials)]);
        end
        popVec=popvecSequences{trialNum};
        numSubstim = size(popVec,1);
        colorVals = linspace(.98,.2,numSubstim);
        for subStimNum = 1: numSubstim
            %colSub=colorVals(subStimNum)*ones(1,3);
            if corrIncorrVals(trialNum)
                colSub=[colorVals(subStimNum) 1 colorVals(subStimNum)]; % ones(1,3); %green hues
            else
                colSub=[1 colorVals(subStimNum) colorVals(subStimNum)]; 
            end
            scatter(popVec(subStimNum,1), popVec(subStimNum,2), 30, 'filled', 'CData', colSub);hold on
            if subStimNum < numSubstim
               plot([popVec(subStimNum,1), popVec(subStimNum+1,1)], ...
                   [popVec(subStimNum,2), popVec(subStimNum+1,2)],'Color',colSub, 'LineWidth', 1);
            end
        end
    end
    plot([-850 0 850 0 -850], [0 850 0 -850 0],'k','LineWidth', 2); %the square
    title([animalName ' ' strrep(dateStr, '_', '/') ' all ' num2str(numTrials) ' trials'])
    grid on; axis equal;
    axis([-850 850 -850 850]);
    set(gca,'XTick',[-850:170:850], 'YTick', [-850: 170: 850])




%****************************
%% Plot psth(s) for jasmin 11/21/14
%from laoding saving prep load neuronalAnalysis from that date folder and plot:
%(Run code at 'Calculate all full-trial psths')
cd 'C:\Users\Eric\Dropbox\IR6\IR-4\StimRecord\analysis\dataCombined'
%print -depsc2 jasminPSTHs

%% Plot contour plots for same
cd 'C:\Users\Eric\Dropbox\IR6\IR-4\StimRecord\analysis\dataCombined'
%print -depsc2 jasminContour




%% Response by number stim


%% Response by number stim
%allRespByNumC  (numChannels x 5)
[mnZByNum, semZByNum] = mean_sem(respByNumCombined);
[pAll, tableAll, statsAll] = anova1(respByNumCombined);
multcompare(statsAll);
%
%%
bar([0:4], mnZByNum, 'k');hold on
plot_errorbars([0:4], mnZByNum, semZByNum, 0, 'k', 1);
axis([-0.5 4.5 -.1 5.7])
set(gca,'YTick',[0:5])
xlabel('Number Channels Active')
ylabel('Response Magnitude (z)')
title('Treatment animals')
shg

%%



%respByNumCombined (numChannels x 5)
%alize
negInds =[];
posInds = [];
posRespByNum =[];
posRespByNum =[];
negRespByNum = [];
negRespByNum = [];
for chanNum = 1: numChannelsTotal
    tmpDat = respByNumCombined(chanNum,:);
    minVal = min(tmpDat);
    if (minVal <0) && (abs(minVal) > max(tmpDat))
        negRespByNum = [negRespByNum; tmpDat];
        negRespByNum = [negRespByNum; tmpDat/minVal];
        negInds = [negInds ;chanNum];
    else
        posRespByNum = [posRespByNum; tmpDat];
        posRespByNum = [posRespByNum; tmpDat/max(tmpDat)];
        posInds = [posInds ;chanNum];
    end
end
numPos= size(posRespByNum,1);
numNeg = size(negRespByNum,1);
percNeg = 100*(numNeg/numChannelsTotal);
percPos = 100*(numPos/numChannelsTotal);
disp(['Positive: ' num2str(percPos) '% (' num2str(numPos) '/' num2str(numChannelsTotal) ')'])
disp(['Negative: ' num2str(percNeg) '% (' num2str(numNeg) '/' num2str(numChannelsTotal) ')'])
allRespByNumCorrected = [posRespByNum; -negRespByNum];
[mnPos, semPos] = mean_sem(posRespByNum);
[mnNeg, semNeg] = mean_sem(negRespByNum);
[mnAll, semAll] = mean_sem(allRespByNumCorrected);
[mnPos, semPos] = mean_sem(posRespByNum);
[mnNeg, semNeg] = mean_sem(negRespByNum);
[mnAll, semAll] = mean_sem([posRespByNum;negRespByNum]);


%% significant differences?
[pPos, tablePos, statsPos] = anova1(posRespByNum);
multcompare(statsPos);
[pNeg, tableNeg, statsNeg] = anova1(-negRespByNum);
multcompare(statsNeg);
[pAll, tableAll, statsAll] = anova1(allRespByNumCorrected);
multcompare(statsAll);

%% bar plot of resp vs number (raw z scores)
bar([0:4], mnAll, 'k');hold on
plot_errorbars([0:4], mnAll, semAll, 0, 'k',1);hold on
xlabel('Num active')
ylabel('z score');
title(['All responses ' num2str(numChannelsTotal) ' channels'])
set(gca,'YLim', [0 10]);

%% subdivide resp by overall/+/- responders
subplot(3,1,1)
bar([0:4], mnAll, 'k');hold on
plot_errorbars([0:4], mnAll, semAll, 0, 'k',1);hold on
xlabel('Num active')
ylabel('z score');
title(['All responses ' num2str(numChannelsTotal) ' channels'])
set(gca,'YLim', [0 12]);

subplot(3,1,2)
bar([0:4], mnPos, 'k');hold on
plot_errorbars([0:4], mnPos, semPos, 0, 'k',1)
xlabel('Num active')
ylabel('z score');
title(['Positive response ' num2str(percPos) '% (' num2str(numPos) '/' num2str(numChannelsTotal) ')'])
set(gca,'YLim', [0 12]);

subplot(3,1,3)
bar([0:4], -mnNeg, 'k');hold on; %- to make positive
plot_errorbars([0:4], -mnNeg, semNeg, 0, 'k',1); %- to make positive
xlabel('Num active')
ylabel('z score');
title(['Negative response ' num2str(percNeg) '% (' num2str(numNeg) '/' num2str(numChannelsTotal) ')'])
set(gca,'YLim', [0 0.6])

%print -depsc2 respByNumber

%% examine baselines positive vs negative responders
[baselineMn, baselineSem] = mean_sem(baselineMeansAll);
[baselineNegMn, baselineNegSem] = mean_sem(baselineMeansAll(negInds)); %.56/.054
[baselinePosMn, baselinePosSem] = mean_sem(baselineMeansAll(posInds));  %.79/.06
[~, pBaseMnTest] = ttest2(baselineMeansAll(negInds), baselineMeansAll(posInds))



