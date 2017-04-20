function plot_stimulus_population( popvecSequences, animalName, maxFreq, dateStr)
% Use to plot stimulus population vector
% popvecSequences = all vectors
% animalName = animal name
% maxFreq = maximum frequency of stimulation
% printOn = should graph be printed after generation
% dateStr = date
% 4/5/2017 Nicolelis Lab
disp('Plotting all population vector trajectories...this may take a minute.')
figure;

numTrials = length(popvecSequences);

for trialNum = 1: numTrials
    popVec=popvecSequences{trialNum};
    numSubstim = size(popVec,1);
    colorVals = linspace(.98,.2,numSubstim);
    
    if numSubstim %don't try to plot if there are no substim
        for subStimNum = 1: numSubstim
            colSub=colorVals(subStimNum)*ones(1,3);  %for paper
            %                 if corrIncorrVals(trialNum)
            %                     colSub=[colorVals(subStimNum) 1 colorVals(subStimNum)]; % ones(1,3); %green hues
            %                 else
            %                     colSub=[1 colorVals(subStimNum) colorVals(subStimNum)];
            %                 end
            scatter(popVec(subStimNum,1), popVec(subStimNum,2), 25, 'filled', 'CData', colSub);hold on
            
            if subStimNum < numSubstim
                plot([popVec(subStimNum,1), popVec(subStimNum+1,1)], ...
                    [popVec(subStimNum,2), popVec(subStimNum+1,2)],'Color',colSub);
            end
            
        end %for subStimNum
    end %if numSubstim
end
plot([-maxFreq*2 0 maxFreq*2 0 -maxFreq*2], [0 maxFreq*2 0 -maxFreq*2 0],'k','LineWidth', 1); %the square
title([animalName ' ' strrep(dateStr, '_', '/') ' all ' num2str(numTrials) ' trials'])
grid on; axis equal;
axis([-maxFreq*2 maxFreq*2 -maxFreq*2 maxFreq*2]);
set(gca,'XTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
    'XTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
set(gca,'YTick',[-maxFreq*2:maxFreq: maxFreq*2], ...
    'YTickLabel', {'-fmax*2', '-fmax', '0', 'fmax', 'fmax*2'});
gridfix([.85 .85 .85])
set(gcf,'Renderer','OpenGL')
