function plot_3x3_popvec (animalName, maxFreq,  popVec, numSubstimMin, dateStr )
% plots first and last three stimuli locations
% animalName = animal name
% maxFreq = maximum frequency of stimulation
% popVec = array of vectors 
% numSubstimMin = minimal number of substimuli per trial
% dateStr = session date
% Nicolelis lab 4/7/17
figure;
colorVals = linspace(.75,0,numSubstimMin);
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
