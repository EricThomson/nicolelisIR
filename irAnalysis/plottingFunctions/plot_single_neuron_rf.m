function plot_single_neuron_rf (channelNum, channelName, maxFreq,contourHeights,xMaxVal,xVals,yMaxVal,yVals,rfValsFiltered)
% plots receptive field for a single channel
% channelNum = number of current channel
% channelName = name of current channel
% maxFreq = maximum frequency of stimulation
% contourHeights = height representing strength of receptive field preference
% xMaxVal = maximum value, x coordinate
% yMaxVal = maximum value, y coordinate
% xVals = coordinates of receptive field lines, x coordinate
% yVals = coordinates of receptive field lines, y coordinate
% rfValsFiltered = values of receptive field 
% Nicolelis lab 4/7/17

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

