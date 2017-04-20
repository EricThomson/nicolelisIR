function plot_number_substim(proportionChannelsActivated, responseByNum_indexBased,animalName,dateStr)
% plot number of substim at each interval and z-scores
% proportionChannelsActivated = % of channels activated
% responseByNum_indexBased = response of certain channel
% animalName = name of animal
% dateStr = date of session
% Nicolelis lab 4/7/17

[mnRespByNum_indexBased, semRespByNum] = mean_sem(responseByNum_indexBased);

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

