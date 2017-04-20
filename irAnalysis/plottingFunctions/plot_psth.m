function plot_psth (channelNumber, chanNamesMultiunit,  fullMeanPsth )
% plot psth
% channelNumber = channel number
% chanNamesMultiUnit = array of channel names
% fullMeanPsth = psth data
% Nicolelis Lab 4/7/17
    channelName = chanNamesMultiunit{channelNumber};
    numChannels = length(chanNamesMultiunit);
    if numChannels <= 16
        subplot(4,4,channelNumber);
    else
        subplot(5,4,channelNumber);
    end
    tickBins = [26:25:242];  %-500, -250, 0
    bar(fullMeanPsth, 1, 'k');axis tight
    set(gca,'XTick', tickBins); %[76 89 102 135 148 161]);
    set(gca,'XTickLabel', [-0.5:0.25:1.5])%{'S1', 'S2', 'S3', 'SN-2', 'SN-1', 'S_N'});
    shg
    axis tight
    title(channelName)
    shg;
end



