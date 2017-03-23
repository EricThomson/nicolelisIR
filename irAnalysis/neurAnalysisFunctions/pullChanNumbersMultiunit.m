function channelNumbers = pullChanNumbersMultiunit(channelNames)
%channelNumbers = pullChanNumbersMultiunit(channelNames)
%Pulls which channel numbers are multiunit
%
%Not really general yet

channelNumbersString = cellfun(@(x) strrep((strrep(x,'sig','')), 'i',''), ...
    channelNames,'UniformOutput', false);
channelNumbers = cellfun(@(x) str2num(x), channelNumbersString);