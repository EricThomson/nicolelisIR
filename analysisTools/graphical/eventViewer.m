function eventViewer(allEvents, eventNames, timeWindow, colors)
%eventViewer(allEvents, eventNames, timeWindow, colors)
%
%Plot event rasters
%
%allEvents: 1xnumEvents cell array of events (each is a numeric array in each cell)
%eventNames: cell array of event names (1xnumEvents cell array of chars)
%timeWindow: 1x2 min and max times
%Plots them using colors (1xnumEvents cell array of colors, e.g., 'r' or [0 1 0])

numEvents = length(allEvents);

for eventNum = 1: numEvents
    tmpEvents = allEvents{eventNum};
    tmpEventsWindowInds = find(tmpEvents >= timeWindow(1) & tmpEvents <= timeWindow(2));
    tmpEvents = tmpEvents(tmpEventsWindowInds);
    numTimes = length(tmpEvents);
    tmpEvents = reshape(tmpEvents, 1, numTimes);
    %yLims = [eventNum-1 eventNum] #used for debugging
    line([tmpEvents; tmpEvents], [zeros(1,numTimes)+eventNum-1; eventNum*ones(1,numTimes)], ...
        'Color', colors{eventNum}, 'LineWidth', 1); 
end

axis([timeWindow 0 numEvents])
xlabel('Time')
ylabel('Event')
set(gca, 'YTick',0.5:(numEvents-1)+0.5);
set(gca,'YTickLabel',eventNames, 'YDir', 'reverse')
shg
  