function [rfVals, rfFunction, mnResponse, xGrid, yGrid] = rfCalculate(data, rfSpacing, limVal, plot_on)
%[rfVals, rfFunction, mnResponse, xGrid, yGrid] = rfCalculate(data, rfSpacing, limVal, plot_on)
%
%Calculates receptive field of data, assuming N x 3 data with first
%two columns coordinates.
%binCenters = [-limVal + rfSpacing/2 :rfSpacing: limVal - rfSpacing/2];
%Right now, rfSpacing and limVal are same for x- and y-, so it cannot
%handle weird axis differences....
%%RF parameters (ultimately want rfSpacing and limVal to be parameter)
%
%Returns rf=receptive field, the mean response as a function of x/y
%location
%
binCenters = [-limVal + rfSpacing/2 :rfSpacing: limVal - rfSpacing/2];
[xGrid,yGrid] = meshgrid(binCenters, binCenters);

xydat = data(:,1:2);
uniqueXy = unique(xydat, 'rows');
numUniqueXy = size(uniqueXy, 1);
mnResponse = zeros(numUniqueXy, 1);
for xyValInd = 1: numUniqueXy
    xyValTmp = uniqueXy(xyValInd,:);
    xyInds = find(ismember(xydat, xyValTmp, 'rows'));
    mnResponse(xyValInd) = mean(data(xyInds, 3));
end



%interpolate with tri (uses linear)
rfFunction = TriScatteredInterp(uniqueXy(:,1), uniqueXy(:,2), mnResponse);
rfVals = rfFunction(xGrid,yGrid);

%pRf = uniform2dRfTest(rfVals)

if plot_on
    figure;
    h= contourf(xGrid, yGrid, rfVals);hold on
    %shading interp;
    xlabel('x')
    ylabel('y')
    zlabel('Magnitude');
    axis equal; axis tight;
end
