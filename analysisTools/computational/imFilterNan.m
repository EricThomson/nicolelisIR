function dataFiltered = imFilterNan(data, imageFilter, plot_on)
%Filtering an image that includes NaNs, returns dataFIltered same size as data. 
%   dataFiltered = imFilterNan(data, imageFilter, plot_on)
%Wrapper for nanconv function from Matlab File Exchange:
%   http://stackoverflow.com/a/29897777/1886357


%Pretty much didn't need a new function
dataFiltered = nanconv(data, imageFilter, 'nanout');


if plot_on
    figure;
    %original data
    subplot(1,2,1);
    contourf(data); set(gca,'YDir', 'normal');hold on
    axis equal; axis tight;
    title('Original Data')
    %filtered data
    subplot(1,2,2)
    contourf(dataFiltered); set(gca,'YDir', 'normal');hold on
    axis equal;
    axis tight
    title('Filtered Data')
end






%:)
