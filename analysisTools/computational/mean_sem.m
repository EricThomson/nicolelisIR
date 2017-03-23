function [mn, sem] = mean_sem(data)
%[mn, sem] = meanSem(data)
%
%Returns mean and sem of data, which is assumed to be Nxp (N data points in p dimensions)
%Mean and sem are 1xp arrays

numDims=size(data,2); %number of dimensions of each vector
numPoints=size(data,1); %number of vectors

if numPoints==1
	
   mn=data;
   sem=zeros(size(data));
   return
else
    mn=nanmean(data);
end

sem= zeros(1, numDims);
divisor = sqrt(numPoints);
for i=1:numDims
    sem(i)=nanstd(data(:,i))/sqrt(numPoints);
end