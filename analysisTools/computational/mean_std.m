function [mn, std_dev]=mean_std(data)
%[mn, std_dev]=mean_std(data)
%
%returns mean and sem of data, which is assumed to be Nxp (N data points in p dimensions)
%Mean and sem are 1xp

if size(data,1)==1
    warning('In function mean_std: data has only one row--check to see if you should transpose data.')
end


num_dims=size(data,2); %number of data points
mn=nanmean(data);
std_dev=zeros(1,num_dims);

for i=1:num_dims
    std_dev(i)=nanstd(data(:,i));
end