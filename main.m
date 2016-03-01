load('aMFG_TimeSeries.mat','Data');
for cc=1:size(Data,3)
	cc
	results{cc} = parcor(Data(:,:,cc),struct('lambda',[],'visible',0));
end
save('aMFG_TimeSeries.mat','results','-append');
clear;

load('rest_TimeSeries.mat','Data');
for cc=1:size(Data,3)
	cc
	results{cc} = parcor(Data(:,:,cc),struct('lambda',[],'visible',0));
end
save('rest_TimeSeries.mat','results','-append');
clear;