% % basepath = '../datasets/PNASresting5.16/'  % files are RESTING*5.16
% % basepath = 'PNASData/pMFG_bilat_4.14'
% basepath = 'PNASData/pMFG_left_5.11/'
% listfiles = dir(sprintf('%s/*pMFG*',basepath));

load('pcnets_options.mat')

for ii=1:length(opts.conditions)
	basepath = opts.basepath{ii};
	listfiles = dir([opts.inputFiles{ii}]);
	listfiles = listfiles(find([listfiles.bytes]~=0));
	[m p] = size(dlmread(sprintf('%s/%s',basepath,listfiles(1).name),'\t',1,0));
	n = length(listfiles);
	Xdata = zeros(m,p,n);
	for cc=[1:n]
		cc
		Xdata(:,:,cc) = dlmread(sprintf('%s/%s',basepath,listfiles(cc).name),'\t',1,0);
	end
	nonzeroroiIdx = find(squeeze(sum(sum(Xdata,1),3))~=0);
	% Xdata = Xdata(:,nonzeroIdx,:);
	nonzerosubIdx = find(squeeze(sum(sum(Xdata,1),2))~=0);
	% Xdata = Xdata(:,:,nonzeroIdx);

	roinames = readtable('roinames.csv', 'ReadVariableNames',1);
	roinames = roinames.Properties.VariableNames;
	rmcerebellar_idx = find(cellfun(@isempty,strfind(roinames,'cerebellum')));
	Xdata = Xdata(:,rmcerebellar_idx,:);

	[m p n] = size(Xdata);

	Xdata = bsxfun(@minus,Xdata,mean(Xdata,1));

	tmpname = regexp(listfiles(1).name,{'[+0-9]_'},'split')
	filename = tmpname{1}(2);
	Data = Xdata;
	roinames = cell2struct(roinames,'name',30)
	save([opts.outputFiles{ii} '.mat'],'Data','roinames','listfiles','nonzeroroiIdx','nonzerosubIdx');
end

isVisible = 0; % Set to opts.isVisible;

if(isVisible)
	figure(1)
	subjno =2;
	n_regions = 10;
	for ii=1:n_regions
		subplot(n_regions,1,ii);
		plot(Xdata(:,ii,subjno));
	end

	figure(2)
	subplot(2,2,1)
	imagesc(corr(Xdata(:,:,subjno)')); title('Temporal Correlation'); colorbar;
	subplot(2,2,2)
	imagesc(corr(Xdata(:,:,subjno+1)')); title('Temporal Correlation'); colorbar;
	subplot(2,2,3)
	imagesc(corr(Xdata(:,:,subjno+2)')); title('Temporal Correlation'); colorbar;
	subplot(2,2,4)
	imagesc(corr(Xdata(:,:,subjno+4)')); title('Temporal Correlation'); colorbar;

	figure(3)
	subplot(2,2,1)
	imagesc(corr(Xdata(:,:,subjno))); title('Spatial Correlation'); colorbar;
	subplot(2,2,2)
	imagesc(corr(Xdata(:,:,subjno+1))); title('Spatial Correlation'); colorbar;
	subplot(2,2,3)
	imagesc(corr(Xdata(:,:,subjno+2))); title('Spatial Correlation'); colorbar;
	subplot(2,2,4)
	imagesc(corr(Xdata(:,:,subjno+4))); title('Spatial Correlation'); colorbar;
end
