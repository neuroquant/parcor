% % basepath = '../datasets/PNASresting5.16/'  % files are RESTING*5.16
% % basepath = 'PNASData/pMFG_bilat_4.14'
% basepath = 'PNASData/pMFG_left_5.11/'
% listfiles = dir(sprintf('%s/*pMFG*',basepath));


load('pcnets_options.mat')

% For manual/visual checking of data, uncomment below
saveFiles = 0;
usePause = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(opts.conditions)
%for ii=5
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
		if(saveFiles)
			save([opts.outputFiles{ii} '.mat'],'Data','roinames','listfiles','nonzeroroiIdx','nonzerosubIdx');
	  end

	 isVisible = 1; % Set to opts.isVisible;
	 addpath('../packages/export_fig');
	 if(isVisible)
	 	figure(1)
	 	subjno =2;
	 	n_regions = 10;
	 	for rr=1:n_regions
	 		subplot(n_regions,1,rr);
	 		plot(Xdata(:,rr,subjno));
	 	end

	 if(strfind(opts.conditions{ii},'aMFG'))
		 new_idx = find(cellfun(@isempty,strfind({roinames.name},'AMFG')));
	 elseif(strfind(opts.conditions{ii},'pMFG'))
		 new_idx = find(cellfun(@isempty,strfind({roinames.name},'PMFG')));
	 end
	  [h1 h2] = plot_correlations(Xdata(:,new_idx,:));
	 	figure(h1.figure)
	 	export_fig([opts.outputFiles{ii} '_datacheck_subjects.png'])
	 	figure(h2.figure)
	 	export_fig([opts.outputFiles{ii} '_datacheck_regions.png'])
	 end

	 if(usePause)
	 	pause(5)
	end
   close all;
end
