% load('aMFG_TimeSeries.mat','Data');
% for cc=1:size(Data,3)
% 	cc
% 	results{cc} = parcor(Data(:,:,cc),struct('lambda',[],'visible',0));
% end
% save('aMFG_TimeSeries.mat','results','-append');
% clear;
%
%
% load('rest_TimeSeries.mat','Data');
% for cc=1:size(Data,3)
% 	cc
% 	results{cc} = parcor(Data(:,:,cc),struct('lambda',[],'visible',0));
% end
% save('rest_TimeSeries.mat','results','-append');
% clear;


% Analytic, optimal ridge penalty.
% addpath('../datasets')
%filename = 'SP_pMFG_4.2.mat';
%filename = 'SP_aMFG_4.1.mat';
%filename = 'SP_pMFG_left_5.11.mat';
load('pcnets_options.mat')

for ii=1:length(opts.conditions)
	filename = opts.outputFiles{ii};
	load(filename,'WhData');

	Data = WhData(:,:,find([squeeze(sum(sum(WhData,2),1))]));

	for cc=1:size(Data,3)
		cc
		results{cc} = parcor(Data(:,:,cc),struct('lambda',[],'visible',0,'ridgeType',1));
	end

	% Display global mean network
	global_network = zeros(size(Data,2),size(Data,2));
	useFisherZ = 0;
	for cc=1:length(results)
		cc
		if(useFisherZ)
			global_network = global_network+fisherZ(results{cc}.Rho);
		else
			global_network = global_network+(results{cc}.Rho);
		end
		if(~isreal(global_network))
			disp('Partial Correlation has imaginary values')
		end
	end
	global_network = global_network/length(results);

	% Matrix plot
	% imagesc((global_network)); colormap redBlueCmap(100,100); colorbar; axis image;
	% set(gca,'YTick',[1:30],'YTickLabel',roinames); set(gca,'XTick',[1:30],'XTickLabel',roinames);
	% XYrotalabel

	% % Alternatively
	% load(filename,'roinames');
	% if(isstruct(roinames))
	% 	roinames = {roinames.name};
	% end
	% roinames = regexprep(roinames,'_', ' '); roinames = regexprep(roinames,'50', '');
	% roinames = regexprep(roinames,'campus', ''); roinames = regexprep(roinames,'noSFG gz', '');
	%h = schemanet(fisherZ(global_network).*(abs(global_network)>.1),roinames,[0 1],'summer',2,{[1:12]'; [13:18]'; [19:30]'}); set(gcf,'Color',[.9 .9 .9])
	%print('-dpng','-r600','Data/test_global_network_aMFG')

	save([filename '.mat'],'results','global_network','-append');
	nii_struct = make_nii(global_network);
	save_nii(nii_struct,[filename '.nii']);
	clear Data WhData results;
end
