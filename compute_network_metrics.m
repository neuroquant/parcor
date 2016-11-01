
% Add packages
addpath(genpath('../packages/BCT'))
addpath('../continuity_netsci')
addpath('../continuity_netsci/netfun')


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Pick Metrics %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
bct_funs = {@betweenness_bin, ... %1
			@clustering_coef_bu, ... 	%2
			@efficiency_bin, ...
			@eigenvector_centrality_und, ...
			@rich_club_bu, ...
			@betweenness_wei, ...
			@clustering_coef_wu, ...
			@efficiency_wei, ...
			@current_flow_metrics, ...
			@rand_hits, ...
			@rich_club_wu
			 }; 			%			@rich_club_bd  %			@rich_club_wd

bct_num = [4,8,9];
isWeighted = 1; % Eventually choose this and network metrics with error checking
isScaled = 0;
if(strfind(func2str(bct_funs{bct_num(1)}),'wei')|strfind(func2str(bct_funs{bct_num(1)}),'wu'))
	isWeighted =1;
elseif(isWeighted)
	warning('isWeighted setting likely incorrect')
	if(strfind(func2str(bct_funs{bct_num(1)}),'bin')|strfind(func2str(bct_funs{bct_num(1)}),'bu'))
		isWeighted=0;
	end
end
disp(['Metric Chosen is ' func2str(bct_funs{bct_num(1)}) ', isWeighted: ' num2str(isWeighted)])

%%%%%%%%%%%%%%%%%%%%%%%
n_thresh = 5; % Granularity of shrinkage or thresholding parameter
isLogscale = 1; % Spacing between thresholds logarithmic
%%%%%%%%%%%%%%%%%%%%%%%


% Set options:
% - opts.bct_funs
% - opts.richness_level
% - opts.nthresh (needed for binary networks)
%%%%%%%% Customize %%%%%%%%%%
load('pcnets_options.mat')
netopts.bct_funs = bct_funs;
netopts.n_thresh = n_thresh;
netopts.logscale = isLogscale;
netopts.tau_start = .01;
netopts.tau_stop = .25;
netopts.isWeighted = isWeighted;
netopts.bct_num = bct_num;
netopts.isScaled = isScaled;
netopts.isDirected = 0;
netopts.normalizeCentralization = 0;
netopts.date = datestr(now,'mmm-dd-yyyy-HHMM');
warning('off', 'MATLAB:table:ModifiedVarnames')
%%%%%%%%%%%%%%%%%%%%%%%
tableMetrics = table();
for ii=1:length(opts.conditions)
	filename = opts.outputFiles{ii};
	load(filename,'results');
	for cc=1:length(results)
		cc
		tmpA = results{cc}.Rho;
		% call_netsci
		[tmp_metrics tmp_opts] = call_netsci(tmpA,netopts);
		% save results
		results{cc}.metrics = tmp_metrics;
		results{cc}.netopts = tmp_opts;
		% export results
		%%%%%%% Export to Tables %%%%%%%%%
		condition = struct('SubjectID',{['Subject_' opts.subjIDs{ii}{cc}]}, 'idxThreshold', [1], ....
										'StimLabel',opts.conditions(ii),'Resample','Resample0','gitcommit',thiscommit);
	  tableMetrics = vertcat(tableMetrics,exportNetworkMetrics2Table(results{cc}.metrics,condition));
	end
		% clear temporary variables
		clear tmpA tmp_metrics tmp_opts
	try
		gitlog = evalc('git log -1');
	  thiscommit = {gitlog(16:25)};
	catch
		[unixs unixr] = unix('git log -1')
		thiscommit = unixr;
	end
	save([filename '.mat'],'results', 'thiscommit', '-append');
end
save(['Data/NetworkMetrics_' netopts.date '.mat'],'tableMetrics', 'thiscommit');

currdate = 'Jun-16-2016-1722' ;
load(['Data/NetworkMetrics_' currdate '.mat'],'tableMetrics');
grammplot_metrics

% for cc=3
% 	for ii=1:(length(opts.conditions)-1)
% 			filename = opts.outputFiles{ii};
% 			load(filename,'results');
% 	end
% end


% Remove packages
rmpath(genpath('../packages/BCT'))
rmpath('../continuity_netsci')
rmpath('../continuity_netsci/netfun')
