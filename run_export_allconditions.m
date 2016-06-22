
load('pcnets_options.mat')
% EXPORT TABLE ONLY
tableMetrics = table();
for ii=1:length(opts.conditions)
	filename = opts.outputFiles{ii};
	load(filename,'results');
	netopts = results{1}.netopts;
	for cc=1:length(results)
		%%%%%%% Export to Tables %%%%%%%%%
		condition = struct('SubjectID',{['Subject_' opts.subjIDs{ii}{cc}]}, 'idxThreshold', [1], ....
										'StimLabel',opts.conditions(ii),'Resample','Resample0','gitcommit',thiscommit);
		tableMetrics = vertcat(tableMetrics,exportNetworkMetrics2Table(results{cc}.metrics,condition));
	end
end
save(['Data/NetworkMetrics_' netopts.date '.mat'],'tableMetrics', 'thiscommit');
writetable(tableMetrics,['Data/NetworkMetrics' netopts.date '.csv']);
