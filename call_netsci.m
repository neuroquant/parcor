function [metrics, opts] = call_netsci(A,opts)


	bct_funs = opts.bct_funs;
	n_thresh = opts.n_thresh ;
	isWeighted = opts.isWeighted;
	isScaled = opts.isScaled;
	isDirected = opts.isDirected;
	Ci = opts.Ci;
	alt_Ci = opts.alt_Ci;
	alt_Ci_names = opts.alt_Ci_names;
	

	% Test binary or weighted

	% Check metric name

	% If weighted: Check distanceo or similarity

	% Check directed or undirected

	p=size(A,1);
	tau_start = opts.tau_start;
	tau_stop =  opts.tau_stop;

	if(opts.logscale)
		taus = logspace(log10(tau_start),log10(tau_stop),n_thresh);
	else
		taus = linspace(tau_start,tau_stop,n_thresh);
  end
	opts.taus = taus;

	% Initialize output
	metrics = {};
	metrics.global = [];
	metrics.nodal  = [];
	metrics.centralization = [];
	metrics.name = '';
	metrics.number = 1;

for metric_no = 1:length(opts.bct_num) %*sum([opts.number{opts.bct_num}]))
 
	if(metric_no==1)
		curr_metric_size = length(metrics);
	else
		curr_metric_size = length(metrics)+1;
	end
	bct_num = opts.bct_num(metric_no);

	if(isWeighted)
			netstat_global = [];
			netstat_nodes = zeros(p,n_thresh*opts.number{bct_num});
		for tau=1:length(taus)
			Sighat = A; Sighat(find(eye(p))) = 0; Sighat = abs(Sighat);
			softSig = triu((abs(Sighat)-taus(tau)),1);
			threshSig = softSig.*(softSig>=0) + 0.*(softSig<0);
			softthreshSig = sign(Sighat).*threshSig;
			softthreshSig = softthreshSig + softthreshSig';
			assert(sum(diag(softthreshSig)<=0)~=0,'Negative values on diagonal');
			affinitySig = eye(p)+exp(-(softthreshSig).^2/.01);
			switch func2str(bct_funs{bct_num})
				case 'betweenness_wei'
					tmp_stats = feval(bct_funs{bct_num},affinitySig);
					if(isScaled)
						tmp_stats = tmp_stats/((p-1)*(p-2));
					end
				case 'current_flow_centrality'
					tmp_stats = feval(@callnetworkx,affinitySig,1,0);
					if(isScaled==1)
						tmp_stats(:,1) = tmp_stats(:,1)/((p-1)*(p-2));
					elseif(isScaled==2)
						tmp_stats(:,1) = (tmp_stats(:,1)-min(tmp_stats(:,1)))./(max(tmp_stats(:,1))-min(tmp_stats(:,1))); 
						tmp_stats(:,2) = (tmp_stats(:,2)-min(tmp_stats(:,2)))./(max(tmp_stats(:,2))-min(tmp_stats(:,2))); 						
					end
				case 'efficiency_wei'
					tmp_stats = feval(bct_funs{bct_num},affinitySig,1);
					
				case 'community_participation'
					assert(~isempty(alt_Ci),'Community affiliation is empty in opts.Ci'); 
					assert(length(alt_Ci)==p,'Community affiliation not specified for all nodes');
					tmp_stats = zeros(size(alt_Ci));
					for cc=1:size(alt_Ci,2);
						tmp_stats(:,cc) = feval(@community_participation,abs(softthreshSig),1*alt_Ci(:,cc),0);
						% if(isScaled)
						% 	tmp_stats(:,cc) = (tmp_stats(:,cc)-min(tmp_stats(:,cc)))./(max(tmp_stats(:,cc))-min(tmp_stats(:,cc)));
						% end
					end
				case 'participation_coef'
					assert(~isempty(Ci),'Community affiliation is empty in opts.Ci'); 
					assert(length(Ci)==p,'Community affiliation not specified for all nodes');
					tmp_stats = feval(bct_funs{bct_num},abs(softthreshSig),Ci,0);
					if(isScaled)
						tmp_stats = (tmp_stats-min(tmp_stats))./(max(tmp_stats)-min(tmp_stats)); 
					end
				case 'eigenvector_centrality_und'
					tmp_stats = feval(bct_funs{bct_num},softthreshSig);
					if(isScaled)
						tmp_stats = (tmp_stats-min(tmp_stats))./(max(tmp_stats)-min(tmp_stats)); 
					end
				otherwise
					tmp_stats = feval(bct_funs{bct_num},affinitySig);
			end
			if(~isreal(tmp_stats))
				warning('Network metrics are complex. Taking absolute values');
				tmp_stats = abs(tmp_stats);
			end
			netstat_nodes(:,tau) = tmp_stats(:,1);
			if(strcmp('current_flow_centrality',func2str(bct_funs{bct_num})))
				netstat_nodes(:,length(taus)+tau) = tmp_stats(:,2);
			end
			if(strcmp('community_participation',func2str(bct_funs{bct_num})))
				for cc=2:size(tmp_stats,2);
					netstat_nodes(:,(cc-1)*length(taus)+tau) = tmp_stats(:,cc);
				end
			end
		end
	else

		netstat_global = zeros(1, n_thresh);
		netstat_nodes = zeros(p,n_thresh*opts.number{bct_num});

		for tau=1:length(taus)
			Sighat = A; Sighat(find(eye(p))) = 0; Sighat = abs(Sighat);
			assert(sum(diag(Sighat)<=0)~=0,'Negative values on diagonal');
			switch func2str(bct_funs{bct_num})
			case 'betweenness_bin'
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)));
				if(isScaled)
					tmp_stats = tmp_stats/((p-1)*(p-2));
				end
			case 'current_flow_centrality'
				tmp_stats = feval(@callnetworkx,1*(abs(Sighat)>taus(tau)),0,0);
				if(isScaled==1)
					tmp_stats(:,1) = tmp_stats(:,1)/((p-1)*(p-2));
				elseif(isScaled==2)
					tmp_stats(:,1) = (tmp_stats(:,1)-min(tmp_stats(:,1)))./(max(tmp_stats(:,1))-min(tmp_stats(:,1))); 
					tmp_stats(:,2) = (tmp_stats(:,2)-min(tmp_stats(:,2)))./(max(tmp_stats(:,2))-min(tmp_stats(:,2))); 						
				end
			case 'efficiency_bin'
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)),1);
				
			case 'community_participation'
				assert(~isempty(alt_Ci),'Community affiliation is empty in opts.Ci'); 
				assert(length(alt_Ci)==p,'Community affiliation not specified for all nodes');
				tmp_stats = zeros(size(alt_Ci));
				for cc=1:size(alt_Ci,3);
					tmp_stats(:,cc) = feval(@community_participation,1*(abs(Sighat)>taus(tau)),alt_Ci(:,cc),0);
					% if(isScaled)
					% 	tmp_stats(:,cc) = (tmp_stats(:,cc)-min(tmp_stats(:,cc)))./(max(tmp_stats(:,cc))-min(tmp_stats(:,cc)));
					% end
				end
					
			case 'participation_coef'
				assert(~isempty(Ci),'Community affiliation is empty in opts.Ci'); 
				assert(length(Ci)==p,'Community affiliation not specified for all nodes');
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)),Ci,0);
				if(isScaled)
					tmp_stats = (tmp_stats-min(tmp_stats))./(max(tmp_stats)-min(tmp_stats)); 
				end
			case 'eigenvector_centrality_und'
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)));	
			otherwise
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)));
			end
			if(~isreal(tmp_stats))
				warning('Network metrics are complex. Taking absolute values');
				tmp_stats = abs(tmp_stats);
			end
			netstat_nodes(:,tau) = tmp_stats(:,1);
			if(strcmp('current_flow_centrality',func2str(bct_funs{bct_num})))
				netstat_nodes(:,length(taus)+tau) = tmp_stats(:,2);
			end
			if(strcmp('community_participation',func2str(bct_funs{bct_num})))
				for cc=2:size(tmp_stats,2);
					netstat_nodes(:,(cc-1)*length(taus)+tau) = tmp_stats(:,cc);
				end
			end
		end

	end

		metrics(curr_metric_size).global = netstat_global;
		metrics(curr_metric_size).nodal = netstat_nodes(:,1:length(taus));
		if(strcmp('current_flow_centrality',func2str(bct_funs{bct_num})))
			metrics(curr_metric_size+1).nodal = netstat_nodes(:,length(taus)+1:2*length(taus));
		end
		if(strcmp('community_participation',func2str(bct_funs{bct_num})))
			for cc=2:size(alt_Ci,2);
				metrics(curr_metric_size+cc-1).nodal = netstat_nodes(:,(cc-1)*length(taus)+1:cc*length(taus));
				disp('More community metrics added')
			end
		end

		switch func2str(bct_funs{bct_num})
		case 'betweenness_bin'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case 'efficiency_bin'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case 'clustering_coef_bu'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case	'eigenvector_centrality_und'
			% To Do. add Normalization option to opts
			tmpcentralization = zeros(size(netstat_global));
			for tau=1:length(taus)
			 tmpcentralization(:,tau)= centrality2centralization(netstat_nodes(:,tau),'eigenvector',opts.normalizeCentralization);
			end
			metrics(curr_metric_size).centralization = tmpcentralization;
			metrics(curr_metric_size).name = 'EigenvectorCentrality';
		case	'rich_club_bu'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case	'betweenness_wei'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case	'clustering_coef_wu'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case	'efficiency_wei'
			metrics(curr_metric_size).name = 'WeightedEfficiency';
			for tau=1:length(taus)
			 tmpcentralization(:,tau)= centrality2centralization(netstat_nodes(:,tau),'eigenvector',0);
			end
			metrics(curr_metric_size).centralization = tmpcentralization;
			clear tmpcentralization;
		case	'current_flow_centrality'
			metrics(curr_metric_size).name = 'RandomWalkBetweenness';
			tmpcentralization = zeros(size(netstat_global));
			for tau=1:length(taus)
			 		tmpcentralization(:,tau)= centrality2centralization(netstat_nodes(:,tau),'betweenness',opts.normalizeCentralization);
			end
			metrics(curr_metric_size).centralization = tmpcentralization;
			clear tmpcentralization

			metrics(curr_metric_size+1).name = 'RandomWalkCloseness';
			tmpcentralization = zeros(size(netstat_global));
			for tau=1:length(taus)
					tmpcentralization(:,tau)= centrality2centralization(netstat_nodes(:,length(taus)+tau),'closeness',opts.normalizeCentralization);
			end
			metrics(curr_metric_size+1).centralization = tmpcentralization;

		case	'rand_hits'
			metrics(curr_metric_size).name = 'RegularizedHITS';
		case	'rich_club_wu'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case 'participation_coef'
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
		case 'community_participation'
			metrics(curr_metric_size).name = [func2str(bct_funs{bct_num}) '_' alt_Ci_names{1}];		
			for cc=2:length(alt_Ci_names)
				metrics(curr_metric_size+cc-1).name = [func2str(bct_funs{bct_num}) '_' alt_Ci_names{cc}];
			end
		otherwise
			metrics(curr_metric_size).name = func2str(bct_funs{bct_num});
	end

end
