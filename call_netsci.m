function [metrics, opts] = call_netsci(A,opts)


	bct_funs = opts.bct_funs;
	n_thresh = opts.n_thresh ;
	isWeighted = opts.isWeighted;
	isScaled = opts.isScaled;
	isDirected = opts.isDirected;

	% Test binary or weighted

	% Check metric name

	% If weighted: Check distanceo or similarity

	% Check directed or undirected

	p=size(A,1);
	tau_start = .01;
	tau_stop = .3*max(max(triu(abs(A),1)));

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

for metric_no = 1:length(opts.bct_num)

	bct_num = opts.bct_num(metric_no);

	if(isWeighted)
			netstat_global = 0;
			netstat_nodes = zeros(p,1);
		for tau=1:length(taus)
			Sighat = A; Sighat(find(eye(p))) = 0; Sighat = abs(Sighat);
			softSig = triu((abs(Sighat)-taus(tau)),1);
			threshSig = softSig.*(softSig>=0) + 0.*(softSig<0);
			softthreshSig = sign(Sighat).*threshSig;
			softthreshSig = softthreshSig + softthreshSig';
			assert(sum(diag(softthreshSig)<=0)~=0,'Negative values on diagonal');
			affinitySig = exp(-(softthreshSig).^2/.15);
			switch func2str(bct_funs{bct_num})
			case 'betweenness_wei'
				tmp_stats = feval(bct_funs{bct_num},affinitySig);
				if(isScaled)
					tmp_stats = tmp_stats/((p-1)*(p-2));
				end
			case 'efficiency_wei'
				tmp_stats = feval(bct_funs{bct_num},affinitySig,1);
			otherwise
				tmp_stats = feval(bct_funs{bct_num},softthreshSig);
			end
			if(~isreal(tmp_stats))
				warning('Network metrics are complex. Taking absolute values');
				tmp_stats = abs(tmp_stats);
			end
			netstat_nodes(:,tau) = tmp_stats;
		end
	else

		netstat_global = zeros(1, n_thresh);
		netstat_nodes = zeros(p,n_thresh);

		for tau=1:length(taus)
			Sighat = A; Sighat(find(eye(p))) = 0; Sighat = abs(Sighat);
			assert(sum(diag(Sighat)<=0)~=0,'Negative values on diagonal');
			switch func2str(bct_funs{bct_num})
			case 'betweenness_bin'
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)));
				if(isScaled)
					tmp_stats = tmp_stats/((p-1)*(p-2));
				end
			case 'efficiency_bin'
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)),1);
			otherwise
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)));
			end
			if(~isreal(tmp_stats))
				warning('Network metrics are complex. Taking absolute values');
				tmp_stats = abs(tmp_stats);
			end
			netstat_nodes(:,tau) = tmp_stats;
		end

	end

		metrics(metric_no).global = netstat_global;
		metrics(metric_no).nodal = netstat_nodes;

		switch func2str(bct_funs{bct_num})
		case 'betweenness_bin'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		case 'efficiency_bin'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		case 'clustering_coef_bu'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		case	'eigenvector_centrality_und'
			% To Do. add Normalization option to opts
			tmpcentralization = zeros(size(netstat_global));
			for tau=1:length(taus)
			 tmpcentralization(:,tau)= centrality2centralization(netstat_nodes(:,tau),'eigenvector',opts.normalizeCentralization);
			end
			metrics(metric_no).centralization = tmpcentralization;
			metrics(metric_no).name = 'EigenvectorCentrality';
		case	'rich_club_bu'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		case	'betweenness_wei'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		case	'clustering_coef_wu'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		case	'efficiency_wei'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		case	'current_flow_betweenness'
			metrics(metric_no).name = 'RandomWalkBetweenness';
			tmpcentralization = zeros(size(netstat_global));
			for tau=1:length(taus)
			 tmpcentralization(:,tau)= centrality2centralization(netstat_nodes(:,tau),'betweenness',opts.normalizeCentralization);
			end
			metrics(metric_no).centralization = tmpcentralization;
		case	'rand_hits'
			metrics(metric_no).name = 'RegularizedHITS';
		case	'rich_club_wu'
			metrics(metric_no).name = func2str(bct_funs{bct_num});
		otherwise
			metrics(metric_no).name = func2str(bct_funs{bct_num});
	end

end
