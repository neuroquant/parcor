function result = parcor(X,opts)
	% parcor.m
	% Takes a multivariate time-series data matrix X and computes sample partial correlation matrix
	% INPUT
	%
	% X:		n_samples x n_regions data matrix. Time-series centered to be zero-mean.
	%
	% opts
	% opts.lambda: 	ridge penalty. 0 by default if n_samples > 3*n_regions else set to .1
	% 
	% OUTPUT
	%
	% result
	% result.Rho: 	n_regions x n_regions partial correlation matrix
	% result.Theta:	n_regions x n_regions inverse covariance  matrix
	% result.Var:	n_regions x 1 ; diagonal matrix of partial variances
	% result.opts:	options parameter
	% result.dfe:	naive degrees of freedom for t-statistic for standard partial correlation. Does not apply if ridge penalty is used. 
	% result.Zstat: If no ridge penalty, apply fisher z-transform and standardize using sample variance. 
	% result.date: 	date,provenance of results
	% 
	% 
	% Examples
	% result = parcar(X,[]) % Use default options


	[m p] = size(X);

	assert(length(size(X))==2,'X has wrong dimensions');
	assert(m>2*p,'Sample size is inadequate for this function')
	try
		assert(abs(sum(mean(X(:,:,1))))<1e-10,'Not zero-mean centered. Please center time-series')
	catch
		X = bsxfun(@minus,X,mean(X,1));
		assert(abs(sum(mean(X(:,:,1))))<1e-10,'Unsuccessful centering time-series')	
	end

	Sighat = cov(X);
	[V D] = eig(Sighat);
	
	try
		assert(sum(diag(D)<=0)==0,'Zero or Negative Eigenvalues');
	catch
		disp('Using ridge penalization. Covariance matrix is not positive definite');
	end
	
	if(isempty(opts)|isempty(opts.lambda))	
		if(min(diag(D))>0)
			if(m>=3*p)
				opts.lambda = 0;
			elseif(m>=2*p)
				opts.lambda = .01;
			else
				opts.lambda = .05;
			end
		else
			opts.lambda = min(diag(D))+.1;
		end
	end
		
	useRidge = (opts.lambda~=0);
	
	if(~useRidge)
		Theta = inv(Sighat);
	else
		Theta = inv(Sighat+opts.lambda*eye(p));
	end
	
	diagTheta = diag(diag(Theta));
	Rho = -inv(sqrt(diagTheta))*triu(Theta,1)*inv(sqrt(diagTheta));
	Rho = eye(p)+Rho+Rho';
	
	% Z-transformed and 
	if(~useRidge)
		df = m-3-p;
		samplevar = 1/sqrt(df);
		Zstat = triu(.5*log((1+Rho)./(1-Rho))./samplevar,1);
	end
	
	result.Rho = Rho;	
	result.Theta = Theta;	
	result.Var = diag(diagTheta); 
	result.opts = opts;
	if(useRidge)
		result.dfe = 'NA';
		result.Zstat = [];
	else
		result.dfe = df;	
		result.Zstat = Zstat;
	end
	result.date = datestr(now);
	
	if(~isempty(opts.visible))		
		if(opts.visible)
			figure('visible','on')
			subplot(1,2,1);
			imagesc(Rho); axis image; title('Partial Correlation Matrix');
			if(~useRidge)
				subplot(1,2,2);
				imagesc(Zstat); axis image; title('Z-Test Statistic Map');
			end
		end
	end
		
end