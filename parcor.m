function result = parcor(X,opts)
	% parcor.m
	% Takes a multivariate time-series data matrix X and computes sample partial correlation matrix
	% INPUT
	%
	% X:		n_samples x n_regions data matrix. Time-series centered to be zero-mean.
	%
	% opts
	% opts.lambda: 	ridge penalty. 0 by default if n_samples > 3*n_regions else set to .1
	% opts.ridgeType: no penalty, ridge type 1 or ridge type 2
	%
	% OUTPUT
	%
	% result
	% result.Rho: 		n_regions x n_regions partial correlation matrix
	% result.Theta:		n_regions x n_regions inverse covariance  matrix
	% result.Var:		n_regions x 1 ; diagonal matrix of variances
	% result.InvVar:	n_regions x 1 ; diagonal matrix of partial variances
	% result.opts:		options parameter
	% result.dfe:		naive degrees of freedom for t-statistic for standard partial correlation. Does not apply if ridge penalty is used.
	% result.Zstat: 	If no ridge penalty, apply fisher z-transform and standardize using sample variance.
	% result.date: 		date,provenance of results
	%
	%
	% Examples
	% result = parcar(X,[]) % Use default options


	[m p] = size(X);

	nan_idx = find(sum(isnan(X))~=0);
	if(~isempty(nan_idx))
		X(:,nan_idx) = 0;
	end

	assert(length(size(X))==2,'X has wrong dimensions');
	assert(m>10*log(p),'Sample size is inadequate for this function')
	try			
		meanX_i = zeros(1,p); 
		assert(abs(sum(mean(X(:,:,1))))<1e-10,'Not zero-mean centered. Please center time-series');
	catch
		disp('Centering time-series')
		if(abs(sum(mean(X(:,:,1))))>1e-10)
			meanX_i = mean(X,1);
			X = bsxfun(@minus,X,mean(X,1));
			assert(abs(sum(mean(X(:,:,1))))<1e-10,'Unsuccessful centering time-series')
			meanX_i = mean(X,1);
		else
			meanX_i = zeros(1,p); 
		end	
	end

	Sighat = cov(X);
	[V D] = eig(Sighat);

	try
		assert(sum(diag(D)<=0)==0,'Zero or Negative Eigenvalues');
	catch
		disp('Covariance matrix is not positive definite');
	end
	
	% Rescale Sighat to be correlation matrix
	diagW = diag(diag(Sighat));
	try
		assert(sum(diag(diagW)<=0)==0,'Zero or Negative Variances');
	catch
		disp('Some Variances Close to 0');
		sum(diag(diagW)<=0)
		diagW = diag(Sighat); 
		diagW(diagW==0) = 1;
		diagW = diag(diagW);
	end
	Sighat = inv(sqrt(diagW))*Sighat*inv(sqrt(diagW));
	Xcorr = bsxfun(@rdivide,X,diag(sqrt(diagW))');




	if(isempty(opts))
		Target = [];
	elseif(isfield(opts,'Target'))
		Target = opts.Target;
	else
		switch opts.ridgeType
		case 0
			Target = [];			
		case 1
			Target = diag(diag(Sighat)); % Diagonal matrix with variances from Sighat.
		case 2
			Target = eye(p);
		case 3
			meancorr = sum(sum(triu(Sighat,1)));
			meancorr = meancorr/nchoosek(p,2);
			Target  = eye(p) + ~eye(p)*meancorr;	% Ledoit Wolf with constant-correlation off diagonals 
		end
	end
		

	if(isempty(opts)|isempty(opts.lambda))
		% if(min(diag(D))>.001)
		% 	% Dumb defaults, replace with analytic or CV based thresholds.
		if(m>=10*p)
				opts.lambda = 0;
			% elseif(m>=2*p)
			% 	opts.lambda = .01;
			% else
			% 	opts.lambda = .05;
			% end

		else
			% Dumb defaults, replace with analytic or CV based thresholds.
			opts.lambda = min(diag(D))+.1;

			% Shaffer-Strimmer analytical lambda
			% opts.lambda_ss = sum(sig_var) - cov(target,sigma) - bias(sigma)(target-sigma)/sum((target-sigma).^2);
			% opts.lambda = max(0,min(1,opts.lambda_ss));
			if(opts.ridgeType==1|opts.ridgeType==3)

				% Shaffer-Strimmer, analytical formula for type D
				for ii=1:p
					for jj=ii+1:p
						W_kij(:,ii,jj) = (Xcorr(:,ii)-meanX_i(ii)).*(Xcorr(:,jj)-meanX_i(jj));
						W_kij(:,jj,ii) = W_kij(:,ii,jj);
					end
				end
				if(m>p)
					assert(sum(sum(triu(Sighat-squeeze(mean(W_kij,1)),1).^2))/nchoosek(p,2)< 1e-3,'Verifying mean W_kij is empirical covariance');
				else
					disp('Frob norm error between empirical cov and sample cov')
					sum(sum(triu(Sighat-squeeze(mean(W_kij,1)),1).^2))/nchoosek(p,2)
				end
				var_Sighat = m/(m-1)^3*squeeze(sum(bsxfun(@minus,W_kij,reshape(Sighat, [1 p p])).^2,1));
				opts.lambda_ss = sum(sum(triu(var_Sighat,1)))/sum(sum(triu(Sighat.^2,1)));
				opts.lambda = max(0,min(1,opts.lambda_ss));

			elseif(opts.ridgeType==2)

				% Shaffer-Strimmer, analytical formula for type A
				for ii=1:p
					for jj=ii:p
						W_kij(:,ii,jj) = (Xcorr(:,ii)-meanX_i(ii)).*(Xcorr(:,jj)-meanX_i(jj));
						W_kij(:,jj,ii) = W_kij(:,ii,jj);
					end
				end
				disp('Verifying mean W_kij is empirical covariance');
				assert(sum(sum(triu(Sighat-squeeze(mean(W_kij,1)),1).^2))/nchoosek(p,2)< 1e-3,'Verifying mean W_kij is empirical covariance');
				var_Sighat = m/(m-1)^3*squeeze(sum(bsxfun(@minus,W_kij,reshape(Sighat, [1 p p])).^2,1));
				opts.lambda_ss = sum(sum(var_Sighat))/(2*(sum(sum(triu(Sighat.^2,1))))+ sum((diag(Sighat)-ones(p,1)).^2));
				opts.lambda = max(0,min(1,opts.lambda_ss));

		end
	end


	useRidge = (opts.lambda~=0)*opts.ridgeType;

	% if(~useRidge)
	% 	Theta = inv(Sighat);
	% else
	% 	Theta = inv(Sighat+opts.lambda*eye(p)); % Type 1, Ridge Threshold
	% end

	switch useRidge

	case 0
		%Target = [];
		Sigma = Sighat;
		Theta = pinv(Sigma);
	case 1
		disp('Using Ridge Penalization, Type I')
		%Target = diag(diag(Sighat)); % Diagonal matrix with variances from Sighat.
		Sigma = (1-opts.lambda)*Sighat + opts.lambda*Target;
		Theta = pinv(Sigma); % Type I, Ledoit-Wolf Ridge Penalty.
		
	case 2
		disp('Using Ridge Penalization, Type II')
		%Target = eye(p);
		% Shaffer-Strimmer, analytical formula for type A
		Sigma = Sighat+opts.lambda*eye(p);
		Theta = pinv(Sigma); % Type II, Ridge Penalty
	case 3
		disp('Using Ridge Penalization, Type I w. Constant Correlation Target')
		Sigma = (1-opts.lambda)*Sighat + opts.lambda*Target;
		opts.lambda
		Theta = pinv(Sigma); % Type I, Ledoit-Wolf Ridge Penalty.
	end

	diagTheta = diag(diag(Theta));
	Rho = -inv(sqrt(diagTheta))*triu(Theta,1)*inv(sqrt(diagTheta));
	Rho = eye(p)+Rho+Rho';
	% Rescaling back Theta, since Sighat is correlation matrix
	Theta = inv(sqrt(diagW))*Theta*inv(sqrt(diagW));

	% Z-transformed
	if(~useRidge)
		df = m-3-p;
		samplevar = 1/sqrt(df);
		Zstat = triu(.5*log((1+Rho)./(1-Rho))./samplevar,1);
	end

	result.Rho = Rho;
	result.Sigma = Sigma;
	result.Theta = Theta;
	result.Var = diag(diagW); % Variances
	result.InvVar = diag(diagTheta); % Inverse Variances used to normalize partial correlation matrix.
	result.opts = opts;
	result.Target = Target;
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
