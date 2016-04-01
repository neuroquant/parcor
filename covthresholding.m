% covthresholding.m
% function [results] = CovThresholding(X,n_1,varargin)
%
% Implements Bickel and Levina Covariance Thresholding Procedure along with Stability Selection on Covariance Thresholded Matrices. 
% Split data into training and test splits.  On the training split, threshold a covariance matrix, on the test split estimate a sample covariance matrix. 
% Then compute some loss or error metric comparing the training versus test covariance estimates. Minimize the loss function to obtain the optimal threshold estimate. 
% To obtain a stable network structure (so that the edges eliminated from the graph do not vary in location too much), 
% instead choose threshold value that also minimizes overall graph stability. 
% 
% At this time, does not incorporate whitening of time-dependent observations. Uses block sub-sampling to create training/test splits 
%
% Manjari Narayan 
% Copyright (c) 2016
%  
% INPUTS:
% X: 		Data matrix of t observations x p regions or voxels
% n_1:		specify the size of training split as a proportion. Takes values in (0,1)
% N:		Number of sample splits. Default is 50
% thresh:	vector of thresholds, optional 
% loss_type:	'stability' by default. Optional. Can be frobenius/squared error loss or 'l2'
% 
% OUTPUTS; 
% results.Sighat 			- Correlation Thresholded Estimate
% results.loss  			- vector of Frobenium norm loss values
% results.Wdiag 			- Diagonal Variance. Can be used to convert to Covariance version. as sqrt((Wdiag))*Sighat*sqrt((Wdiag))
% results.Pi 				- Stability of Support of Thresholded Covariance Matrices;
% results.optimal_thresh & optimal_idx 	- Give you optimal value of threshold and corresponding index.

function [results] = CovThresholding(X,n_1,varargin)
	
	% Initialize
	[t p] = size(X);
	
	% Re-scale covariance matrices to correlation matrices. 
	S = cov(X);
	Wdiag = diag(diag(S)); 
	R = sqrt(inv(Wdiag))*S*sqrt(inv(Wdiag));
	
	maxcorr = max(max(triu(R,1)));
	
	switch nargin
	case 2
		N = 50; 
		thresh = linspace(.2*sqrt(log(p)/t),maxcorr,45);
		loss_type = 'stability'
		
	case 3
		N  = varargin{1};
		thresh = linspace(.2*sqrt(log(p)/t),maxcorr,45);
		loss_type = 'stability'
		
	case 4
		N  = varargin{1};
		thresh = varargin{2};
		loss_type = 'stability'
		
	case 5
		N  = varargin{1};
		thresh = varargin{2};
		loss_type = varargin{3}
	end


	
	if(n_1>1)
		error('n_1 needs to be a proportion');
	else
		n_2 = 1-n_1;
	end
	
	n_thresh = length(thresh);
	
	% Block splits (in the absence of whitening), Subsampling
	blen = ceil(1*sqrt(t)); % length of sliding window
	B = floor(t/blen); % subsamples;
	blocks = zeros(B,blen);
	for bb=1:B
		blocks(bb,:) = [(bb-1)*blen+1:(bb-1)*blen+blen];
	end

	
	% % Block splits (in the absence of whitening), Bootstraps
	% blen = ceil(3*sqrt(t)); % length of sliding window
	% B = t-blen+1; % no. of blocks for boostraps
	% blocks = zeros(B,blen);
	% for bb=1:B
	% 	blocks(bb,:) = [1+bb-1:blen+bb-1]; % bootstrap
	% end


	
	
	l2_loss = zeros(n_thresh,N);
	quad_loss = zeros(n_thresh,N);
	stability_loss = zeros(n_thresh,1); cum_stability = stability_loss;
	Pi = zeros(p,p,n_thresh);
	
	for trial_no = 1:N
		for thresh_no = 1:n_thresh
		
			tr_idx = randsample(1:B,ceil(n_1*B),0); 
			ts_idx = setdiff(1:B,tr_idx);
			TrCov = zeros(p,p); 
			TsCov = zeros(p,p); 

		
			% Threshold TrCov
			resamples_idx = reshape(blocks(tr_idx,:),[1 length(tr_idx)*blen]);
			TrCov = cov(X(resamples_idx,:));
			Wd = Wdiag;			
			%Wd = diag(diag(TrCov)); 
			TrCor = sqrt(inv(Wd))*TrCov*sqrt(inv(Wd));
			clear Wd resamples_idx;
			TrCor = TrCor.*(abs(TrCor)>thresh(thresh_no));
			Pi(:,:,thresh_no) = Pi(:,:,thresh_no)+(TrCor~=0);
		
			% Unthresholded Sample Covariance
			resamples_idx = reshape(blocks(ts_idx,:),[1 length(ts_idx)*blen]);
			TsCov = cov(X(resamples_idx,:));
			Wd = Wdiag;			
			%Wd = diag(diag(TsCov)); 
			TsCor = sqrt(inv(Wd))*TsCov*sqrt(inv(Wd));
			clear Wd resamples_idx;
		
			l2_loss(thresh_no,trial_no) = sum(sum((TrCor-TsCor).^2))/p^2;
			if(length(tr_idx)*blen>100*p)
				quad_loss(thresh_no,trial_no) = sum(sum((TsCor*inv(TrCor) - eye(p)).^2))/p^2;	
			end		
		end		
	end
	
	for thresh_no = 1:n_thresh
		tmp_Pi= squeeze(Pi(:,:,thresh_no)/N);
		stability_loss(thresh_no) = sum(sum(triu(2*tmp_Pi.*(1-tmp_Pi),1)))/nchoosek(p,2);
	end
	stability_loss = fliplr(stability_loss);
	cum_stability = [stability_loss(1) max(stability_loss(1:end-1),stability_loss(2:end))'];
	
	
	switch loss_type
		
	case 'l2'
		mean(l2_loss,2)';
		[minval minidx] = min(mean(l2_loss,2));
	case 'quad'
		mean(quad_loss,2)';
		[minval minidx] = min(mean(quad_loss,2));
	case 'stability'
		stable_idx = (find(cum_stability < .1))
		minidx = min(setdiff(stable_idx, [1]));
		minval = stability_loss(minidx);
		
		stability_loss = fliplr(stability_loss);
		minidx = n_thresh-minidx+1;	
		
	end
	
	disp('Best Threshold')
	thresh(minidx)
	
	disp('Min. Loss:')
	minval
	
	R = R.*(abs(R)>thresh(minidx));
	Sighat = R;
	%Sighat = sqrt(Wdiag)*R*sqrt(Wdiag);
	loss.l2 = mean(l2_loss,2);
	loss.quad = mean(quad_loss,2);
	loss.stability = fliplr(cum_stability);
	
	results.Sighat = Sighat;
	results.optimal_thresh = thresh(minidx);
	results.optimal_idx = minidx;
	results.loss = loss;
	results.thresh = thresh;
	results.N = N; 
	results.Wdiag = Wdiag;
	results.Pi = Pi/N;
	