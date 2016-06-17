load('pcnets_options.mat')

for ii=1:length(opts.conditions)
	filename = opts.outputFiles{ii};
	load(filename,'Data');

	Data = Data(:,:,find([squeeze(sum(sum(Data,2),1))]));

  	tmpKGGM = KGGM(Data);
	tmpKGGM.Data = Data;
  	tmpKGGM.verbose=0;
	
	% Block-Banding Matrix
	t = size(Data,1);
	TRlen = 8;
	numBlocks = 10;
	nlags = 1;

	if(strcmp(opts.conditions(ii),'Resting'))
		
	  	tmpKGGM.fit_banding(nlags,1);
	else	
		inv_supp = zeros(t,t);
		inv_blk = zeros(TRlen,TRlen);
		for l=1:nlags
			if (mod(l,2)==1)
				 inv_blk = inv_blk +  ...
				 				(-1)^(l)*(triu(ones(TRlen,TRlen),l)  + tril(ones(TRlen,TRlen),-l)) +  ...
								(-1)^(l+1)*(triu(ones(TRlen,TRlen),l+1) + tril(ones(TRlen,TRlen),-(l+1)));

			else
				inv_blk = inv_blk + (-1)^(l-1)*(triu(ones(TRlen,TRlen),l)  + tril(ones(TRlen,TRlen),-(l))) + ...
								(-1)^(l)*(triu(ones(TRlen,TRlen),l+1)  + tril(ones(TRlen,TRlen),-(l+1)));

			end
		end
		inv_sup = repmat(inv_blk,[numBlocks numBlocks]);

	  	tmpKGGM.fit_banding(nlags,1,inv_sup);	
	end
	
	Omeghat = tmpKGGM.Theta1;
	[V D] = eig(Omeghat); 
	% Omeghat = VDV';
	sqrtOmeg = V*sqrt(D)*V'; 
	
	WhData = zeros(size(Data)); 
	for cc=1:size(Data,3)
		% Cholesky whitening of temporal observations
		% Replace X(i) with Omeg^(-1/2)X(i)
		WhData(:,:,cc) = sqrtOmeg*Data(:,:,cc); 
	end

	if(tmpKGGM.verbose)
		tmpKGGM2 = KGGM(WhData);
		tmpKGGM2.Data = WhData;
		tmpKGGM2.verbose = 0
		tmpKGGM2.fit_banding(1,1);
		figure;
		subplot(1,2,1); imagesc(tmpKGGM.Sigma1); axis image;
		subplot(1,2,2); imagesc(tmpKGGM2.Sigma1); axis image;
		figure;
		subplot(1,2,1); imagesc(corr(Data(:,:,1))); axis image;
		subplot(1,2,2); imagesc(corr(WhData(:,:,1))); axis image;
	end		
	save([filename '.mat'],'WhData','-append');
end
