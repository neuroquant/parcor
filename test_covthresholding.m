close all; clear;

% load dataset for a single subject
load('aMFG_TimeSeries.mat','Data');
X = Data(:,:,1); 


% Call covariance thresholding on X, with half-half training and test split sizes, and N monte-carlo trials. 
results = covthresholding(X,.3,500);



% Plot loss functions as a function of threshold. 
figure(1); 
fontsz = 20;
loss_types = fieldnames(results.loss);
for ll=1:length(loss_types)
	loss = getfield(results.loss,loss_types{ll});
	plot(results.thresh,loss,'linewidth',3,'color',[.7+.3*(ll/length(loss_types)) 0 .8-.6*(ll/length(loss_types))]);
	hold on;
end
hold off;
ylim([0 .6]);
xlabel('Thresholds','fontsize',fontsz); ylabel('Loss Value','fontsize',fontsz); 
set(gca,'Fontsize',fontsz); 
legend(loss_types); 


% Plot best adjacency matrix (vi thresholded stability) and corresponding correlation values. 
figure(2); 
% Keep all edges that have > .95 stability. 
subplot(1,2,1); title('Stability Thresholded Graph','fontsize',fontsz); 
imagesc(results.Sighat.*(results.Pi(:,:,results.optimal_idx)>.95)); axis image; colorbar;
subplot(1,2,2); title('Correlation Thresholded Graph','fontsize',fontsz); 
imagesc((results.Sighat)); axis image; colorbar;


% Optimally threshold based on stability scores rather than actual correlation values. This guarantees the best consistency across subjects as well. Consistency across subjects is not guaranteed with just correlation thresholds.  


% NOTES: 

% -  Note quadratic loss doesn't work unless large no. of observations are available. 
% - In case of too much variations, Do this at the group level instead. Get estimates of results.loss for all subjects. Pick the minimum across all subjects. 