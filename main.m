disp('Creating analysis parameters')
create_opts;
disp('Collecting datasets')
create_dataset;
% Add whitening
% Add resampling
disp('Estimate partial correlation networks')
run_parcor;
disp('Compute network metrics for all partial correlation networks')
compute_network_metrics;

% exit;
