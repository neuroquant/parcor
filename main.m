disp('Creating analysis parameters')
create_opts;
disp('Collecting datasets')
create_dataset;
% Add whitening
disp('Remove dependencies')
run_kggm_whitening
% Add resampling
disp('Estimate partial correlation networks')
run_parcor;
disp('Compute network metrics for all partial correlation networks')
compute_network_metrics;

% Export data
% run_export_allconditions

% Plot (exploratory)
% grammplot_metrics

% Conduct inference (TBD)

% Correct for multiple testing (TBD)

% grammplot_metrics

% exit;
