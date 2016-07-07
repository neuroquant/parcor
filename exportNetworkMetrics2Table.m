function [metrictbl] = exportNetworkMetrics2Table(metrics, condition)
% exportNetworkMetrics2Table
%
% Input
% metrics: At most metrics data has 4 n features x m metrics x s subjects x b resamples
% metricstype: edgel, nodal, global
% metricsnames: Name of metric
% condition: a struct containing fieldnames and values that will be added to the table, and remain constant across all
% observations/rows of the table. (This exists to aid in creating a global table of metrics across conditions)
%
% Output
% - table of observations x  (metrics , subject_id, resample_no)
%
% Description:
% Exports results of pcnets or any other package into tabular form suitable for import into statistical inference packages
% For example global metrics have n=1 feature while nodal metrics have n=#nodes
%
% function [table] = exportNetworkMetrics2Table(metrics)

% Loop through the following
% length(condition.SubjectID)
% length(condition.Resample)

metrictbl = exportTableHelper(metrics,condition);

end



function [metrictbl] = exportTableHelper(metrics,condition)
  % Input
  % - metrics: Has at most 4 n features  x s subjects x b resamples for a given metric in metricfield
  %
  %
  % Output
  % - table of observations x  ($metricfield , subject_id, resample_no)
  %

% if more than 1 resample per subject
% Helps unwrap subjects x resamples into table format
metricfields = fieldnames(metrics);
n_metrics = length(metrics);
metrictbl = table();

if(sum(strcmp(fieldnames(condition),'idxThreshold')))
  if(~isempty(getfield(condition,'idxThreshold')))
    thresh_idx = getfield(condition,'idxThreshold');
  else
    thresh_idx = 1;
  end
end

for mm=1:n_metrics
  % Find base metric name. example: EigenvectorCentrality
  if(sum(strcmp(metricfields,'name')))
    basemetric = getfield(metrics(mm),'name');
    otherfields = setdiff(metricfields,'name');
  else
    basemetric= ['metric_' num2str(mm)];
    otherfields = metricfields;
  end

  % Store all other fieldnames and values in table according to EigenvectorCentrality_<fieldname>
  for ff=1:length(otherfields)
    if(~isempty(getfield(metrics(mm),otherfields{ff})))
      tmpVariableName = {strcat(basemetric, '_', otherfields{ff})};
      tmpValue = getfield(metrics(mm),otherfields{ff});
      % Check if nodal metrics have 2 dimensions.
      if(sum(size(tmpValue)>1)>=2)
        % Store each node's metric as a variable
        for nn=1:size(tmpValue)
          tmpVariableName2= {strcat(tmpVariableName{1},'_',num2str(nn))};
          metrictbl = setfield(metrictbl,tmpVariableName2{1},tmpValue(nn,thresh_idx));
        end
        % Alternatively, store as matrix
        %tmpValue = reshape(tmpValue, [1 size(tmpValue)]);
      else
        metrictbl = setfield(metrictbl,tmpVariableName{1},tmpValue(thresh_idx));
      end
    end
  end


end

% Specify name of the subject, resample and condition, common to all metrics.
metrictbl = setfield(metrictbl,'SubjectID',{getfield(condition,'SubjectID')});
metrictbl = setfield(metrictbl, 'StimLabel',{getfield(condition,'StimLabel')});
metrictbl = setfield(metrictbl,'Resample',getfield(condition,'Resample'));
metrictbl = setfield(metrictbl,'gitcommit',getfield(condition,'gitcommit'));


end
