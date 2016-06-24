saveFiles = 0;
addGlobal = 0;
if(addGlobal)
  % Add global network metrics from nodal ones %
  bw_idx = strfind(tableMetrics.Properties.VariableNames,'RandomWalkBetweenness_nodal');
  bw_idx = find(~cellfun(@isempty,bw_idx));
  bw_val = mean(table2array(tableMetrics(:,bw_idx)),2);
  cl_idx = strfind(tableMetrics.Properties.VariableNames,'RandomWalkCloseness_nodal');
  cl_idx = find(~cellfun(@isempty,cl_idx));
  cl_val = mean(table2array(tableMetrics(:,cl_idx)),2)
  eff_idx = strfind(tableMetrics.Properties.VariableNames,'WeightedEfficiency_nodal');
  eff_idx = find(~cellfun(@isempty,eff_idx));
  eff_val = mean(table2array(tableMetrics(:,eff_idx)),2)
  eig_idx = strfind(tableMetrics.Properties.VariableNames,'EigenvectorCentrality_nodal');
  eig_idx = find(~cellfun(@isempty,eig_idx));
  eig_val = mean(table2array(tableMetrics(:,eig_idx)),2)
  tmp_table = table( bw_val, cl_val, eff_val, eig_val);
  tmp_table.Properties.VariableNames = {'RandomWalkBetweenness_global','RandomWalkCloseness_global','WeightedEfficiency_global','EigenvectorCentrality_global'};
  tableMetrics = horzcat(tableMetrics,tmp_table);
  clear tmp_table;
end
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Create Design Matrix %%%%%%
nodal_no = [];
ytitle = 'Centralization'
%ytitle = [roinames(nodal_no).name ' Nodal Centrality Metrics'];
figuretitle = 'Centralization changes by Site';
subjects = unique(tableMetrics.SubjectID);
conditions = unique(tableMetrics.StimLabel); conditions{5} = ['_' conditions{5}];
metricOI = {['WeightedEfficiency_centralization'],'StimLabel'};
% metricOI = {['WeightedEfficiency_nodal_' num2str(nodal_no)],'StimLabel'};
tableDiffMetrics = table();
%contrasts = [1 0 0 0 -1; 0 1 0 0 -1; 0 0 1 0 -1; 0 0 0 1 -1; 1 -1 0 0 0; 0 0 1 -1 0];
contrasts = [1 0 0 0 -1; 0 1 0 0 -1; 0 0 1 0 -1; 0 0 0 1 -1];
for contrastno = 1:size(contrasts,1)
  contrastLabel = strcat(conditions{find(contrasts(contrastno,:))})
  tmpMetric = []; tmpSubjects = {};
  for ss=1:length(subjects)
          subj_metric = tableMetrics(strcmp(tableMetrics.SubjectID,subjects{ss}),metricOI);
          if(height(subj_metric)==size(contrasts,2))
              tmpMetric(ss) = contrasts(contrastno,:)*table2array(subj_metric(:,1));
              tmpSubjects(ss) = {subjects{ss}};
          end
  end
  if(contrastno==1)
    tmptable = cell2table(tmpSubjects'); tmptable.Properties.VariableNames{1} = 'SubjectID';
    tableDiffMetrics = setfield(tmptable,contrastLabel,tmpMetric');
  else
    tableDiffMetrics =  setfield(tableDiffMetrics,contrastLabel,tmpMetric');
  end
end
%%%%%%% Plot Results %%%%%%%%%%%
h = figure('position', [100 100 1200 900]); set(gcf,'Renderer','Painters');
%%%%%%%%%%%%%%%%%%%%%%%%
fontsz = 24;
subj_NA = find(sum(table2array(tableDiffMetrics(:,2:end)),2)~=0);
ds = table2dataset(tableDiffMetrics(subj_NA,:));

contrastLabels = ds.Properties.VarNames(2:end);

diff_lme = {}; x = []; y = []; group = {};
for cc=1:length(contrastLabels)
    diff_lme{cc} = fitlme(ds,[contrastLabels{cc} ' ~ 1 '])
end

% fitlme(ds,[contrastLabels{1} ' ~  1 + ' contrastLabels{2}])
% fitlme(ds,[contrastLabels{3} ' ~  1 + ' contrastLabels{4}])

x = reshape(repmat(contrastLabels,[length(subj_NA) 1]), [length(subj_NA)*length(contrastLabels) 1]);
x_abbrv = regexprep(x,'_Resting','');
y = reshape(table2array(tableDiffMetrics(subj_NA,2:end)), [length(subj_NA)*length(contrastLabels) 1]);
g = gramm('x',regexprep(x,'_Resting',''), ....
                  'y',y, ....;
                  'color',x);
% g.geom_point()
g.stat_summary('geom',{'bar','black_errorbar'}, 'dodge',.5,'width',.6)
% g.geom_jitter('width',0.05,'height',0.05);
% g.geom_line();
g.set_names('x','Stimulation Site','y',ytitle);
g.set_title(figuretitle);
%%%
% Do the actual drawing
g.draw();
%%%%%%%%%%%% Fix Figure %%%%%%%%%%
legh = get(g.legend_axe_handle,'Children');set(legh(1:2:end),'fontsize',.75*fontsz,'fontweight','bold');
tith = get(g.title_axe_handle,'Children'); set(tith,'fontsize',fontsz)
set(g.facet_axes_handles,'fontsize',fontsz,'FontWeight','bold');
ax(1) = get(g.facet_axes_handles,'XLabel'); ax(2) = get(g.facet_axes_handles,'YLabel');
set(ax(1),'FontSize',.8*fontsz,'FontWeight','bold'); set(ax(2),'FontSize',.8*fontsz,'FontWeight','bold')

if(saveFiles)
  export_fig(['Data/' metricOI{1} '_' date], '-png','-transparent','-q101','-depsc','-nocrop','-nofontswap')
  matlab2tikz('filename',['Data/' metricOI{1} '_' date '.tex'], ...
                    'floatFormat', '%.3f','externalData', false, 'standalone', false, ...
                    'height', '.3\textwidth','width', '.3\textwidth', ...
                    'extraTikzpictureOptions',{'baseline','trim axis left', 'trim axis right'});
end
