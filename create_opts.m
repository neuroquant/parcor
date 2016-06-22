% Create options structure for this project
%

opts = {};

opts.basepath = {'PNASData/pMFG_bilat_4.14/',...
                          'PNASData/pMFG_left_5.11/',...
                          'PNASData/aMFG_bilat_4.14/', ...
                          'PNASData/aMFG_left_5.11/', ...
                          'PNASResting/'};
opts.inputFiles = {'PNASData/pMFG_bilat_4.14/*pMFG*',...
                          'PNASData/pMFG_left_5.11/*pMFG*',...
                          'PNASData/aMFG_bilat_4.14/*aMFG*', ...
                          'PNASData/aMFG_left_5.11/*aMFG*', ...
                          'PNASResting/*RESTING*'};
opts.conditions = {'R_pMFG','L_pMFG','R_aMFG','L_aMFG','Resting'};
opts.outputFiles = strcat('Data/',opts.conditions, '_', datestr(now,'mmm-dd-yyyy'))
ii=1
basepath = opts.basepath{ii};
listfiles = dir([opts.inputFiles{ii}]);
listfiles = listfiles(find([listfiles.bytes]~=0));
subjID_tmp = cell2table(regexp({listfiles(:).name},'_SP','split')');
subjIDs = subjID_tmp.Var1(:,1);
subjIDs = repmat({subjIDs}, [1 4]);
ii=5
basepath = opts.basepath{ii};
listfiles = dir([opts.inputFiles{ii}]);
listfiles = listfiles(find([listfiles.bytes]~=0));
subjID_tmp = cell2table(regexp({listfiles(:).name},'_RESTING','split')');
subjIDs = cat(2,subjIDs,{subjID_tmp.Var1(:,1)});
opts.subjIDs = subjIDs;
save('pcnets_options','opts')
