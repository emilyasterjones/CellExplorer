%% set session
animal = 'BaggySweatpants';
session_name = '20220214_BaggySweatpants_DY01';
basepath = ['Z:\WT_Sequences\2022_winter\Preprocessed_Data\Spikes\g1\',animal,'\Ecephys\',...
    session_name,'\catgt_',session_name,'_g1\',session_name,'_g1_imec0'];
cd(basepath)

%% calculate
tic
% % set params
session = sessionTemplate(basepath);
session.extracellular.fileName = [session_name,'_g1_tcat.imec0.ap.bin'];
session.extracellular.srLfp = 2500;
% session = gui_session(session);
% validateSessionStruct(session);

% generate metrics file
cell_metrics = ProcessCellMetrics('session', session,...
    'metrics',{'monoSynaptic_connections'},...
    'includeInhibitoryConnections',true,...
    'manualAdjustMonoSyn',false,...
    'getWaveformsFromDat',false,...
    'showWaveforms',false,...
    'sessionSummaryFigure',false,...
    'showGUI',false);
toc

% %% view
% cell_metrics = CellExplorer('metrics',cell_metrics);
% 
% %% view from file
% animal = 'BaggySweatpants';
% session_name = '20220214_BaggySweatpants_DY01';
% basepath = ['Z:\WT_Sequences\2022_winter\Preprocessed_Data\Spikes\g1\',animal,'\Ecephys\',...
%     session_name,'\catgt_',session_name,'_g1\',session_name,'_g1_imec0'];
% cd(basepath)
% session = loadSession;
% cell_metrics = loadCellMetrics;
% cell_metrics = CellExplorer('metrics',cell_metrics);