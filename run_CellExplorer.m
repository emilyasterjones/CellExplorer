addpath(genpath("C:\Users\Niflheim\Documents\GitHub\External\CellExplorer"))

% set session
sessions = readtable('D:\WT_Sequences\all_sessions.csv');
start_gate = 1;
probe = 0;
for s = 1:height(sessions)
    task = string(sessions{s,'Task'});
    rec_error = string(sessions{s,'Recording_Error'});
    if strcmp(task, 'X Maze') && strcmp(rec_error, 'FALSE')

        % generate the path to the directory containing the ap.bin file
        base_dir = string(sessions{s,'Base_Directory'});
%         base_dir = strrep(base_dir, 'Z', 'D'); %while running off local drive
        base_dir = strrep(base_dir, '/', '\');
        ecephys_path = strcat(base_dir, '\Preprocessed_Data\Spikes\g1');
    
        rec_file_stem = split(string(sessions{s,'File'}),'/');
        rec_file_stem = convertStringsToChars(rec_file_stem(2));
        rec_file_path = sprintf('%s\\%s\\Ecephys\\%s\\catgt_%s\\%s_imec%d',...
                                ecephys_path, string(sessions{s,'Animal'}),...
                                rec_file_stem(1:end-3), rec_file_stem,...
                                rec_file_stem, probe);
        cd(rec_file_path)
        
        % set params
        tic
        session = sessionTemplate(rec_file_path);
        % pass the ap.bin file name (not the path)
        session.extracellular.fileName = sprintf('%s_tcat.imec%d.ap.bin', rec_file_stem, probe);
        
        % generate metrics file
%         session.extracellular.srLfp = 2500; %only necessary if running LFP metrics

        % monosynaptic connections only (fast)
%         cell_metrics = ProcessCellMetrics('session', session,...
%             'metrics',{'monoSynaptic_connections'},...
%             'includeInhibitoryConnections',true,...
%             'manualAdjustMonoSyn',false,...
%             'getWaveformsFromDat',false,...
%             'showWaveforms',false,...
%             'sessionSummaryFigure',false,...
%             'showGUI',false);
        % all interneuron metrics (slow_
        cell_metrics = ProcessCellMetrics('session', session,...
            'metrics',{'monoSynaptic_connections','waveform_metrics','acg_metrics'},...
            'includeInhibitoryConnections',true,...
            'manualAdjustMonoSyn',false,...
            'showWaveforms',false,...
            'sessionSummaryFigure',false,...
            'showGUI',false, ...
            'keepCellClassification',false);
        % display GUI
%         cell_metrics = CellExplorer('metrics',cell_metrics);
        toc
    end
end

%% view from file
% cd(rec_file_path)
% session = loadSession;
% cell_metrics = loadCellMetrics;
% cell_metrics = CellExplorer('metrics',cell_metrics);