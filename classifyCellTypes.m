addpath(genpath("C:\Users\Niflheim\Documents\GitHub\External\CellExplorer"))

%% Build dataset of all identified monosynaptic connections
% initialize
troughToPeak = [];
ab_ratio = [];
monoSynType = {};

% set session
sessions = readtable('D:\WT_Sequences\all_sessions.csv');
probe = 0;
for s = 1:height(sessions)
    task = string(sessions{s,'Task'});
    rec_error = string(sessions{s,'Recording_Error'});
    if strcmp(task, 'X Maze') && strcmp(rec_error, 'FALSE')

        % generate the path to the directory containing the ap.bin file
        base_dir = string(sessions{s,'Base_Directory'});
        base_dir = strrep(base_dir, 'Z', 'D'); %while running off local drive
        base_dir = strrep(base_dir, '/', '\');
        ecephys_path = strcat(base_dir, '\Preprocessed_Data\Spikes\g1');
    
        rec_file_stem = split(string(sessions{s,'File'}),'/');
        rec_file_stem = convertStringsToChars(rec_file_stem(2));
        rec_file_path = sprintf('%s\\%s\\Ecephys\\%s\\catgt_%s\\%s_imec%d',...
                                ecephys_path, string(sessions{s,'Animal'}),...
                                rec_file_stem(1:end-3), rec_file_stem,...
                                rec_file_stem, probe);
        cd(rec_file_path)
        
        % load session
        session = loadSession;
        cell_metrics = loadCellMetrics;

        % find putative monosynaptic connections
        exc_ind = sort(unique(cell_metrics.putativeConnections.excitatory(:,1)));
    	inh_ind = sort(unique(cell_metrics.putativeConnections.inhibitory(:,1)));
        % remove overlaps from excitatory
        [~,overlap_exc, overlap_inh] = intersect(exc_ind, inh_ind);
        exc_ind = exc_ind(~ismember(1:length(exc_ind),overlap_exc));
        inh_ind = inh_ind(~ismember(1:length(inh_ind),overlap_inh));

        % add to variables
        troughToPeak = [troughToPeak, cell_metrics.troughToPeak(exc_ind), cell_metrics.troughToPeak(inh_ind)];
        ab_ratio = [ab_ratio, cell_metrics.ab_ratio(exc_ind), cell_metrics.ab_ratio(inh_ind)];
        for i = 1:length(exc_ind)
            monoSynType{end+1} = 'Excitatory';
        end
        for i = 1:length(inh_ind)
            monoSynType{end+1} = 'Inhibitory';
        end
    end
end

%% Train model
% fit an SVM
svm_model = fitcsvm([troughToPeak', ab_ratio'], monoSynType');
save('D:\WT_Sequences\Analysis\MEC_cell_type_SVM_2021+2022_XMazeonly.mat', 'svm_model')
% Linear predictor coefficients
beta = svm_model.Beta; 
bias = svm_model.Bias;

% plot
figure
hold on
h = gscatter(troughToPeak',ab_ratio',monoSynType',['r','c'],['o','o'],[4,4]);
%set(h(1), 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5);
%set(h(2), 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.5);
X1 = linspace(min(troughToPeak),max(troughToPeak),100);
X2 = -(beta(1)/beta(2)*X1)-bias/beta(2);
plot(X1,X2,'-')
xlim([min(troughToPeak),max(troughToPeak)])
ylim([min(ab_ratio),max(ab_ratio)])

%% For each session, classify neurons
% set session
sessions = readtable('D:\WT_Sequences\all_sessions.csv');
probe = 0;
for s = 22%1:height(sessions)
    task = string(sessions{s,'Task'});
    rec_error = string(sessions{s,'Recording_Error'});
    if strcmp(task, 'X Maze') && strcmp(rec_error, 'FALSE')

        % generate the path to the directory containing the ap.bin file
        base_dir = string(sessions{s,'Base_Directory'});
        base_dir = strrep(base_dir, 'Z', 'D'); %while running off local drive
        base_dir = strrep(base_dir, '/', '\');
        ecephys_path = strcat(base_dir, '\Preprocessed_Data\Spikes\g1');
    
        rec_file_stem = split(string(sessions{s,'File'}),'/');
        rec_file_stem = convertStringsToChars(rec_file_stem(2));
        rec_file_path = sprintf('%s\\%s\\Ecephys\\%s\\catgt_%s\\%s_imec%d',...
                                ecephys_path, string(sessions{s,'Animal'}),...
                                rec_file_stem(1:end-3), rec_file_stem,...
                                rec_file_stem, probe);
        cd(rec_file_path)
        
        % load session
        session = loadSession;
        cell_metrics = loadCellMetrics;

        % classify cells by which side of the SVM line they are on
        which_side = beta(1)*cell_metrics.troughToPeak + beta(2)*cell_metrics.ab_ratio + bias;
        cell_metrics.MizusekiType = cell(size(cell_metrics.putativeCellType));
        [cell_metrics.MizusekiType{which_side<0}] = deal('Excitatory');
        [cell_metrics.MizusekiType{which_side>0}] = deal('Inhibitory');
        % save to file to read into Python
        cellTypeFileName = sprintf('%s//%s_tcat.imec%d.ap.MizusekiType.csv', rec_file_path, rec_file_stem, probe);
        cellTypeTable = cell2table(cell_metrics.MizusekiType', 'VariableNames', ["cell_type"]);
        writetable(cellTypeTable,cellTypeFileName)

        % add putative monosynaptic connection types
        % find putative monosynaptic connections
        exc_ind = sort(unique(cell_metrics.putativeConnections.excitatory(:,1)));
    	inh_ind = sort(unique(cell_metrics.putativeConnections.inhibitory(:,1)));
        % remove overlaps
        [~,overlap_exc, overlap_inh] = intersect(exc_ind, inh_ind);
        exc_ind = exc_ind(~ismember(1:length(exc_ind),overlap_exc));
        inh_ind = inh_ind(~ismember(1:length(inh_ind),overlap_inh));
        % assign to cells
        cell_metrics.MonoSynType = cell(size(cell_metrics.putativeCellType));
        [cell_metrics.MonoSynType{exc_ind}] = deal('Excitatory');
        [cell_metrics.MonoSynType{inh_ind}] = deal('Inhibitory');
        % save to file to read into Python
        cellTypeFileName = sprintf('%s//%s_tcat.imec%d.ap.MonoSynType.csv', rec_file_path, rec_file_stem, probe);
        cellTypeTable = cell2table(cell_metrics.MonoSynType', 'VariableNames', ["cell_type"]);
        writetable(cellTypeTable,cellTypeFileName)

        % save updated cell metrics to file
        cellMetricsFileName = sprintf('%s//%s_tcat.imec%d.ap.cell_metrics.cellinfo.mat', rec_file_path, rec_file_stem, probe);
        saveCellMetrics(cell_metrics,cellMetricsFileName)
    end
end

%% Archived Code
% % Plot CellExplorer types
% uniqueGroups = unique(cell_metrics.putativeCellType); 
% colors = ['b','c','r','m','y','k','g'];
% figure
% view(3)
% grid on
% hold on
% for k = 1:length(uniqueGroups)
% 	ind = strcmp(uniqueGroups{k},cell_metrics.putativeCellType); 
% 	scatter3(cell_metrics.troughToPeak(ind),cell_metrics.ab_ratio(ind),cell_metrics.acg_tau_rise(ind),100,colors(k),'.'); 
% end
% legend('Narrow Interneuron', 'Narrow Interneuron - MonoSynInh', 'Pyramidal Cell', 'Pyramidal Cell - MonoSynExc', 'Wide Interneuron', 'Wide Interneuron - MonoSynExc', 'Wide Interneuron - MonoSynInh')
% xlabel('troughToPeak')
% ylabel('ab_ratio')
% zlabel('acg_tau_rise')
% 
% % Calculate and plot K means
% X = NaN(length(cell_metrics.ab_ratio),3);
% X(:,1) = cell_metrics.troughToPeak;
% X(:,2) = cell_metrics.ab_ratio;
% X(:,3) = cell_metrics.acg_tau_rise;
% [idx, centroids] = kmeans(X,2,'Start',[0.3, 0.4, 3; 0.6, 0, 0]);
% 
% num_clusters = 2;
% colors = ['b','r'];
% figure
% view(3)
% grid on
% hold on
% for k = 1:num_clusters
% 	ind = (idx==k);
% 	scatter3(cell_metrics.troughToPeak(ind),cell_metrics.ab_ratio(ind),cell_metrics.acg_tau_rise(ind),200,colors(k),'.'); 
% end
% plot(centroids(:,1),centroids(:,2),'kx')
% legend('Inhibitory', 'Excitatory')
% xlabel('troughToPeak')
% ylabel('ab_ratio')
% zlabel('acg_tau_rise')
