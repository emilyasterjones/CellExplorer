function cell_metrics = FR2017(cell_metrics,session,spikes,parameters)
    % This is an example template for creating your own calculations
    %
    % INPUTS
    % cell_metrics - cell_metrics struct
    % session - session struct with session-level metadata
    % spikes_intervalsExcluded - spikes struct filtered by (manipulation) intervals
    % spikes - spikes cell struct
    %   spikes{1} : all spikes
    %   spikes{2} : spikes excluding manipulation intervals
    % parameters - input parameters to ProcessCellExplorer
    %
    % OUTPUT
    % cell_metrics - updated cell_metrics struct
    
    sr = session.extracellular.sr;
    spikes = spikes{1};
    if any(spikes.total==0)
        cell_indexes = find(spikes.total>0);
    else
        cell_indexes = 1:spikes.numcells;
    end
    
%     acg_mean = nan(1,spikes.numcells);
    burst_index = nan(1,spikes.numcells);
    refrac = nan(1,spikes.numcells);
    for i = cell_indexes
        acg = CCG(spikes.times{i},ones(size(spikes.times{i})),'binSize',0.001,'duration',.1,'norm','rate','Fs',1/sr);
        acg = acg(51:end); % subset to after spike
%         acg_mean(i) = mean(acg);
        if max(acg(1:10))==0
            burst_index(i) = 0;
        else
            burst_index(i) = (max(acg(1:10)) - mean(acg(end-10:end)))/max(acg(1:10));
        end
        
        [~, max_ind] = max(acg);
        for m = 2:max_ind
          slope(m) = (acg(m) - acg(1))/(m - 1);
        end
        slope_sd = std(slope);
        exceed_1SD = find(slope > slope_sd);
        refrac(i) = exceed_1SD(1) - 1;
    end
    
%     cell_metrics.acg_mean_FR = acg_mean;
    cell_metrics.burst_index_FR = burst_index;
    cell_metrics.refrac_period_FR = refrac;
end