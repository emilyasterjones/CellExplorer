function [trans, trans_fr_corr, trans_prior_corr] = ce_GetTransProb(rawCCG,n,fr,binSize,conv_win,varargin)
% Extract the baseline corrected CCG + spike trans probe from raw CCG

%rawCCG = spike count between reference and target spike train
% n = number of reference spikes
% bin size = the binning of the CCG (in seconds)
% conv_win = slow, network comodulation time scale (in seconds)
% (optional input) = intwin = time bins in which synapse should inject
% excess synchrony


% define integration window
if ~isempty(varargin)
    intwin = varargin{1};
else
    intwin = round(length(rawCCG)/2) + round([.0008:binSize:.0048]/binSize);
end

%subtract low freq. network comodulation and normalize to # reference spikes
[ ~, pred ] = ce_cch_conv( rawCCG, round(conv_win/binSize));
prob = (rawCCG(:) - pred)/n;
prob(prob<0) = 0;

%integrate
trans = sum(prob(intwin));

% EAJ added 4/19/2023: deal with cases where there are >1 spikes in the
% integration window (so probability can't exceed 1)
trans = min(trans, 1);

%repeat with: subtract baseline FR
prob = rawCCG(:)/n;
trans_fr_corr = sum(prob(intwin)) - fr*(.0048-.0008);
trans_fr_corr = min(trans_fr_corr, 1);

%repeat with: subtract firing rate prior to reference spike
pre_intwin = round(length(rawCCG)/2) - round([.0008:binSize:.0048]/binSize);
trans_prior_corr = sum(prob(intwin)) - sum(prob(pre_intwin));
trans_prior_corr = min(trans_prior_corr, 1);

end