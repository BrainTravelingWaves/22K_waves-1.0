function [IndMax, ValMax, ind_m, events_idx] =music(kwave_ind, Data, G3_dir, channel_idx, ...
     f_low, f_high, corr_thresh, Fs, pre_T, post_T, signal_dim_thr)

% -------------------------------------------------------------------------
% Spike localization with RAP-MUSIC
% -------------------------------------------------------------------------
% INPUTS:
%   spike_ind -- time stamps with spiking events
%   Data -- brainstorm structure with artifact corrected maxfiltered MEG data
%   G3 -- brainstorm structure with forward operator
%   channel_type -- channels used ('grad' or 'mag')
%   f_low, f_high -- the bands for prefiltering before fitting
%   spikydata -- indicatior, showing whether you want to fit
%                   dipoles for the raw data (0) or to the ICA composed data (1)  
%   picked_components -- if spikydata == 1, the picked ICA components
%                           timeseries
%   picked_comp_top -- if spikydata == 1, the picked ICA components
%                           topographies
%   corr_thresh -- subcorr threshold level
%   RAP -- 'RAP' to run RAP-MUSIC and something else for fast one-round
%       MUSIC
%   
% OUTPUTS:
%   IndMax -- locations of fitted dipoles
%   ValMax -- values of subspace correlation
%   ind_m -- indices of spikes survived the subcorr threshold
%   spikeind -- time indices of spikes in RAP-MUSIC
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

G3 = load(G3_dir);

% 2D forward operator
[G2, ~] = G3toG2(G3, channel_idx); 
clear G3

[b,a] = butter(3, [f_low f_high]/(Fs/2)); 
Ff = Data;
clear Data
for nch = 1:size(Ff,1)
    Ff(nch,:) = filtfilt(b, a, Ff(nch,:)')';
end


[b,a] = butter(3, [48 52]/(Fs/2), 'stop'); 
for nch = 1:size(Ff,1)
    Ff(nch,:) = filtfilt(b, a, Ff(nch,:)')';
end
[b,a] = butter(3, [96 104]/(Fs/2), 'stop'); 
for nch = 1:size(Ff,1)
    Ff(nch,:) = filtfilt(b, a, Ff(nch,:)')';
end


ValMax = [];
IndMax = [];

T = size(Ff,2);
for j = 1:length(kwave_ind((kwave_ind<T-post_T)&(kwave_ind>pre_T)))   
     kwave = Ff(:,(kwave_ind(j)-pre_T):(kwave_ind(j)+post_T));
     [U,S,V] = svd(kwave); %maybe do it econ
     h = cumsum(diag(S)/sum(diag(S)));
     
     %dimension of signal space
     n = find(h>=signal_dim_thr);
     corr = MUSIC_scan(G2, U(:,1:n(1)));
     [ValMax(j), IndMax(j)] = max(corr);
        j
end


figure
histogram(ValMax)

ind_m = find((ValMax > corr_thresh));
IndMax = IndMax(ind_m);
events_idx = kwave_ind(ind_m);

disp(['Subcorr threshold: ', num2str(corr_thresh), ' Number of spike found: ', ...
    num2str(length(ind_m))]);

end