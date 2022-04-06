
%%
%data loading
clear all
close all
subj_name = '0026';

protocol_dir = '/home/ksasha/anaconda3/Documents/brainstorm_db/F0026/';
save_dir = ['/home/ksasha/projects/K_waves/results/', subj_name, '/'];

channel_type = 'grad'; % channels you want to analyse ('grad' or 'mag')

pre_T = 20;
post_T = 30;
corr_thresh = 0.85;
signal_dim_thr = 0.95;

%BAND for filtering for localization
f_low = 4;
f_high = 40;


%clustering
Nmin = 8;
thr_dist = 0.02;
draw = 1;

cortex = load(strcat([protocol_dir, 'anat/', subj_name,'/tess_cortex_pial_high.mat']));



channels = load(strcat([protocol_dir, 'data/', subj_name, '/block001/channel_vectorview306_acc1.mat']));

G3_dir=  strcat([protocol_dir, 'data/', subj_name, '/block001/headmodel_surf_os_meg.mat']);


ch_arr = {channels.Channel.Type};
ch_idx = boolean(zeros(1,length(ch_arr)));

% channels: gradiometers or magnetometers
if strcmp(channel_type, 'grad') == 1
    full_ch_type = 'MEG GRAD';
elseif strcmp(channel_type, 'mag') == 1
    full_ch_type = 'MEG MAG';
end
for i = 1:length(ch_arr)
    ch_idx(i) = strcmp(ch_arr{i}, full_ch_type);
end
ch_idx = find(ch_idx);


Data = load(strcat([protocol_dir, 'data/', subj_name, '/block001/data_210511_0421_copy.mat']));

channel_idx = ch_idx(Data.ChannelFlag(ch_idx)~=-1);
Fs = 1/(Data.Time(2)-Data.Time(1));


% Evaluate until this point to plot the waves (plot_wave.m)

%%



%base_sp_ind stores indices of manually marked spikes (all events with
%names, differerent from  "BAD"
kwave_ind = zeros(1,0);
if length(Data.Events) > 0 
    addit = Data.Time(1);
    for i = 1:size(Data.Events,2)
        if ~strcmp(Data.Events(i).label,'BAD')
            for j = 1:length(Data.Events(i).times)
                kwave_ind = [kwave_ind, int32((Data.Events(i).times(1,j)-addit)*Fs)+1];
            end
        end
    end
end

%%
%DATA filtering






% 3. MUSIC dipole fitting


[IndMax, ValMax, ind_m] = music(kwave_ind, Data.F(channel_idx,:), G3_dir, channel_idx, ...
    f_low, f_high, corr_thresh, Fs, pre_T, post_T, signal_dim_thr);



cluster_ind = clustering(Nmin, IndMax, thr_dist, draw, cortex);

%reduction of non-clustered events
ind_m(cluster_ind == 0) = [];
IndMax(cluster_ind == 0) = [];
cluster_ind(cluster_ind == 0) =[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save('/home/ksasha/projects/adapted_TW/results/jkan_first_record','IndMax','ValMax','ind_m','spikeind');
% load('/home/ksasha/projects/adapted_TW/results/jkan_first_record');

%IndMax - location
%ValMax - subspace corr
%ind_m - indices of events
save([save_dir, 'time_and_sources.mat'],'IndMax', 'ValMax', 'ind_m','cluster_ind','kwave_ind');