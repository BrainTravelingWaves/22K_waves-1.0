
load([save_dir, 'time_and_sources.mat'])

Vertices_sm = tess_smooth_(cortex.Vertices, 1, 500, cortex.VertConn, 1);



%ValMax % - number of vertex
%kwave_ind % - index of  sample

nhb_degree = 1; %1 - first order neighbours
f_low = 4;
f_high = 40;
PARAMS = struct();
PARAMS.draw_paths = 0;
PARAMS.draw_wave = 1;
PARAMS.max_distance = 0.04; 
PARAMS.wave.duration = 0.02; % (s) time of wave duration
PARAMS.nspoints = 30; % number of spatial points in wave
PARAMS.wave.speeds = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]; % (mm/ms = m/s)
PARAMS.sampling_rate = Fs;
%faces = cortex.Faces;
vertices = cortex.Vertices;


%% Forward model matrix
G = gain_orient(G3_dir,channel_idx); % TO DO: write it without braintorm function


clear Ndir
clear dirstore
PARAMS.VertConn = cortex.VertConn;
wavelin_cluster = [];
dirstore = [];
cumdir = 1;
n_strts = [];


path = {};
path_sm= {};
cum_n= 1;
for npoint = 1:length(IndMax)
    central_strt = IndMax(npoint);
    cum_nhb = sparse(1,length(vertices));
    cum_nhb(1,central_strt) = 1;
    cum_nhb = cum_nhb>0;
    for n_nhb = 1:nhb_degree
        cum_nhb = cum_nhb*PARAMS.VertConn+cum_nhb;
        cum_nhb = cum_nhb > 0;
    end
    find_nhb = find(cum_nhb);
    %find new vertices
    
    local_pathes = {};
    local_pathes_sm = {};
     for n_start = 1:length(find_nhb)



        PARAMS.seed_vertices = find_nhb(n_start);
        [waves,directions,PATH, PATH_sm] = wave_on_sensors_simple(cortex, PARAMS, G, Vertices_sm);

        wavelin = [];
        Ndir(cum_n) = size(waves,1);
        
        substrat = cell2mat(PATH);
        buf = zeros([size(PATH),3]);
        buf(:,:,:,1) = substrat(:,1:3:90,:);
        buf(:,:,:,2) = substrat(:,2:3:90,:);
        buf(:,:,:,3) = substrat(:,3:3:90,:);
        
        substrat_sm = cell2mat(PATH_sm);
        buf_sm = zeros([size(PATH_sm),3]);
        buf_sm(:,:,:,1) = substrat_sm(:,1:3:90,:);
        buf_sm(:,:,:,2) = substrat_sm(:,2:3:90,:);
        buf_sm(:,:,:,3) = substrat_sm(:,3:3:90,:);
        
        local_pathes_sm{n_start} = buf_sm;
        local_pathes{n_start} = buf;


        dirstore = [dirstore, directions];
        for s = 1:length(PARAMS.wave.speeds)
            wavesspeed = waves(:,s);%wave of certain speed
            clear wavesspeedlin
            for d = 1:Ndir(cum_n) % for each direction
                tmp = wavesspeed{d}; % matrix (num sensors) by (num time points)
                wavesspeedlin(d,:) = tmp(:); % (num dir) by (num sensors)x(num time points)
            end
            wavelin = [wavelin, wavesspeedlin];
        end
        wavelin_cluster = [wavelin_cluster; wavelin];%Nspikes * Ndirineachspike x nT * Nsensors

        cumdir = cumdir + Ndir(cum_n);
        cum_n = cum_n+1;
     end
    path{npoint} = local_pathes;
    path_sm{npoint} = local_pathes_sm;
    n_strts = [n_strts, length(find_nhb)];
    npoint
end

csvwrite(join([save_dir 'directions.csv']), dirstore)%cl5
csvwrite(join([save_dir 'wave.csv']), wavelin_cluster)
csvwrite(join([save_dir 'ndir.csv']), Ndir) %info about num of directions
csvwrite(join([save_dir 'nstrts.csv']),n_strts)
save([save_dir 'path.mat'],'path')
save([save_dir 'path_sm.mat'],'path_sm')
clear wavelin_cluster



[b,a] = butter(3, [f_low f_high]/(Fs/2)); % butterworth filter before ICA
Ff = filtfilt(b, a, Data.F(channel_idx,:)')';


[b,a] = butter(3, [48 52]/(Fs/2), 'stop'); % butterworth filter before ICA
Ff = filtfilt(b, a, Ff')';
[b,a] = butter(3, [96 104]/(Fs/2), 'stop'); % butterworth filter before ICA
Ff = filtfilt(b, a, Ff')';


DataLin = [];
wavelength = size(waves{1, 1}, 2);
k = 1;
%number of slides
R = post_T+pre_T - size(waves{1, 1}, 2);
for n_point = 1:length(kwave_ind(ind_m))
    spike_ts = Ff(:,(kwave_ind(ind_m(n_point))-pre_T):(kwave_ind(ind_m(n_point))+post_T));
    range = 1:wavelength;
    for r = 1:R % sliding time window
        tmp = spike_ts(:, range); % time interval with wave
        DataLin(:, k) = tmp(:);
        range = range+1;
        k = k+1;
    end
end



%datalin - wavedur *nsensors x Nspikes*
csvwrite(join([save_dir 'spike.csv']), DataLin)










f_low = 1;
f_high = 40;
[b,a] = butter(3, [f_low f_high]/(Fs/2)); % butterworth filter before ICA
Ff = filtfilt(b, a, Data.F(channel_idx,:)')';


[b,a] = butter(3, [48 52]/(Fs/2), 'stop'); % butterworth filter before ICA
Ff = filtfilt(b, a, Ff')';
[b,a] = butter(3, [96 104]/(Fs/2), 'stop'); % butterworth filter before ICA
Ff = filtfilt(b, a, Ff')';


DataLinRaw = [];
wavelength = size(waves{1, 1}, 2);
k = 1;
%number of slides
R = pre_T+post_T - size(waves{1, 1}, 2);
for n_point = 1:length(kwave_ind(ind_m))
    spike_ts = Ff(:,(kwave_ind(ind_m(n_point))-200):(kwave_ind(ind_m(n_point))+200));
  
   
        tmp = spike_ts; % time interval with wave
        DataLinRaw(:, k) = tmp(:);
      
        k = k+1;
end




%datalin - wavedur *nsensors x Nspikes*
csvwrite(join([save_dir 'spike_raw.csv']), DataLinRaw)


