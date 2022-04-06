close all

%parameters of the wave from "create_basic_waves.m" and "wave_localization"
ntpoints = 20;
nspoints = 30;
pre_t = 20;

%function plot_wave(V,F,pathes,activation,coefs)
num_cl = 1;
%number of wave in cluster to plot (from 1)
num_wave = 3;
nhb_degree = 1;


data_plot_range = [-200,200];
n_channels = 100;

load([save_dir, 'time_and_sources.mat'])
load([save_dir, 'wave_prop_cl', num2str(num_cl), '.mat'])
load([save_dir, 'R2val_cl', num2str(num_cl), '.mat'])
load([save_dir, 'path_sm.mat'])

Vertices_sm = tess_smooth_(cortex.Vertices, 1, 500, cortex.VertConn, 1);



%choose best channels for vizualilzation
channels_locs = zeros(3,length(channel_idx));

for i = 1:length(channel_idx)
    channels_locs(:,i) = mean(channels.Channel(channel_idx(i)).Loc,2);
end







%define strt vertice
cl_idc = find(cluster_ind == num_cl);

central_strt = IndMax(cl_idc(num_wave));
cum_nhb = sparse(1,length(cortex.Vertices));
cum_nhb(1,central_strt) = 1;
cum_nhb = cum_nhb>0;
for n_nhb = 1:nhb_degree
    cum_nhb = cum_nhb*cortex.VertConn+cum_nhb;
    cum_nhb = cum_nhb > 0;
end
find_nhb = find(cum_nhb);

strt_vert = find_nhb(best_strt(num_wave)+1);


ch_distances = zeros(length(channel_idx));
for i = 1:length(channel_idx)
    ch_distances(i) = norm(cortex.Vertices(strt_vert,:)-channels_locs(:,i));
end

[~,sort_channels_idx] = sort(ch_distances);
best_ch_idx = sort_channels_idx(1:n_channels);



 
coefs = all_coefs{num_wave};


pathes = squeeze(path_sm{cl_idc(num_wave)}{best_strt(num_wave)+1}(:,:,bestind(num_wave)+1,:));


 t = 0:(ntpoints-1);
 n = 0:nspoints-1;
 
activation = zeros(ntpoints,nspoints);
    
 for i = t
    activation((i+1),:) = (1 + cos(2*pi * (n - i) / nspoints));
 end%wave - t (0 ... T) x 


 
 
[b,a] = butter(3, [2 40]/(Fs/2)); 
Ff = Data.F(channel_idx(best_ch_idx),:);
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
 
 
 
filename = [save_dir, '/vizualization/surf_cl', num2str(num_cl), '_num_wave', num2str(num_wave) ,'.gif'];
 
fig = figure('Position',[100,400,1300,500])%('Renderer', 'painters', 'Position', [700 500 500 300]);
subplot(1,2,1)
event_idx= kwave_ind(ind_m(cl_idc(num_wave)));


plot([0:data_plot_range(2)-data_plot_range(1)],Ff(:,event_idx+data_plot_range(1):event_idx+data_plot_range(2))','k')
hold on
yl = ylim;

event_wave_idx = -data_plot_range(1)- pre_T+int32(bestshifts(num_wave,bestind(num_wave)));
plot([event_wave_idx,event_wave_idx],yl,'r')
hold on
plot([event_wave_idx,event_wave_idx]+ntpoints,yl,'b')
xlabel('time, ms')


subplot(1,2,2)%('Renderer', 'painters', 'Position', [100 100 500 300]);% 1100 800]);
ctoon = patch('Faces',cortex.Faces,'Vertices',Vertices_sm,'FaceColor',[248/255 214/255 197/255],'FaceAlpha',1.0,'EdgeAlpha',0.1);
hold on
%ctoon.FaceVertexCData = vertex_wave;
%ctoon.FaceColor = 'interp';


xlabel('x')
ylabel('y')
zlabel('z')



normal = cortex.VertNormals(strt_vert,:);

view(normal)
set(gca, 'CameraPosition', Vertices_sm(strt_vert,:)*1.2);


t = 1;
for t = 1:size(activation,1)
    subplot(1,2,1)
    %figure(1)
 
  
    plot([0:data_plot_range(2)-data_plot_range(1)],Ff(:,event_idx+data_plot_range(1):event_idx+data_plot_range(2))','k')
  
    hold on
    event_wave_idx = -data_plot_range(1)- pre_T+int32(bestshifts(num_wave,bestind(num_wave)));
    plot([event_wave_idx,event_wave_idx],yl,'r')
    hold on
    plot([event_wave_idx,event_wave_idx]+ntpoints,yl,'b')
    xlabel('time, ms')

    hold on 
    plot([event_wave_idx,event_wave_idx]+t-1,yl,'r--')
    
     



    
    subplot(1,2,2)
    %figure(2)
    for i = 1:length(coefs)
    if coefs(i) >0
    patch(pathes(i,:,1),pathes(i,:,2),pathes(i,:,3),coefs(i)*activation(t,:),...
                              'EdgeColor','interp','FaceAlpha',0,'Faces',...
                              [1:size(pathes,2)-1;2:size(pathes,2)]','LineWidth',3);
                          

    axis equal
    end
    end
    drawnow
    frame = getframe(fig);
    [A,map] = rgb2ind(frame2im(frame),256);
    if t == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2)%,'ScreenSize',[1000 700]);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2)%,'ScreenSize',[1000 700]);
    end
    
end
%xlim([min(pathes(:,:,1),[],'all') max(pathes(:,:,1),[],'all')])
%ylim([min(pathes(:,:,2),[],'all') max(pathes(:,:,2),[],'all')])
%zlim([min(pathes(:,:,3),[],'all') max(pathes(:,:,3),[],'all')])



%figure
%plot3(pathes(1,1:end-1,1),pathes(1,1:end-1,2),pathes(1,1:end-1,3))



%%\


