function [cluster] = clustering(Nmin, src_idx, thr_dist, draw, cortex)

% -------------------------------------------------------------------------
% Spatial clustering of localized spikes
% -------------------------------------------------------------------------
% INPUTS:
%   spike_ind -- indices of all detected spikes
%   G3 -- brainstorm structure with forward operator
%   Nmin -- minimal number of sources in one cluster
%   ValMax -- subcorr values from RAP-MUSIC
%   IndMax -- locations of sources from RAP-MUSIC
%   ind_m -- indices of sources survived after the subcorr thresholding
%   draw -- draw plot or not 
%   cortex -- cortical structure from brainstorm
%   RAP -- 'RAP' to use complete RAP-MUSIC procedure, smth else for one-round 
%   spikeind -- timeindices from RAP-MUSIC procedure
% 
% OUTPUTS:
%   cluster -- structure [length(ind_m)x4], first column -- source
%           location, second column -- spike timestamp, third -- the
%           subcorr value, fourth -- index of the spike from the whole set
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com


locs = cortex.Vertices;
dist = zeros(length(src_idx),length(src_idx));
for i = 1:length(src_idx) % distances between each vertex
    for j = 1:length(src_idx)
        dist(i,j) = norm(locs(src_idx(i),:)-locs(src_idx(j),:));
    end
end

cluster = zeros(1,length(src_idx));
fl = 1; 
k = 1;
idc = 1:length(src_idx);

while fl == 1
    dst = sum(dist < thr_dist, 2); 
    [val, ind] = max(dst); % vertex with the highest number of close neighbours
    if val >= Nmin
        ind_nhb = find(dist(ind,:) < thr_dist); % neighbours
        cluster(idc(ind_nhb)) = k;
        %src_idx = src_idx(:,setdiff(1:size(dist,1), ind_nhb));
        idc = idc(1,setdiff(1:size(dist,1),ind_nhb));
        dist = dist(setdiff(1:size(dist,1), ind_nhb),setdiff(1:size(dist,1), ind_nhb));
        
        k = k + 1;
        fl = 1;
    else
        fl = 0;
    end
end

if draw == 1
    cortex_lr = cortex;
    cortex_hr = cortex;
    
 
     c =[flipud(lines(7)); 0, 0, 0; ...
         0, 0.5, 0.5; 0.5,1.0,0.5; 0.6,0.2,0.5; 0.5,0.5,0.5;0.9,0.8,0.1;...
         0.8, 0.08, 0.5; 0.9, 0.5, 0.5;lines(7); 0.15, 0.15, 0.15; ...
         0, 0.5, 0.5; 0.3, 0.18, 0.4; ...for src_hemi in fwd_fixed['src']]
         0.8, 0.08, 0.5; 0.9, 0.5, 0.5; 0 0 0; 0 0 0 ; 0 0 0; 0 0 0 ; 0 0 0; 0 0 0; 0 0 0]; 
    data_lr = ones(length(cortex_lr.Vertices),1);
    mask_lr = zeros(size(data_lr));
    figure
    plot_brain_cmap2(cortex_lr, cortex_lr, [], data_lr, ...
        mask_lr, 0.3)
    hold on

    for i = 1:k-1
%         if length(c(:,1))==i
%             print('too many clusters');
%             break;
%         end
        ind = src_idx(i == cluster);
        scatter3(cortex.Vertices(ind,1), cortex.Vertices(ind,2), ...
            cortex.Vertices(ind,3), 100, 'filled', 'MarkerEdgeColor','k',...
                'MarkerFaceColor',c(i,:));
     
    end
    %view(270, 90)
%     subplot(1,2,2)
%     plot_brain_cmap2(cortex_lr, cortex_lr, [], data_lr, ...
%         mask_lr, 0.3)
%     hold on
%     for i = 1:length(cluster)
% %         if length(c(:,1))==i
% %             print('too many clusters');
% %             break;
% %         end
%         ind = cluster{1,i}(1,:);
%         scatter3(cortex.Vertices(ind,1), cortex.Vertices(ind,2), ...
%             cortex.Vertices(ind,3), 100, 'filled', 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor',c(i,:));
%     end
     view(180, 0)

end
end
