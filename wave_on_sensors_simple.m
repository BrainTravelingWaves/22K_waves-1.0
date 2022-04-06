function [sensor_waves,directions,PATH, PATH_sm] = wave_on_sensors_simple(cortex,  PARAMS, G, vertices_sm)

    strt =PARAMS.seed_vertices;% 160578;%171473;%PARAMS.seed_vertices;
    vertices = cortex.Vertices;
    faces = cortex.Faces;
    speed = PARAMS.wave.speeds;
    duration = PARAMS.wave.duration;
    nspoints = PARAMS.nspoints;
    
    VertConn = cortex.VertConn;
    
    indn1 = find(VertConn(strt,:)); %indices for the first neighbours
    max_step = 120;  %?????????????????????? WAS 9
    num_dir = length(indn1); %number of directions
    indn = indn1;

    
    for n = 1:num_dir
        % find the direction vector
        IND(n,1) = strt; %first vertex
        ind0 = indn(n); %second vertex (for given direction)
        IND(n,2) = ind0; %second vert into vert array
        norm0 = mean(cortex.VertNormals(indn,:));
        norm0 = norm0/norm(norm0); %mean normal for the set of vertices
        Pnorm0 = eye(3)-norm0'*norm0; %projector out of the norm (on the plane)
        dr0 = (cortex.Vertices(ind0,:) - cortex.Vertices(strt,:));
        dr0 = dr0*Pnorm0';
        dr0 = dr0/norm(dr0);% define direction vector
        directions(:,n) = dr0;
        d = 3; %number of the step
        
        while(d <= max_step)
            indn1 = find(VertConn(ind0,:));%new neighbours
            clear cs;
            for n1 = 1:length(indn1)
                dr1 = (cortex.Vertices(indn1(n1),:) - cortex.Vertices(ind0,:)); %next minus prev
                dr1 = dr1*Pnorm0';
                dr1 = dr1/norm(dr1);

                cs(n1) = dr1*dr0';
            end
            if sum(cs < 0) ~= length(cs) || (d < 3)
                [csmaxval,csmaxind] = max(cs);%neighbour with max incrementio 
                ind0 = indn1(csmaxind);
                IND(n,d) = ind0;
            else
                dr0  =  (cortex.Vertices(IND(n,d-1),:) - cortex.Vertices(IND(n,d-2),:));
                dr0= dr0*Pnorm0';
                dr0 = dr0/norm(dr0);
                for n1 = 1:length(indn1)
                  cs(n1) = dr1*dr0';
                end
                [csmaxval,csmaxind] = max(cs);%neighbour with max incrementio 
                ind0 = indn1(csmaxind);
                IND(n,d) = ind0;
            end
           
            %if d>=3
            %    if IND(n,d) == IND(n,d-2)
            %        cs(csmaxind) = -1;
            %        [csmaxval,csmaxind] = max(cs);%neighbour with max incrementio 
            %        ind0 = indn1(csmaxind);
            %        IND(n,d) = ind0;
            %    end
            %end
            
            
            d = d+1;
        end
    end
        
    DIST = zeros(num_dir, (max_step-1));
    for d = 1:num_dir
        for i = 2:max_step
            DIST(d,(i-1)) = norm(vertices(IND(d,i),:)-vertices(IND(d,(i-1)),:));
        end
    end
        
    %G = from_bst_get_gain_matrix(PARAMS.forward.name, PARAMS); 
    
    SR = PARAMS.sampling_rate;
    ntpoints = round(SR*duration);%number of points in time for wave
    PATH = cell(num_dir, nspoints, length(speed));
    PATH_sm = cell(num_dir, nspoints, length(speed));
    FM = cell(num_dir, nspoints, length(speed));
    tstep = (1/SR);
    
    alpha = 0;
    for s = 1:length(speed)
        l = speed(s)*tstep;
        for d = 1:num_dir
            PATH{d,1,s} = vertices(strt,:); %path - direct x num of step x speed
            PATH_sm{d,1,s} = vertices_sm(strt,:);
            FM{d,1,s} = G(:,strt); %corresponding forward coef
            res = 0;
            v1 = 1;%prev vertex
            v2 = 2;%next vertex
            for t = 2:nspoints
                if (strt == 218301) && (s == 5) && (d == 4) && (t ==5)
                    bbb = 1;
                end
                if (alpha > 1) || (alpha < 0)
                    aaa = 1;
                end
                if l < res
                   alpha = 1-l/res;
                   PATH{d,t,s} = alpha*PATH{d,(t-1),s}+(1-alpha)*vertices(IND(d,v2),:);
                   PATH_sm{d,t,s} = alpha*PATH_sm{d,(t-1),s}+(1-alpha)*vertices_sm(IND(d,v2),:);
                   FM{d,t,s} = alpha*FM{d,(t-1),s}+ (1-alpha)*G(:,IND(d,v2));
                   res = res-l;
                elseif l > res
                    if res == 0
                        if l < DIST(d,(v2-1))
                            alpha = 1-l/DIST(d,(v2-1));
                            PATH{d,t,s} = alpha*vertices(IND(d,v1),:)+(1-alpha)*vertices(IND(d,v2),:);
                            PATH_sm{d,t,s} = alpha*vertices_sm(IND(d,v1),:)+(1-alpha)*vertices_sm(IND(d,v2),:);
                            FM{d,t,s} = alpha*G(:,IND(d,v1)) + (1-alpha)*G(:,IND(d,v2));
                            res = DIST(d,(v2-1))-l;
                        elseif l == DIST(d,(v2-1))
                            PATH{d,t,s} = vertices(IND(d, v2),:);
                            PATH_sm{d,t,s} = vertices_sm(IND(d, v2),:);
                            FM{d,t,s} = G(:, IND(d, v2));
                            v1 = v1 + 1;
                            v2 = v2 + 1;
                        else
                            l2 = l-DIST(d,(v2-1));
                            v1 = v1 + 1;
                            v2 = v2 + 1;
                            res = DIST(d,(v2-1))-l2;
                            while res < 0
                                l2 = -res;
                                v1 = v1 + 1;
                                v2 = v2 + 1;
                                res = DIST(d,(v2-1))-l2;
                            end
                            alpha = 1-l2/DIST(d,(v2-1));
                            PATH{d,t,s} = alpha*vertices(IND(d,v1),:)+(1-alpha)*vertices(IND(d,v2),:);
                            PATH_sm{d,t,s} = alpha*vertices_sm(IND(d,v1),:)+(1-alpha)*vertices_sm(IND(d,v2),:);
                            FM{d,t,s} = alpha*G(:, IND(d,v1)) + (1-alpha)*G(:,IND(d,v2));
                            
                        end
                    else
         
                        l2 = l-res;
                        v1 = v1 + 1;
                        v2 = v2 + 1;
                        res = DIST(d,(v2-1))-l2;
                        
                        while res < 0
                            l2 = -res;
                            v1 = v1 + 1;
                            v2 = v2 + 1;
                            res = DIST(d,(v2-1))-l2;
                        end
                        alpha = 1-l2/DIST(d,(v2-1));
                        PATH{d,t,s} = alpha*vertices(IND(d,v1),:)+(1-alpha)*vertices(IND(d,v2),:);
                        PATH_sm{d,t,s} = alpha*vertices_sm(IND(d,v1),:)+(1-alpha)*vertices_sm(IND(d,v2),:);
                        FM{d,t,s} = alpha*G(:,IND(d,v1))+ (1-alpha)*G(:,IND(d,v2));
                              
                    end
                else %l == res
                    PATH{d,t,s} = vertices(IND(d, v2),:);
                    PATH_sm{d,t,s} = vertices_sm(IND(d, v2),:);
                    FM{d,t,s} = G(:,IND(d, v2));
                    v1 = v1 + 1;
                    v2 = v2 + 1;
                end
            end
        end
    end
            
    
%         
% figure
% 
% %c_hr = repmat(0.3, length(vertices), 1)';
% %c_hr = repmat([0.5 0.5 1],length(vertices),1)';
% 
% %%
% 
% c_hr = [255 249 230]/255;%[0.7 0.9 1];
% trisurf(faces, vertices(:,1),vertices(:,2),vertices(:,3),'FaceColor',c_hr);
%  lighting phong;
%   material dull;
%   view(-12,15)
%   camva(1)
% hold on
% 
% %vert 171473
% a =find(sqrt(sum((vertices -[-0.0143093 -0.0610668 0.0739706]).^2,2))<0.001);
% 
% plot3(vertices(strt,1)*1.001,vertices(strt,2)*1.005,vertices(strt,3)*1.001,'*','Color','red','MarkerSize',20,'LineWidth',2)
% 
% nhb = vertices(VertConn(strt,:),:);
% s = 4;
% for i = 1:length(nhb)
%    hold on 
%    
% plot3(nhb(i,1), nhb(i,2), nhb(i,3),'o','Color','red','MarkerSize',5,'LineWidth',2)
% 
% end
% % 
% %  for j = 1:num_dir
% %          for i = 1:ntpoints
% %              plot3(PATH{j,i,s}(1), PATH{j,i,s}(2), PATH{j,i,s}(3),'ko','MarkerFaceColor',[255/255, 140/255,0])
% %          end%plot path
% %      end
% % 
% % 
% % 
% %  for j = 1:num_dir
% %       a = reshape(cell2mat(PATH(j,1:20,s)),3,20)';
% %              plot3(a(:,1), a(:,2), a(:,3),'Color',[255/255, 140/255,0],'LineWidth',4)%         end%plot path
% %         
% %  end
% 
% axis equal
% axis off
% 
% 
% 
% 
% 
% 
% 
%   
%         
% figure
% 
% %c_hr = repmat(0.3, length(vertices), 1)';
% %c_hr = repmat([0.5 0.5 1],length(vertices),1)';
% 
% %%
% 
% c_hr =[0.9 1 0.98];%[255 225 220]/255;
% trisurf(faces, vertices(:,1),vertices(:,2),vertices(:,3),'FaceColor',c_hr);
%  lighting phong;
%   material dull;
%   view(-12,15)
%   camva(5)
% hold on
% 
% %vert 171473
% %a =find(sqrt(sum((vertices -[-0.01618 -0.038945 0.10206]).^2,2))<0.0001);
% 
% %plot3(vertices(strt,1)*1.001,vertices(strt,2)*1.005,vertices(strt,3)*1.001,'*','Color','red','MarkerSize',20,'LineWidth',2)
% 
% nhb = vertices(VertConn(strt,:),:);
% s = 9;
% for i = 1:length(nhb)
%    hold on 
%    
% %plot3(nhb(i,1), nhb(i,2), nhb(i,3),'o','Color','red','MarkerSize',5,'LineWidth',2)
% 
% end
% 
%  for j = 1:num_dir
%          for i = 1:ntpoints
%              plot3(PATH{j,i,s}(1), PATH{j,i,s}(2), PATH{j,i,s}(3),'ko','MarkerFaceColor','k')%[255/255, 140/255,0])
%          end%plot path
%      end
% 
% 
% 
% %  for j = 1:num_dir
% %       a = reshape(cell2mat(PATH(j,1:20,s)),3,20)';
% %              plot3(a(:,1), a(:,2), a(:,3),'Color',[255/255, 140/255,0],'LineWidth',4)%         end%plot path
% %         
%  for j = 1:num_dir
%       a = reshape(cell2mat(PATH(j,1:20,s)),3,20)';
%              plot3(a(:,1), a(:,2), a(:,3),'Color','k','LineWidth',4)%         end%plot path
%         
%  end
% 
% axis equal
% axis off
% 
% 
% 
% 
%         
% figure
% 
% %c_hr = repmat(0.3, length(vertices), 1)';
% %c_hr = repmat([0.5 0.5 1],length(vertices),1)';
% 
% %%
% 
% c_hr = [0.9 1 0.98];%[255 225 220]/255;
% trisurf(faces, vertices(:,1),vertices(:,2),vertices(:,3),'FaceColor',c_hr);
%  lighting phong;
%   material dull;
%   view(-12,15)
%   camva(5)
% hold on
% 
% %vert 171473
% %a =find(sqrt(sum((vertices -[-0.01618 -0.038945 0.10206]).^2,2))<0.0001);
% 
% %plot3(vertices(strt,1)*1.001,vertices(strt,2)*1.005,vertices(strt,3)*1.001,'*','Color','red','MarkerSize',20,'LineWidth',2)
% 
% nhb = vertices(VertConn(strt,:),:);
% s = 9;
% for i = 1:length(nhb)
%    hold on 
%    
% %plot3(nhb(i,1), nhb(i,2), nhb(i,3),'o','Color','red','MarkerSize',5,'LineWidth',2)
% 
% end
% 
% alphadic = [0 0 0.3 0 1 0 0.6 0];
%  for j = 1:num_dir
%          for i = 1:ntpoints
%              plot3(PATH{j,i,s}(1), PATH{j,i,s}(2), PATH{j,i,s}(3),'ko','Color',[0 0 0 alphadic(j)])%[255/255, 140/255,0])
%          end%plot path
%      end
% 
% 
% 
% %  for j = 1:num_dir
% %       a = reshape(cell2mat(PATH(j,1:20,s)),3,20)';
% %              plot3(a(:,1), a(:,2), a(:,3),'Color',[255/255, 140/255,0],'LineWidth',4)%         end%plot path
% %         
%  for j = 1:num_dir
%       a = reshape(cell2mat(PATH(j,1:20,s)),3,20)';
%              plot3(a(:,1), a(:,2), a(:,3),'Color',[0 0 0 alphadic(j)],'LineWidth',4)%         end%plot path
%         
%  end
% 
% view(-17,48)
% camva(1)
% axis equal
% axis off
% 
% 
% 
% 
% 
% 
% 
% a =1
% 
% %%
% 
% 
% 
% 
% %%%%%%%%%%%
%h = trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3));   
if PARAMS.draw_paths == 1
    close all
    
    figure
    h = trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3));
    set(h,'FaceAlpha',0.5);
    s = 9;%number of speed to plot
    hold on
    for i = 1:num_dir
        hplot = plot3(vertices(IND(i,:),1), vertices(IND(i,:),2), vertices(IND(i,:),3),'ko','MarkerFaceColor','r')
    end%plot vertices
    for j = 1:num_dir
        for i = 1:nspoints
            hplot = plot3(PATH{j,i,s}(1), PATH{j,i,s}(2), PATH{j,i,s}(3),'ko','MarkerFaceColor','g')
        end%plot path
    end
    hplot = plot3(vertices(strt,1),vertices(strt,2),vertices(strt,3),'ko','MarkerFaceColor','b')
end
 
%     
% %%%%%%%%%%%
% 

% 
    t = 0:(ntpoints-1);
    n = 0:nspoints-1;
    
    for i = t
       wave((i+1),:) = (1 + cos(2*pi * (n - i) / nspoints));
    end%wave - t (0 ... T) x 
    if PARAMS.draw_wave
        close all
        figure
        for i = t
            plot(wave(i+1,:))
            hold on
        end
    end

%     figure
%     for i = 1:6
%     plot(t, wave(i,:))
%     hold on
%     end
       
    % waves on sensors  
    for s = 1:length(speed)
        for i = 1:num_dir
            for t = 1:ntpoints
                for k = 1:nspoints
                    FM_s(:,k) = FM{i,k,s};
                end
                sensor_waves{i,s}(:,t) = FM_s*wave(t,:)';
            end
        end
    end
    %sensor waves - {dir x speed} (sensors x time)
end





