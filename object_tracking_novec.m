%输入:rgb图像last_frame,rgb图像this_frame,逻辑值this_frame是否是第一帧is_first_frame,上一帧的边界last_boundary,需要维持的全局变量opticalFlow
%输出:投影变换H,当前帧边界this_contour,需要维持的全局变量opticalFlow
function [H,this_boundary, this_mask, opticalFlow,this_corner_list,estimated_corner_list,flow] = object_tracking_novec(last_frame,this_frame,is_first_frame,last_boundary,opticalFlow)

    uint8_this_frame = im2uint8(rgb2gray(this_frame));
    
    %% 第一帧进行光流初始化
    tic;
    if is_first_frame == true 
        %opticalFlow = opticalFlowLKDoG('NumFrames',3);
        opticalFlow = opticalFlowFarneback;
        %opticalFlow = opticalFlowHS; %opticalFlowLK在很多地方值都为0
        this_corner_list = [];
        estimated_corner_list = [];
        this_mask = [];
        flow = estimateFlow(opticalFlow,uint8_this_frame);
        H = [1,0,0;0,1,0;0,0,1];
        this_boundary = last_boundary;
        return
    end
    
    %% 
    
    uint8_last_frame = im2uint8(rgb2gray(last_frame));
    [img_height,img_width] = size(uint8_this_frame);

    %求出上一帧到这一帧的光流
    flow = estimateFlow(opticalFlow,uint8_this_frame);
    
    last_boundary_list = matrix2list(last_boundary,1);
    [last_boundary_count,~] = size(last_boundary_list);
    
    %膨胀一次，增加candidate数量
    se = ones(8,8);
    candidate_matrix = imdilate(last_boundary,se);
    candidate_list = matrix2list(candidate_matrix,1);
    
    %估计投影变换需要至少4对corner,因此要求harris返回上一帧至少4个corner
    this_corner_list = harris(uint8_this_frame,candidate_list,30); %只搜索candidate中是否有交点
    %this_corner_list = harris(uint8_this_frame); %所有全局的角点
    
    estimated_corner_list = this_corner_list;%找到的角点根据光流法对应出的这一帧的位置
    [count,~] = size(estimated_corner_list);
    
    t = toc; fprintf('%.4f sec for optical flow\n', t); 
    
    for i = 1:count %计算上一帧的角点按照光流法,这一帧应该到哪了
        this_y = this_corner_list(i,1);
        this_x = this_corner_list(i,2);
       
        vx = flow.Vx(this_y,this_x);
        vy = flow.Vy(this_y,this_x);
        
        min_dist = inf;
        min_dist_x = 0;
        min_dist_y = 0;
        neighborhood = 10;
        for dx = -neighborhood:neighborhood
            for dy = -neighborhood:neighborhood
                if (dx==0) && (dy==0)
                    continue
                end
                
                tmp_y = this_y + dy;
                tmp_x = this_x + dx;
                
                if (tmp_y <= 0) || (tmp_y > img_height) || (tmp_x <= 0) || (tmp_x > img_width)
                    continue
                end
                
                tmp_vx = flow.Vx(tmp_y,tmp_x);
                tmp_vy = flow.Vy(tmp_y,tmp_x);
                
                dist = (vx-tmp_vx)^2 + (vy-tmp_vy)^2;
                if dist < min_dist
                    min_dist = dist;
                    min_dist_y = tmp_y;
                    min_dist_x = tmp_x;
                end
            end
        end
        
        estimated_corner_list(i,:) = [min_dist_y,min_dist_x];
    end
    
    %'similarity'效果不好
    %[tform, ~, ~, status] = estimateGeometricTransform(fliplr(this_corner_list),fliplr(estimated_corner_list),'affine');
    [tform, ~, ~, status] = estimateGeometricTransform(fliplr(this_corner_list),fliplr(estimated_corner_list),'projective');
    if status == 2 %如果估计投影变换时inlier太少，认为两帧之间没有运动
        warning('inlier not enough');
        H = [1,0,0;0,1,0;0,0,1];
        this_boundary = last_boundary;
        return
    end
    
    %% 使用vision.PointTracker找两帧之间的变换关系
    
    %初始化vision.PointTracker
    tic;
    imagePoints1 = detectMinEigenFeatures(uint8_last_frame, 'MinQuality', 0.1);  
    tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);  
    imagePoints1 = imagePoints1.Location;  
    initialize(tracker, imagePoints1, uint8_last_frame);
    t = toc; fprintf('%.4f sec for init Point Tracker\n', t); 
   
    %vision.PointTracker
    [imagePoints2, validIdx] = step(tracker, uint8_this_frame);  
    matchedPoints1 = imagePoints1(validIdx, :);  
    matchedPoints2 = imagePoints2(validIdx, :);  
    
    [tform, ~, ~, ~] = estimateGeometricTransform(matchedPoints1,matchedPoints2,'projective');
    

    %% 根据投影变换H，上一帧的boundary位置，估计这一帧的boundary位置
    
    H = tform.T';%从角点位置变换,得到两帧之间的投影变换关系
    
    this_boundary_list = fliplr(last_boundary_list);
    for i = 1:last_boundary_count
        cur_point = this_boundary_list(i,:);
        extended_cur_point = [cur_point';1];
        transformed_extended_cur_point = H*extended_cur_point;
        transformed_x = round(transformed_extended_cur_point(1) / transformed_extended_cur_point(3));
        transformed_y = round(transformed_extended_cur_point(2) / transformed_extended_cur_point(3));
        this_boundary_list(i,:) = [transformed_y,transformed_x];
    end

    this_boundary = list2matrix(this_boundary_list,1,img_height,img_width);
    
    
    %%
    
    %this_boundary = imclose(this_boundary,strel('disk',3)); %修整轮廓
    
    se = ones(8,8);
    this_boundary = imdilate(this_boundary,se);
    this_boundary_list = matrix2list(this_boundary,1);
    
    this_boundary_list_y = this_boundary_list(:,1);
    this_boundary_list_x = this_boundary_list(:,2);
    
    K = convhull(this_boundary_list_x,this_boundary_list_y);
    this_boundary_list = [this_boundary_list_y(K),this_boundary_list_x(K)];
    
    
    
    x_list = [];
    y_list = [];
    
    [count,~] = size(this_boundary_list);
    for i=1:count
        
        cur_y = this_boundary_list(i,1);
        cur_x = this_boundary_list(i,2);
        if i ~= count
            next_y = this_boundary_list(i+1,1);
            next_x = this_boundary_list(i+1,2);
        else
            next_y = this_boundary_list(1,1);
            next_x = this_boundary_list(1,2);
        end
        
        n = max(abs(cur_y-next_y),abs(cur_x-next_x));
        dy = (next_y - cur_y)/n;
        dx = (next_x - cur_x)/n;
            
        for j=0:n-1
            tmp_y = round(cur_y + j*dy);
            tmp_x = round(cur_x + j*dx);
                
            x_list(end+1) = tmp_x;
            y_list(end+1) = tmp_y;
        end
            
    end
    
    this_mask = maxLianTongYu_smallst(imclose(detect_obj_smallst(this_frame, [x_list',y_list']),se)); 
    this_mask = imdilate(this_mask, se);
    this_boundary = edge(this_mask, 'sobel');
    
%     [X,Y] = meshgrid(1:img_height,1:img_width);
%     
%     this_boundary = inpolygon(X,Y,this_boundary_list(:,2),this_boundary_list(:,1));
%     this_boundary = edge(this_boundary, 'sobel');
%     this_boundary_list = matrix2list(this_boundary,1);
%     
%     se = strel('disk',4);
%     this_boundary = maxLianTongYu_smallst(imclose(detect_obj_smallst(this_frame, fliplr(this_boundary_list)),se));   
%    
%     this_boundary = edge(this_boundary, 'sobel');

end
