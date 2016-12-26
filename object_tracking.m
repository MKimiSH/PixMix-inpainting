%使用前需要先addpath(genpath('libs/tracking'));
%输入:uint8灰度图像last_frame,uint8灰度图像this_frame,逻辑值this_frame是否是第一帧is_first_frame,上一帧的边界last_boundary,需要维持的全局变量opticalFlow
%输出:投影变换H,当前帧边界this_contour,需要维持的全局变量opticalFlow
function [H,this_boundary,opticalFlow,last_corner_list,flow] = object_tracking(last_frame,this_frame,is_first_frame,last_boundary,opticalFlow)
    %% 第一帧进行光流初始化
    
    if is_first_frame == true 
        opticalFlow = opticalFlowLK('NoiseThreshold',0.009);
        last_corner_list = [];
        flow = estimateFlow(opticalFlow,this_frame);
        H = [1,0,0;0,1,0;0,0,1];
        this_boundary = last_boundary;
        return
    end
    
    %% 
    
    [img_height,img_width] = size(this_frame);

    %求出上一帧到这一帧的光流
    opticalFlow = opticalFlowLK('NoiseThreshold',0.009);
    flow = estimateFlow(opticalFlow,last_frame);
    flow = estimateFlow(opticalFlow,this_frame);
    
    %flow = estimateFlow(opticalFlow,this_frame);
    max(max(flow.Vx))
    
    %figure,plot(flow,'DecimationFactor',[5 5],'ScaleFactor',10)
    
    last_boundary_list = matrix2list(last_boundary,1);
    [last_boundary_count,~] = size(last_boundary_list);
    
    %膨胀一次，增加candidate数量
    se = ones(4,4);
    candidate_matrix = imdilate(last_boundary,se);
    candidate_list = matrix2list(candidate_matrix,1);
    
    %估计投影变换需要至少4对corner,因此要求harris返回上一帧至少4个corner
    %last_corner_list = harris(last_frame,candidate_list,4); %只搜索candidate中是否有交点
    last_corner_list = harris(last_frame); %所有全局的角点
    
    matchedpoints_last = last_corner_list;%找到的上一帧的角点位置
    
    matchedpoints_this = matchedpoints_last;%找到的角点根据光流法对应出的这一帧的位置
    [count,~] = size(matchedpoints_this);
    for i = 1:count %计算上一帧的角点按照光流法,这一帧应该到哪了
        last_y = matchedpoints_last(i,1);
        last_x = matchedpoints_last(i,2);
       
        dx = round(flow.Vx(last_y,last_x));
        dy = round(flow.Vy(last_y,last_x));
        
        this_y = last_y+dy;
        this_x = last_x+dx;
        
        this_y = max(min(this_y,img_height),1);
        this_x = max(min(this_x,img_width),1);
        
        matchedpoints_this(i,:) = [this_y,this_x];
    end

    [tform, ~, ~, status] = estimateGeometricTransform(fliplr(matchedpoints_last),fliplr(matchedpoints_this),'projective');
    if status == 2 %如果估计投影变换时inlier太少，认为两帧之间没有运动
        warning('inlier not enough');
        H = [1,0,0;0,1,0;0,0,1];
        this_boundary = last_boundary;
        
        return
    end
    
    H = tform.T%从角点位置变换,得到两帧之间的投影变换关系

    %% 根据投影变换H，上一帧的boundary位置，估计这一帧的boundary位置
    
    this_boundary_list = fliplr(last_boundary_list);
    for i = 1:last_boundary_count
        cur_point = this_boundary_list(i,:);
        extended_cur_point = [cur_point';1];
        transformed_extended_cur_point = H'*extended_cur_point;
        transformed_x = round(transformed_extended_cur_point(1) / transformed_extended_cur_point(3));
        transformed_y = round(transformed_extended_cur_point(2) / transformed_extended_cur_point(3));
        this_boundary_list(i,:) = [transformed_y,transformed_x];
    end

    this_boundary = list2matrix(this_boundary_list,1,img_height,img_width);
    this_boundary = imclose(this_boundary,strel('disk',3)); %修整轮廓

end
