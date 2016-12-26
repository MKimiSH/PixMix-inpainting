%% yuzhang
%输入:uint8灰度图像img,n行2列矩阵candidate_list存储可能是角点的点的坐标[y(height),x(width)],要求检测到的最少的角点数min_corner_num
%输入img,candidate_list(可选),min_corner_num(可选)
%输出:标定的角点坐标列表n行2列的corner_list,存储格式[y(height),x(width)]

%% 用例
% img = imread('building.jpg');
% img = im2uint8(rgb2gray(img));
% [height,width] = size(img);
% corner_list = harris(img);
% corner_matrix = list2matrix(corner_list,1,height,width);
% imshow(corner_matrix);

%Harris角点检测对尺度敏感
function [corner_list] = harris(varargin)
    %% 处理变长输入
    input_var_length = length(varargin);
    if input_var_length == 0
        warning('error');
    else %input_var_length >=1
        img = varargin{1};
        candidate_list = []; %candidate_list = []时认为没有candidate约束,全局检查
        min_corner_num = -1; %min_corner_num = -1时认为没有min_corner_num约束,有多少个理想角点,返回多少个角点
        if input_var_length >= 2
            candidate_list = varargin{2};
            if input_var_length >=3
                min_corner_num = varargin{3};
            end
        end
    end
    
    %% 差分+平滑
    
    %可以不用一阶差分，而用Prewitt
    %dx = [-1 0 1;-1 0 1;-1 0 1];  %dx：横向Prewitt差分模版  
    dx = [-1,1];
    Ix = filter2(dx,img);
    Iy = filter2(dx',img);
    Ix2 = Ix.^2;     
    Iy2 = Iy.^2;  
    Ixy = Ix.*Iy;

    h = fspecial('gaussian',[7 7],2);%模板尺寸[7,7],sigma=2
    Ix2 = filter2(h,Ix2);
    Iy2 = filter2(h,Iy2);
    Ixy = filter2(h,Ixy);

    %% 根据差分结果计算各点响应
    
    [img_height, img_width] = size(img);
    Rmax = 0;
    k = 0.06; %常取0.04-0.06
    R = zeros(img_height,img_width); %计算图像中每个点的角点响应
    for i = 1:img_height
        for j = 1:img_width
            M = [Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)]; %偏导数矩阵
            R(i,j) = det(M)-k*(trace(M))^2; %角点响应函数
            if R(i,j) > Rmax
                Rmax = R(i,j);
            end
        end
    end

    if isempty(candidate_list) == 0 %如果指定了candidate_list，可以缩小计算范围
        %% 计算每个点比邻域中多少个邻点响应更大

        tmp = zeros(img_height+2,img_width+2);
        tmp(2:img_height+1,2:img_width+1) = R;

        [count,~] = size(candidate_list);
        extended_candidate_list = zeros(count,4);
        extended_candidate_list(:,1:2) = candidate_list;

        for i = 1:count
            y = extended_candidate_list(i,1);
            x = extended_candidate_list(i,2);
            extended_candidate_list(i,3) = R(y,x);

            satisfied_condition = 0;
            directions = [ -1,-1 ; -1,0 ; -1,1 ; 0,-1 ; 0,1 ; 1,-1 ; 1,0 ; 1,1]; %8邻域
            for j = 1:length(directions) %计算每个节点比directions存储的邻域中多少个邻点更大
                dx = directions(j,1);
                dy = directions(j,2);
                biasx = 1;
                biasy = 1;
                if tmp(y+biasy,x+biasx) > tmp(y+biasy+dy,x+biasx+dx)
                    satisfied_condition = satisfied_condition + 1;
                end
            end
            extended_candidate_list(i,4) =  satisfied_condition;
        end

        %% 升序排序,优先要求角点响应强度尽量是局部极值,其次要求角点响应强度越大越好

        extended_candidate_list = sortrows(extended_candidate_list,[4 3]);
        extended_candidate_list = flipud(extended_candidate_list);%转成降序排序

        %根据排序结果返回一组角点，尽可能满足min_corner_num的要求
        if count < min_corner_num %如果要求的corner数量比所有可能的candidate还多,就把所有的candidate都作为角点
            corner_list = extended_candidate_list(1:count,1:2);
        elseif (extended_candidate_list(min_corner_num,4) == 8) && (extended_candidate_list(min_corner_num,3) > 0.01*Rmax) %如果确保了至少有min_corner_num个优质角点
            sum = 0;
            for i = 1:count
                if (extended_candidate_list(i,4) == 8) && (extended_candidate_list(i,3) > 0.01*Rmax)
                    sum = sum+1;
                end
            end
            corner_list = extended_candidate_list(1:sum,1:2);%最后返回所有的优质角点
        else %如果candidate很多，但优质角点不多，则返回min_corner_num个相对最好的candidate
            corner_list = extended_candidate_list(1:min_corner_num,1:2);
        end
    else %没有指定candidate_list时要对整幅图像计算
        tmp = zeros(img_height+2,img_width+2);
        tmp(2:img_height+1,2:img_width+1) = R;
        img_re = zeros(img_height+2,img_width+2);
        for i = 2:img_height+1
            for j = 2:img_width+1
                if tmp(i,j)>0.01*Rmax &&...
                   tmp(i,j)>tmp(i-1,j-1) && tmp(i,j)>tmp(i-1,j) && tmp(i,j)>tmp(i-1,j+1) &&...
                   tmp(i,j)>tmp(i,j-1) && tmp(i,j)>tmp(i,j+1) &&...
                   tmp(i,j)>tmp(i+1,j-1) && tmp(i,j)>tmp(i+1,j) && tmp(i,j)>tmp(i+1,j+1)
                        img_re(i,j)=1; %3*3极大响应,且比历史最大响应的0.01更大时，认为是角点
                end
            end
        end

        corner_matrix=zeros(img_height,img_width);
        corner_matrix(1:img_height,1:img_width)=img_re(2:img_height+1,2:img_width+1);
        
        corner_list = matrix2list(corner_matrix,1);
    end
end
