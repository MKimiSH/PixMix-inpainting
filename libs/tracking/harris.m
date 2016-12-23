%yuzhang
%输入:uint8灰度图像img,n行2列矩阵candidate_list存储可能是角点的点的坐标,要求检测到的最少的角点数min_corner_num
%输出:标定的角点坐标列表n行2列的corner_list

%用例:
%img=imread('building.jpg');
%img = im2uint8(rgb2gray(img)); 
%imshow(img);
%corner_matrix = harris(img);
%figure,imshow(mat2gray(corner_matrix));
%Harris角点检测对尺度敏感
function [corner_list] = harris(img,candidate_list,min_corner_num)
    [img_height, img_width] = size(img);

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

    Rmax = 0;
    k = 0.06;%常取0.04-0.06
    R = zeros(img_height,img_width);%计算图像中每个点的角点响应
    for i = 1:img_height
        for j = 1:img_width
            M = [Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)];%偏导数矩阵
            R(i,j) = det(M)-k*(trace(M))^2;%角点响应函数
            if R(i,j) > Rmax
                Rmax = R(i,j);
            end
        end
    end

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
        if tmp(y+1,x+1)>tmp(y,x)
            satisfied_condition = satisfied_condition +1;
        end
        if tmp(y+1,x+1)>tmp(y+1,x)
            satisfied_condition = satisfied_condition +1;
        end
        if tmp(y+1,x+1)>tmp(y+2,x)
            satisfied_condition = satisfied_condition +1;
        end
        if tmp(y+1,x+1)>tmp(y,x+1)
            satisfied_condition = satisfied_condition +1;
        end
        if tmp(y+1,x+1)>tmp(y+2,x+1)
            satisfied_condition = satisfied_condition +1;
        end
        if tmp(y+1,x+1)>tmp(y,x+2)
            satisfied_condition = satisfied_condition +1;
        end
        if tmp(y+1,x+1)>tmp(y+1,x+2)
            satisfied_condition = satisfied_condition +1;
        end
        if tmp(y+1,x+1)>tmp(y+2,x+2)
            satisfied_condition = satisfied_condition +1;
        end
        extended_candidate_list(i,4) =  satisfied_condition;
    end
    
    %升序排序,优先要求角点响应强度尽量是局部极值,其次要求角点响应强度越大越好
    extended_candidate_list = sortrows(extended_candidate_list,[4 3]);
    extended_candidate_list = flipud(extended_candidate_list);%转成降序排序
    
    if (extended_candidate_list(min_corner_num,4) == 8) && (extended_candidate_list(min_corner_num,3) > 0.01*Rmax)
        sum = 0;
        for i = 1:count
            if (extended_candidate_list(i,4) == 8) && (extended_candidate_list(i,3) > 0.01*Rmax)
                sum = sum+1;
            end
        end
        corner_list = extended_candidate_list(1:sum,1:2);
    else
        corner_list = extended_candidate_list(1:min_corner_num,1:2);
    end
end
