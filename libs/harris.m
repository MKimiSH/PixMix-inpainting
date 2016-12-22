%用例:
%img=imread('building.jpg');
%img = im2uint8(rgb2gray(img)); 
%imshow(img);
%corner_matrix = harris(img);
%figure,imshow(mat2gray(corner_matrix));

%输入:uint8灰度图像img
%输出:标定的角点，矩阵中为1的地方是角点，为0的地方非角点
%Harris角点检测对尺度敏感
function [corner_list] = harris(img,candidate_list,min_corner_num)
    [m, n]=size(img);

    %可以不用一阶差分，而用Prewitt
    %dx = [-1 0 1;-1 0 1;-1 0 1];  %dx：横向Prewitt差分模版  
    
    dx = [-1,1];
    Ix = filter2(dx,img);
    Iy = filter2(dx',img);
    Ix2 = Ix.^2;     
    Iy2 = Iy.^2;  
    Ixy = Ix.*Iy;

    h=fspecial('gaussian',[7 7],2);%模板尺寸[7,7],sigma=2
    Ix2=filter2(h,Ix2);
    Iy2=filter2(h,Iy2);
    Ixy=filter2(h,Ixy);

    Rmax=0;
    k=0.06;%常取0.04-0.06
    R=zeros(m,n);%计算图像中每个点的角点响应
    for i=1:m
        for j=1:n
            M=[Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)];%偏导数矩阵
            R(i,j)=det(M)-k*(trace(M))^2;%角点响应函数

            if R(i,j)>Rmax
                Rmax=R(i,j);
            end
        end
    end

    tmp=zeros(m+2,n+2);
    tmp(2:m+1,2:n+1)=R;
    
    [height,~] = size(candidate_list);
    extended_candidate_list = zeros(height,4);
    extended_candidate_list(:,1:2) = candidate_list;
    
    for i=1:height
        x = extended_candidate_list(i,1);
        y = extended_candidate_list(i,2);
        extended_candidate_list(i,3) = R(x,y);
        
        count_satisfied_condition = 0;
        if tmp(x+1,y+1)>tmp(x,y)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        if tmp(x+1,y+1)>tmp(x,y+1)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        if tmp(x+1,y+1)>tmp(x,y+2)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        if tmp(x+1,y+1)>tmp(x+1,y)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        if tmp(x+1,y+1)>tmp(x+1,y+2)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        if tmp(x+1,y+1)>tmp(x+2,y)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        if tmp(x+1,y+1)>tmp(x+2,y+1)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        if tmp(x+1,y+1)>tmp(x+2,y+2)
            count_satisfied_condition = count_satisfied_condition +1;
        end
        extended_candidate_list(i,4) =  count_satisfied_condition;
    end
    
    %升序排序,优先要求角点响应强度尽量是局部极值,其次要求角点响应强度越大越好
    extended_candidate_list = sortrows(extended_candidate_list,[4 3]);
    %转成降序排序
    extended_candidate_list = rot90(extended_candidate_list);
    extended_candidate_list = rot90(extended_candidate_list);
    
    
    if (extended_candidate_list(min_corner_num,4) == 8) && (extended_candidate_list(min_corner_num,3) > 0.01*Rmax)
        count = 0;
        for i = 1:height
            if (extended_candidate_list(i,4) == 8) && (extended_candidate_list(i,3) > 0.01*Rmax)
                count = count+1;
            end
        end
        corner_list = extended_candidate_list(1:count,1:2);
    else
        corner_list = extended_candidate_list(1:min_corner_num,1:2);
    end
end
