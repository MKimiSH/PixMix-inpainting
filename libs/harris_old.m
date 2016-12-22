%用例:
%img=imread('building.jpg');
%img = im2uint8(rgb2gray(img)); 
%imshow(img);
%corner_matrix = harris_old(img);
%figure,imshow(mat2gray(corner_matrix));

%输入:uint8灰度图像img
%输出:标定的角点，矩阵中为1的地方是角点，为0的地方非角点
%Harris角点检测对尺度敏感
function [corner_matrix] = harris_old(img)

    [m, n]=size(img);

    tmp=zeros(m+2,n+2);
    tmp(2:m+1,2:n+1)=img;
    Ix=zeros(m+2,n+2);
    Iy=zeros(m+2,n+2);

    Ix(:,2:n)=tmp(:,3:n+1)-tmp(:,1:n-1);%x方向差分
    Iy(2:m,:)=tmp(3:m+1,:)-tmp(1:m-1,:);%y方向差分

    %可以不用一阶差分，而用Prewitt
    %dx = [-1 0 1;-1 0 1;-1 0 1];  %dx：横向Prewitt差分模版  
    %Ix2 = filter2(dx,Image).^2;     
    %Iy2 = filter2(dx',Image).^2;  

    Ix2=Ix(2:m+1,2:n+1).^2;
    Iy2=Iy(2:m+1,2:n+1).^2;
    Ixy=Ix(2:m+1,2:n+1).*Iy(2:m+1,2:n+1);

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

    tmp(2:m+1,2:n+1)=R;
    img_re=zeros(m+2,n+2);
    for i=2:m+1
        for j=2:n+1
            if tmp(i,j)>0.01*Rmax &&...
               tmp(i,j)>tmp(i-1,j-1) && tmp(i,j)>tmp(i-1,j) && tmp(i,j)>tmp(i-1,j+1) &&...
               tmp(i,j)>tmp(i,j-1) && tmp(i,j)>tmp(i,j+1) &&...
               tmp(i,j)>tmp(i+1,j-1) && tmp(i,j)>tmp(i+1,j) && tmp(i,j)>tmp(i+1,j+1)
                    img_re(i,j)=1; %3*3极大响应,且比历史最大响应的0.01更大时，认为是角点
            end
        end
    end

    corner_matrix=zeros(m,n);
    corner_matrix(1:m,1:n)=img_re(2:m+1,2:n+1);
end
