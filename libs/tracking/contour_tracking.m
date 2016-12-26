function [contour_tracking_list] = contour_tracking(binary_img)
    [m n]=size(binary_img);
    
    imgn=zeros(m,n);        %边界标记图像
    ed=[-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0]; %从左上角像素判断
    for i=2:m-1
        for j=2:n-1
            if img(i,j)==1      %如果当前像素是前景像素

                for k=1:8
                    ii=i+ed(k,1);
                    jj=j+ed(k,2);
                    if img(ii,jj)==0    %当前像素周围如果是背景，边界标记图像相应像素标记
                        imgn(ii,jj)=1;
                    end
                end

            end
        end
    end
end
    
figure;
imshow(imgn,[]);
