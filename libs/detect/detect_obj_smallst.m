%输入:彩色图像img,坐标序列landmarks
%输出:0,1边界res
function res = detect_obj_smallst(img, landmarks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sz = size(img);
image = double(img);
res = zeros(sz(2:-1:1));
gamma = 0.95;
l = size(landmarks,1);
U = zeros(l, 3);
m3 = meanfilter_smallst(image, 3);
%m5 = meanfilter(img,5);
m7 = meanfilter_smallst(image,7);

%U = mf(m3, m7, landmarks);
for i = 1:l
    U(i,:)= mf_smallst(m3, m7, landmarks(i,2),landmarks(i,1));
end
%fenlei
C = clustering_smallst(U);
size(C)
% vi
lC = size(C,2);
%vi = zeros(3);
for i = 1:lC
    v(i,:) = std(C{i},[], 1);
end
%for i = 1:size(v,2)
%    vi(i) = max(v(:,i));
%end
vi = max(v,[],1)
minx = min(landmarks(:,2));
miny = min(landmarks(:,1));
maxx = max(landmarks(:,2));
maxy = max(landmarks(:,1));
snookercp = img;
% [x,y] = meshgrid(minx:maxx,miny:maxy);
% x = reshape(x, 1, numel(x));
% y = reshape(y, 1, numel(y));
% in = inpolygon(x,y, landmarks(:,2), landmarks(:,1));
res = (mex_detect_smallst(landmarks(:,2), landmarks(:,1),m3,m7, U,vi,[minx, maxx, miny, maxy],gamma));
a = size(res,1);
b = size(res,2);
for i = 1:a
    for j = 1:b
        if(res(i,j) == 1)
           snookercp(i,j,:) = [255,255,25];
        end
    end
end

% 
% for i = minx:maxx
%     for j = miny:maxy
%        if(inpolygon(i,j,landmarks(:,2), landmarks(:,1)))
%             sumU = 0;
%             for k = 1:l
%                 sumU = sumU + dis(mf(m3,m7,i,j), U(k,:), vi);
%             end
%             %fprintf('%d\n',sumU);
%             if(sumU >= l * gamma)
%                 snookercp(i,j,:) = [255,0,0];
%                 res(j,i) = 1;
%             end
%        else
%             snookercp(i,j,:)=[255,255,255];
%        end
%        
%     end
% end

% for i = minx: maxx
%     for j = miny: maxy
%        % if(pointInPolygon([j,i],landmarks))
%        if(inpolygon(landmarks(:,1), landmarks(:,2), [j],[i]))
%             cccc = cccc + 1;
%             sumU = 0;
%             for k = 1:l
%                 sumU = sumU + dis(mf(m7,m3,i,j), U(k,:), vi);
%             end
%             sumU;
%             if(sumU >= l * gamma)
%                 snookercp(i,j,:) = [0,0,0];
%                 res(j,i) = 1;
%             end
%         else
%            snookercp(i,j,:)=[255,255,255];
%         end
%     end
% end
% cccc/((maxx-minx)*(maxy-miny))
for i = 1 :size(landmarks,1)
    snookercp(landmarks(i,2),landmarks(i,1),:)=  [255,0,0];
end
%figure;
imshow(snookercp);
end

