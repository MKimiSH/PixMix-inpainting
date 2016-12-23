clear;
% 
% landmarks = [83,48; 
%     75,48;
%     70,50;
%     65,55;
%     65,60;
%     65,68;
%     70,70;
%     75,72;
%     80,72;
%     85,68;
%     88,62;
%     86,55];
%test1 = insertMarker(snooker, landmarks, 'circle','color','red');imshow(test1)

global snooker;
global v;
global landmarks;
global fres;

%snooker = imread('s.jpg');
%imshow(snooker);
%hold on;
v = VideoReader('output.mp4');
snooker = readFrame(v);
imshow(snooker);
hold on;
%landmarks = false(size(snooker));
landmarks = [];
tmouse_smallst();