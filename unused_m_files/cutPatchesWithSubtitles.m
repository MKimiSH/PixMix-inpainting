function [] = cutPatchesWithSubtitles(dirname)
dbstop if error 

%% preparation
% dirname = 'G:\caffe\data\movie-0\batch-0\';
cutdirname = [dirname, 'patch_subtitle\'];

if exist(cutdirname, 'dir')
    return;
end

mkdir(cutdirname);

o_ext = '.jpg';
s_ext = '-sub.jpg';
w_ext = '-white.jpg';
b_ext = '-black.jpg';
m_ext = '-mask.jpg';

%% calculate the size to preserve
D = dir(dirname);
imNum = floor((length(D)-3)/4);
cutTL = [10000, 10000];
cutBR = [0,0];
bb = imread([dirname, num2str(1), b_ext]);
bigSz = [size(bb,1), size(bb,2)];

parfor i=1:imNum
    white = imread([dirname, num2str(i), w_ext]);
    white = 255-white;
    black = imread([dirname, num2str(i), b_ext]);
    % both black and white are mostly zeroes!!
    if max(white(:)) == 0 || max(black(:)) == 0
        continue;
    end
    if ndims(white) == 3
        white = rgb2gray(white);
        black = rgb2gray(black);
    end
    [rows, cols] = find(white>0);
    tlw = [min(rows), min(cols)];
    brw = [max(rows), max(cols)];
    
    [rows, cols] = find(black>0);
    tlb = [min(rows), min(cols)];
    brb = [max(rows), max(cols)];
    
    cutTL = min(min(tlw, tlb), cutTL);
    cutBR = max(max(brw, brb), cutBR);
    
end
cutCtr = floor((cutTL+cutBR)/2);
if cutBR(1) ~= 0 || cutBR(2) ~= 0
%     cutTL = max(cutTL - [150, 150], [1,1]);
%     cutBR = min(cutBR + [150, 150], bigSz);
    cutBR = min(cutCtr + [128,128], bigSz);
    %cutTL = max(cutCtr - [64,64], [1,1]);
    cutTL = cutBR - [255,255];
else
    %cutTL = round([bigSz(1)/2, bigSz(2)/3]);
    cutBR = round([bigSz(1), bigSz(2)/3*2]);
    cutTL = cutBR - [255,255];
end
%% cut every image
parfor i = 1:imNum
    orig = imread([dirname, num2str(i), o_ext]);
    imwrite(orig(cutTL(1):cutBR(1), cutTL(2):cutBR(2), :), [cutdirname, num2str(i), o_ext]);
    
    sub = imread([dirname, num2str(i), s_ext]);    
    imwrite(sub(cutTL(1):cutBR(1), cutTL(2):cutBR(2), :), [cutdirname, num2str(i), s_ext]);

    white = imread([dirname, num2str(i), w_ext]);
%     imwrite(white(cutTL(1):cutBR(1), cutTL(2):cutBR(2), :), [cutdirname, num2str(i), w_ext]);
    
    white = 255-white;
    black = imread([dirname, num2str(i), b_ext]);
    % both black and white are mostly zeroes!!
    if max(white(:)) == 0 || max(black(:)) == 0
        continue;
    end
    if ndims(white) == 3
        white = rgb2gray(white);
        black = rgb2gray(black);
    end  
    mask = (white>0) | (black>0);
%     imwrite(black(cutTL(1):cutBR(1), cutTL(2):cutBR(2), :), [cutdirname, num2str(i), b_ext]);
    imwrite(mask(cutTL(1):cutBR(1), cutTL(2):cutBR(2), :), [cutdirname, num2str(i), m_ext]);
end