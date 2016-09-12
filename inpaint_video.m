% inpaint a video.
% the mask should be a single image.
dbstop if error

base = 'videos\';
vidName =  'small_out8.mp4';
maskName = 'mask_more_bin.png'; % mask is a single image.
fileName = [base, vidName];
vidName(end-3:end) = [];
startTime = 0;
endTime = 3;

vidObj = VideoReader(fileName);
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
vidObj.CurrentTime = startTime;

vidStr = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
kk = 1;
while vidObj.CurrentTime <= endTime
%     vidStr(kk).cdata = imresize(readFrame(vidObj), 0.5);
	vidStr(kk).cdata = readFrame(vidObj);
    kk = kk+1;
end
timeSpan = kk-1; % number of frames in vidStr

%% read mask M
maskName = [base, maskName];
M = imread(maskName);
if ~ismatrix(M)
    M = rgb2gray(M);
end
% I = im2single(I);
M = im2single(M);
M = M>0;

[R, C] = size(M);
%% process video frame by frame
outVid = zeros(R,C,3,timeSpan);
figure;
for fr = 1:timeSpan
    I = im2double(vidStr(fr).cdata);
    outVid(:,:,:,fr) = fillOneImage(I,M,0);
    imshow(squeeze(outVid(:,:,:,fr)), 'InitialMagnification', 'fit');
end

