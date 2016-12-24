function [] = inpaintVideo(vidname, outname)
% main function of the project, reads a video by _vidname_, opens a GUI for selection,
% remove the selected object, output the resulting video with _outname_

%% read video

startTime = 0;
endTime = 3;

vidObj = VideoReader(vidname);
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
vidObj.CurrentTime = startTime;

vidStr = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'), 'colormap',[]);
kk = 1;
while vidObj.CurrentTime <= endTime
%     vidStr(kk).cdata = imresize(readFrame(vidObj), 0.5);
	vidStr(kk).cdata = readFrame(vidObj);
    kk = kk+1;
end
timeSpan = kk-1; % number of frames in vidStr

vid = cell(timeSpan, 1);
for i=1:timeSpan
    vid{i} = im2single(vidStr(i).cdata);
end

%% inpainting!!!
% call smallst's code to select object and inpaint first frame
[vid{1}, f0, contour] = inpaintFirstFrame(vid{1});

% Call zhangyu's code to track object, and then forward _f_ to get result
% for second frame
f = f0;
oflow = [];
H = [];
for i=2:timeSpan
    [vid{i}, f, contour, H, oflow] = inpaintSecondFrame(vid{i}, vid{i-1}, f, oflow);
end

%% write video
v = VideoWriter(outname, 'MPEG-4');
open(v);
v.FrameRate = vidObj.FrameRate;
for i=1:timeSpan
    writeVideo(v, vid{i});
end
close(v);

end