function [] = inpaintVideo(vidname, outname, M0, UM0, usrlines, C0)
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
[vid{1}, m0, f0, c0, of0, ofobj] = inpaintFirstFrame(vid{1}, M0, UM0, usrlines, C0);

% Call zhangyu's code to track object, and then forward _f_ to get result
% for second frame
f = f0;
of = of0;
c = c0;
m = m0;
H = [];
for i=2:timeSpan
    [vid{i}, m, f, c, H, of, ofobj] = inpaintSecondFrame(vid{i}, vid{i-1}, f, m, c, of, ofobj);
end

%% write video
if nargin==1
    outname = 'ret.mp4';
end
v = VideoWriter(outname, 'MPEG-4');
open(v);
v.FrameRate = vidObj.FrameRate;
for i=1:timeSpan
    writeVideo(v, vid{i});
end
close(v);

end