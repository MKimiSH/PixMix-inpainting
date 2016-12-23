function [] = cutOneVideo(vidnm, outnm, bt, lt) 
fpath = 'C:\M.K.S.H\Utilities\ffmpeg\bin\ffmpeg.exe';
opath = ['"C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\', outnm, '"'];

cmd = [fpath, ' -i ', vidnm, ' -ss ', sprintf('%02d:%02d:%.2f -t %02d:%02d:%.2f', bt(1), bt(2), bt(3), ...
lt(1), lt(2), lt(3)), ' -acodec aac -vcodec h264 -strict -2 ', opath];
system(cmd);

end