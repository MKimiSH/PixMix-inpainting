

windpath = 'F:\TDDOWNLOAD\uTorrentDown\起风了.2013.台粤日英四语.简繁中字￡CMCT小鱼\[起风了].The.Wind.Rises.2013.BluRay.720p.x264.FLAC.4Audios-CMCT.mkv';
llpath = '"F:\TDDOWNLOAD\uTorrentDown\情书(国粤日) Love.Letter.1995.BluRay.720p.x264.3Audio.AC3-CnSCG.mkv"';

% wind rises 38:06.5 - 38:09
% cutOneVideo(windpath, 'wind1.mp4', [0, 0, 6.5], [0, 0, 2.5]);
% cmd = 'C:\M.K.S.H\Utilities\ffmpeg\bin\ffmpeg.exe -i "F:\TDDOWNLOAD\uTorrentDown\情书(国粤日) Love.Letter.1995.BluRay.720p.x264.3Audio.AC3-CnSCG.mkv" -ss 00:10:49 -t 3.7 -acodec aac -vcodec h264 -strict -2 "C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\ll1.mp4"';
% cmd = [cmd, '-ss 00:00:6.50 -t 00:00:2.50 -acodec aac -vcodec h264 -strict -2 "C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\wind1.mp4"';]
cmd = 'C:\M.K.S.H\Utilities\ffmpeg\bin\ffmpeg.exe -i "F:\TDDOWNLOAD\uTorrentDown\起风了.2013.台粤日英四语.简繁中字￡CMCT小鱼\[起风了].The.Wind.Rises.2013.BluRay.720p.x264.FLAC.4Audios-CMCT.mkv" -ss 00:38:6.50 -t 00:00:2.8 -acodec aac -vcodec h264 -strict -2 "C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\wind1.mp4"';
cmd = [cmd, ' -ss 00:38:10.5 -t 3 -acodec aac -vcodec h264 -strict -2 "C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\wind2.mp4"'];
cmd = [cmd, ' -ss 00:57:42.3 -t 6.7 -acodec aac -vcodec h264 -strict -2 "C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\wind3.mp4"'];
cmd = [cmd, ' -ss 01:05:38.4 -t 9.5 -acodec aac -vcodec h264 -strict -2 "C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\wind4.mp4"'];
cmd = [cmd, ' -ss 01:09:48.4 -t 4.5 -acodec aac -vcodec h264 -strict -2 "C:\M.K.S.H\MATLAB Workspace at SSD\inpainting with PixMix\videos\cut\wind5.mp4"'];

system(cmd);