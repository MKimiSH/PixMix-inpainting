clear;
close all

addpath(genpath('libs'));

global snooker;
global v;
global landmarks;
global fres;

v = VideoReader('fr1.mp4');
%v = VideoReader('rock.mp4');

snooker = readFrame(v);
imshow(snooker);
hold on;
%landmarks = false(size(snooker));
landmarks = [];
tmouse_smallst();