clear;
close all

addpath(genpath('libs'));

global this_frame;
global v;
global landmarks;
global fres;

%v = VideoReader('fr1.mp4');
v = VideoReader('rock.mp4');
%v = VideoReader('snk2c2.mp4');

this_frame = readFrame(v);
imshow(this_frame);
hold on;
%landmarks = false(size(this_frame));
landmarks = [];
tmouse_smallst();