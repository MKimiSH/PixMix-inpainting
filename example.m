clear;
close all

addpath(genpath('libs'));

global snooker;
global v;
global landmarks;
global fres;

%snooker = imread('s.jpg');
%imshow(snooker);
%hold on;
v = VideoReader('fr1.mp4');
snooker = readFrame(v);
imshow(snooker);
hold on;
%landmarks = false(size(snooker));
landmarks = [];
tmouse_smallst();