function ourUIImage(image)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Landmarks = [];
avoidMarks = [];
Lines = [];
startFrame = imread(image);
figure, imshow(startFrame);
hold on;

select = uicontrol('String', 'select', ...
    'Callback', {@tmouse, 'select'}, ...
    'Position', [20 620 60 20]);

avoid = uicontrol('String', 'avoid', ...
    'Callback', {@tmouse, 'avoid'}, ... 
    'Position', [20 580 60 20]);

lines = uicontrol('String', 'lines', ...
    'Callback', @drawLines, ...
    'Position', [20 540 60 20]);

start = uicontrol('String', 'start', 'Enable', 'off',...
    'Callback', {@mystart},...  
    'Position', [20 500 60 20]);

    function tmouse(hObject, eventdata,action)
        if(strcmp(action, 'avoid') || strcmp(action, 'select'))
            hd = struct;
            hd.tAction = action;
            hd.marks = [];   
            guidata(gcbo, hd);
            action = 'start';
        end
                
        hd = guidata(gcbo);
        
        switch(action)
            case 'start'
                set(gcf, 'WindowButtonDownFcn', {@tmouse, 'down'});
            case 'down'
                set(gcf, 'WindowButtonMotionFcn', {@tmouse, 'move'});
                set(gcf, 'WindowButtonUpFcn', {@tmouse, 'up'});
            case 'move'
                currPt = get(gca, 'CurrentPoint');
                x = currPt(1,1);
                y = currPt(1,2);
                if(strcmp(hd.tAction, 'avoid'))
                    plot(x,y,'r.');
                else
                    plot(x,y,'w.');                   
                end

                hd.marks = [hd.marks;floor(x), floor(y)];
                guidata(gcbo, hd);
            case 'up'
                set(gcf, 'WindowButtonMotionFcn', '');
                set(gcf, 'WindowButtonUpFcn', '');
                set(gcf, 'WindowButtonDownFcn', '');
                
                marks = unique(hd.marks, 'rows', 'stable');
                if(strcmp(hd.tAction, 'avoid'))
                    avoidMarks = marks;
                else
                    Landmarks = marks;
                end

                set(start,'Enable','on');
        end
    end
    function drawLines(hObject, eventdata)
        l = imline;
        pos = wait(l)
        Lines = [Lines;reshape(pos, 1, 4)];
    end
    function mystart(hObject, eventdata)
        se = strel('disk',4);
        fres = imclose(detect_obj_smallst(startFrame, Landmarks),se);  
        fres = imdilate(fres, strel('disk', 3));
        this_boundary = edge(fres, 'sobel'); % first frame boundary
        [h,w] = size(fres);
        [x,y] = meshgrid(1:w, 1:h);
        avoidArea = zeros(size(this_boundary), 'logical');
        if ~isempty(avoidMarks)
            avoidArea = inpolygon(x,y, avoidMarks(:,1), avoidMarks(:,2)); % avoid area
        end
%         Lines % line point (n x 4 matrix)
        Lines = round(Lines);
        usrlines = struct('point1', [-1 -1], 'point2', [0 0], 'theta', 0, 'rho', 0,...
            'k', 0, 'isgood', 1);
        for i=1:size(Lines, 1);
            usrlines(i).point1 = Lines(i, [1 3]);
            usrlines(i).point2 = Lines(i, [2 4]);
            [usrlines(i).theta,  usrlines(i).rho] = calcThetaRho(Lines(i,:));
        end
%         inpaintVideo(video, 'out.mp4', fres, avoidArea, usrlines, this_boundary);
        inpaintFirstFrame(startFrame, fres, avoidArea, usrlines, this_boundary);
    end
end

function [t, r] = calcThetaRho(p)
x1 = p(1)-1;
y1 = p(3)-1;
x2 = p(2)-1;
y2 = p(4)-1;

eq1 = sprintf('r = %d*cos(t) + %d*sin(t)', x1, y1);
eq2 = sprintf('r = %d*cos(t) + %d*sin(t)', x2, y2);
[r, t] = solve(eq1, eq2, '-pi/2<=t<pi/2');
r = eval(r); t = 180*eval(t)/pi;

end