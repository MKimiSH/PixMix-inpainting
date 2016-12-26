function ourUI(video)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Landmarks = [];
avoidMarks = [];
Lines = [];
v = VideoReader(video);
startFrame = readFrame(v);
imshow(startFrame);
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
        fres = maxLianTongYu_smallst(imclose(detect_obj_smallst(this_frame, landmarks),se));   
        this_boundary = edge(fres, 'sobel'); % first frame boundary
        [h,w] = size(fres);
        [x,y] = meshgrid(1:w, 1:h);
        avoidArea = inpolygon(x,y, avoidMarks(:,1), avoidMarks(:,2)); % avoid area
        Lines % line point (n x 4 matrix)
    end
end    