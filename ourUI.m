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

select = uicontrol('String', 'select', 'Callback', ...
    {@tmouse, 'select'}, ...
    'Position', [20 20 60 20]);

avoid = uicontrol('String', 'avoid area', ...
    'Callback', {@tmouse, 'avoid'}, ... 
    'Position', [100 20 100 20]);

lines = uicontrol('String', 'lines', ...
    'Callback', @drawLines, ...
    'Position', [220 20 60 20]);

start = uicontrol('String', 'start', 'Enable', 'off',...
    'Callback', {@mystart},...  
    'Position', [300 20 60 20]);

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
                plot(x,y,'w.');
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
        Landmarks
        Lines
        avoidMarks
    end
end    