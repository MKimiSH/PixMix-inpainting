function tmouse(action)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if(nargin == 0)
    action = 'start';
end
global snooker;
global landmarks;
global fres;
switch(action)
  case 'start'
   %axis([0 size(snooker,1) 0 size(snooker,2)]);% 设定图轴范围
  % box on;% 将图轴加上图框
   %title('Click and drag your mouse in this window!');
   set(gcf, 'WindowButtonDownFcn', 'tmouse down');
  case 'down'
   set(gcf, 'WindowButtonMotionFcn', 'tmouse move'); 
   set(gcf, 'WindowButtonUpFcn', 'tmouse up');
  case 'move'
   currPt = get(gca, 'CurrentPoint');
   x = currPt(1,1);
   y = currPt(1,2);
   %line(x, y, 'marker', '.', 'EraseMode', 'xor');
   plot(x,y,'w.');
   landmarks = [landmarks; floor(x),floor(y)];
   %landmarks(floor(x), floor(y),:) = true;
   %fprintf('Mouse is moving! Current location = (%g, %g)\n', currPt(1,1), currPt(1,2));
  case 'up'
   set(gcf, 'WindowButtonMotionFcn', '');
   set(gcf, 'WindowButtonUpFcn', '');
   landmarks = unique(landmarks, 'rows', 'stable');
   fres = (detect_obj(snooker, landmarks));   
end
end

