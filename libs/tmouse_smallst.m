function tmouse_smallst(action)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if(nargin == 0)
    action = 'start';
end
global snooker;
global landmarks;
global fres;
global v;
switch(action)
  case 'start'
   %axis([0 size(snooker,1) 0 size(snooker,2)]);% 设定图轴范围
  % box on;% 将图轴加上图�?
   %title('Click and drag your mouse in this window!');
   set(gcf, 'WindowButtonDownFcn', 'tmouse_smallst down');
  case 'down'
   set(gcf, 'WindowButtonMotionFcn', 'tmouse_smallst move'); 
   set(gcf, 'WindowButtonUpFcn', 'tmouse_smallst up');
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
   fres = (detect_obj_smallst(snooker, landmarks));   
   last_contour = edge(fres, 'sobel');
   
   %object_tracking(last_frame,this_frame,is_first_frame,last_contour,opticalFlow)
   last_frame = im2uint8(rgb2gray(snooker));
   
   snooker = im2uint8(rgb2gray(readFrame(v)));
    i = 0;
    [H,this_contour,opticalFlow] = object_tracking(snooker,last_frame,1,last_contour,[]);
       figure,imshow(this_contour); 
   while hasFrame(v)
      i = i + 1
      if i > 500
        break;
      end
       last_frame = snooker;
       snooker = im2uint8(rgb2gray(readFrame(v)));
       [H,this_contour,opticalFlow] = object_tracking(snooker,last_frame,0,last_contour,opticalFlow);
       imshow(this_contour);
   end
   
end
end

