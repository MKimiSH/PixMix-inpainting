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
   this_boundary = edge(fres, 'sobel');
   
   size(fres)
   
   last_frame = im2uint8(rgb2gray(snooker));
   %snooker = last_frame;
   snooker = im2uint8(rgb2gray(readFrame(v)));
   snooker = snooker(1:2:end,1:2:end);
   
   i = 0;
   [H,this_boundary,opticalFlow,~,~] = object_tracking(last_frame,snooker,true,this_boundary,[]);
   figure,imshow(this_boundary); 
   while hasFrame(v)
      i = i + 1
      if i > 10
        break;
      end
      last_frame = snooker;
      
      snooker = im2uint8(rgb2gray(readFrame(v)));
      snooker = snooker(1:2:end,1:2:end);
      
      last_boundary = this_boundary;
      
      [H,this_boundary,opticalFlow,last_corner_list,flow] = object_tracking(last_frame,snooker,false,this_boundary,opticalFlow);
      
      change = isequal(last_boundary, this_boundary)
      
      this_corner_list = harris(snooker);
      this_boundary_list = matrix2list(this_boundary,1);
      
      display = insertMarker(snooker, fliplr(last_corner_list), '+','color','blue');
      display = insertMarker(display, fliplr(this_corner_list), '+','color','red');
      display = insertMarker(display, fliplr(this_boundary_list), '+','color','green');
      imshow(display);
      hold on
      plot(flow,'DecimationFactor',[1 1],'ScaleFactor',10)
      
      figure
   end
   
end
end

