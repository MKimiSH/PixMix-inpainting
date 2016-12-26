function tmouse_smallst(action)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if(nargin == 0)
    action = 'start';
end
global this_frame;
global landmarks;
global fres;
global v;
switch(action)
  case 'start'
   %axis([0 size(this_frame,1) 0 size(this_frame,2)]);% 设定图轴范围
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
   se = strel('disk',4);
   fres = maxLianTongYu_smallst(imclose(detect_obj_smallst(this_frame, landmarks),se));   
   
   %fres = fres(1:4:end,1:4:end);
   
   this_boundary = edge(fres, 'sobel');
   
   
   last_frame = this_frame;
   %last_frame = im2uint8(rgb2gray(this_frame));
   %this_frame = last_frame;
   
   %this_frame = this_frame(1:4:end,1:4:end);
   
   i = 0;
   [H,this_boundary,opticalFlow,~,~,~] = object_tracking(last_frame,this_frame,true,this_boundary,[]);
   figure,imshow(this_boundary); 
   while hasFrame(v)
      i = i + 1
      last_frame = this_frame;
      
      this_frame = readFrame(v);
      %this_frame = im2uint8(rgb2gray(readFrame(v)));
      %this_frame = this_frame(1:4:end,1:4:end);
      
      [H,this_boundary,opticalFlow,this_corner_list,estimated_corner_list,flow] = object_tracking(last_frame,this_frame,false,this_boundary,opticalFlow);
      this_boundary_list = matrix2list(this_boundary,1);
      
      %ͼƬ���
      if i==50
          display = insertMarker(this_frame, fliplr(this_corner_list), '+','color','yellow');
          display = insertMarker(display, fliplr(estimated_corner_list), 'circle','color','red');
          display = insertMarker(display, fliplr(this_boundary_list), '+','color','green');
          imshow(display);
          hold on
          plot(flow,'DecimationFactor',[8 8],'ScaleFactor',800)
          figure
      end
   end
   
end
end

