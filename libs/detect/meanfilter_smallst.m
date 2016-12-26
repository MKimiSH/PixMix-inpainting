function RGB_Average  = meanfilter_smallst(img, num)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 %average = fspecial('average',num);
 average = fspecial('gaussian',1,0.5);
 RGB_Average = imfilter(img, average);
end

