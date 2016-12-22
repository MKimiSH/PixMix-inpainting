%yuzhang
%输入:list(n行2列),复原值target_value,复原的矩阵的尺寸[height,width];
%输出:matrix(height行width列,在list标记的位置取值为target_value,其余位置取值为0)
function [matrix] = list2matrix(list,target_value,height,width)
    matrix = zeros(height,width);
    [list_height,~] = size(list);
    for i = 1:list_height
        x = list(i,1);
        y = list(i,2);
        matrix(x,y) = target_value;
    end
end