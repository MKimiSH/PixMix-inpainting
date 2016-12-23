%yuzhang
%输入:list(n行2列),复原值target_value,复原的矩阵的尺寸[height,width];
%输出:matrix(height行width列,在list标记的位置取值为target_value,其余位置取值为0)
function [matrix] = list2matrix(list,target_value,height,width)
    matrix = zeros(height,width);
    [count,~] = size(list);
    for i = 1:count
        y = list(i,1);
        x = list(i,2);
        if (y <= height) && (y >= 1) && (x <= width) && (x >= 1)
            matrix(y,x) = target_value;
        end
    end
end