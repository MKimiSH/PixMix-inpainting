%yuzhang
%输入:matrix,查询目标target_value
%输出:list(n行2列,存储matrix中所有值为traget_value的点的坐标)
function [list] = matrix2list(matrix,target_value)
    [height,length] = size(matrix);
    count = 0;
    for i = 1:height
        for j = 1:length
            if matrix(i,j) == target_value
                count = count + 1;
            end
        end
    end
    
    list = zeros(count,2);
    index = 1;
    for i = 1:height
        for j = 1:length
            if matrix(i,j) == target_value
                list(index,:) = [i,j];
                index = index +1;
            end
        end
    end
end