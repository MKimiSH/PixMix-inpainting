%yuzhang
%输入:matrix,查询目标target_value
%输出:list(n行2列,存储matrix中所有值为traget_value的点的坐标[y(height),x(width)])
function [list] = matrix2list(matrix,target_value)
%     [height,width] = size(matrix);
    [row, col] = find(matrix == target_value);
    list = [row, col];
end

% function [list] = matrix2list(matrix,target_value)
%     [height,width] = size(matrix);
%     count = sum(matrix(:)==target_value);
%     list = zeros(count,2);
%     index = 1;
%     for y = 1:height
%         for x = 1:width
%             if matrix(y,x) == target_value
%                 list(index,:) = [y,x];
%                 index = index +1;
%             end
%         end
%     end
% end