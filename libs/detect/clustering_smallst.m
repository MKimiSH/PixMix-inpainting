function res = clustering_smallst(array)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    temp = array;
    E = 70;
    i = 0;
    while(isempty(temp) == 0)
        i = i + 1;
        l = size(temp,1);
        ri = randi(l);
        res{i} = temp(ri,:);
        m = temp(ri, :);
        temp(ri,:) = [];
        leni = 1;
        j = 1;
        while(j <= size(temp, 1))
             if(norm(m- temp(j,:)) < E)
                res{i} = [res{i};temp(j,:)];
                m = (m * leni + temp(j,:))/(leni+1);
                leni = leni + 1;
                temp(j,:) = [];
                j = j - 1;
            end
            j = j + 1;
        end
    end
end

