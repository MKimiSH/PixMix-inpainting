function d = dis_smallst(Uk, Uj, v)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    d = 0;
    for i = 1:size(Uk,2)
        if(abs(double(Uk(i)) - double(Uj(i))) > v(i))
            d = 1;
            break;
        %else
           % fprintf('%d-%d,%d\n', Uk(i),Uj(i), v(i));
        end
        
    end
end

