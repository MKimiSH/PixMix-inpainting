function f = initializeMap(orif, M)
fprintf('start iniMap\n');
tic
[R, C] = size(M);
f = zeros(R, C, 2, 'int32'); % 1 for row mapping, 2 for column
% f = int32(f);
% dbM = repmat(M, [1 2]);
% dbM = reshape(dbM, [size(M,1), size(M,2), 2]);

if isempty(orif)
%     f(:,:,1) = round(rand(R,C)*(R-1)) + 1;
%     f(:,:,2) = round(rand(R,C)*(C-1)) + 1;
    f(:,:,1) = randi(R, [R,C]);
    f(:,:,2) = randi(C, [R,C]);
    
    for i=1:R
        for j=1:C
            while M(i,j)==1 && M(f(i,j,1), f(i,j,2))==1
                f(i,j,1) = randi(R,1);
                f(i,j,2) = randi(C,1);
            end
            if M(i,j)==0
                f(i,j,:) = [i,j];
            end
        end
    end
    return;
end

[R1, C1, ~] = size(orif);
for i=1:R1
    for j=1:C1
        % shamelessly easy interpolation
        f(i*2-1, j*2-1, :) = orif(i,j,:)*2 - 1;
        f(i*2-1, j*2, :) = squeeze(f(i*2-1, j*2-1, :)) + int32([0,1]');
        f(i*2, j*2-1, :) = squeeze(f(i*2-1, j*2-1, :)) + int32([1,0]');
        f(i*2, j*2, :) = squeeze(f(i*2-1, j*2-1, :)) + int32([1,1]');
        
    end
end
f = f(1:R, 1:C, :);
% if there should be any violation of the mask, randomizing is the simplest
% solution.
cnt = 0;
for i=1:R
    for j=1:C
%         fprintf('p:%d, %d, f(p): %d, %d\n', i,j,f(i,j,1), f(i,j,2));
        while M(i,j)==1 && (f(i,j,1)>= R || f(i,j,2)>=C || M(f(i,j,1), f(i,j,2))==1)
            f(i,j,1) = randi(R,1);
            f(i,j,2) = randi(C,1);
            cnt = cnt+1;
        end
        if M(i,j)==0
            f(i,j,:) = [i,j];
        end
    end
end
toc
fprintf('%d violations and %d total mask size', cnt, nnz(M));

end