function [f] = vec_initMap(orif, M)
fprintf('start vecinimap\n');
tic
[R, C] = size(M);

if isempty(orif)
    ind1 = find(M==1);
    ind0 = find(M==0);
    l0 = length(ind0);
    l1 = length(ind1);
    f = zeros(R,C,2, 'int32');
    f(:,:,1) = meshgrid(1:R, 1:C)';
    f(:,:,2) = meshgrid(1:C, 1:R);
    [row1, col1] = ind2sub(size(M), ind1);
    [row0, col0] = ind2sub(size(M), ind0);
    randfill = randi(l0,l1,1);
    f_ind1 = row0(randfill);
    f_ind2 = col0(randfill);
    f(sub2ind(size(f), row1, col1, ones(l1,1))) = f_ind1;
    f(sub2ind(size(f), row1, col1, 2*ones(l1,1))) = f_ind2;
    return;
end
[R1, C1, ~] = size(orif);
tf = zeros(R1*2, C1*2, 2);
tf(1:2:R1*2, 1:2:C1*2, :) = orif*2 - 1; % tl
tf(1:2:R1*2, 2:2:C1*2, 1) = tf(1:2:R1*2, 1:2:C1*2, 1); % tr
tf(1:2:R1*2, 2:2:C1*2, 2) = tf(1:2:R1*2, 1:2:C1*2, 2) + 1; 
tf(2:2:R1*2, 1:2:C1*2, 1) = tf(1:2:R1*2, 1:2:C1*2, 1) + 1; % bl
tf(2:2:R1*2, 1:2:C1*2, 2) = tf(1:2:R1*2, 1:2:C1*2, 2);
tf(2:2:R1*2, 2:2:C1*2, 1) = tf(1:2:R1*2, 1:2:C1*2, 1) + 1; % br
tf(2:2:R1*2, 2:2:C1*2, 2) = tf(1:2:R1*2, 1:2:C1*2, 2) + 1;

f = tf(1:R, 1:C, :);
[row1, col1] = find(M==1 & (f(:,:,1)>=R | f(:,:,2)>=C));
[row0, col0] = find(M==0);
l1 = length(row1); 
l0 = length(row0);
randfill = randi(l0, l1, 1);
f_ind1 = row0(randfill);
f_ind2 = col0(randfill);
f(sub2ind(size(f), row1, col1, ones(l1,1))) = f_ind1;
f(sub2ind(size(f), row1, col1, 2*ones(l1,1))) = f_ind2;
toc
end