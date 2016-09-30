function [filledI, retf] = fillOneLevel(initf, I, M, level, useLineConstr)
% initf is initial value of the mapping f
% Use random sampling and propagation to fill the missing parts of I with
% the known pixels.

% use 8-neighbor for cost_spatial, 
% use 5x5 patch for cost_apperance.

global f;
global R;
global C;
global im;
global mask;
f = initf;
numIter = 5 ;
[R, C] = size(M);
im = I;
mask = M;
D = bwdist(~M);
alphaSp = 0.05;
alphaAp = 0.95;
wrs = 0.5; 
rrs = min(R,C); % weight random search and range random search
%% fill the image with initial guess of f
for i=1:R
    for j=1:C
        if M(i,j) == 1
            im(i,j,:) = im(f(i,j,1), f(i,j,2), :);
        end
    end
end

%% the propagation and random search
for iter = 1:numIter
    for i=1:R
        for j=1:C
            if M(i,j) == 1
                im(i,j,:) = im(f(i,j,1), f(i,j,2), :);
            end
        end
    end
%     imshow(im);
    if mod(iter,2) == 1
        for i=1:R
            for j=1:C
                %% propagation
                if M(i,j)==1
                    fijOld = squeeze(f(i,j,:));
                    costLeft = Inf;
                    costUp = Inf;
                    spcOld = calcSpatialCost(i,j, fijOld);
                    apcOld = calcAppearanceCost(i,j, fijOld);
                    costOld = spcOld*alphaSp + apcOld*alphaAp;
                    if M(i-1, j)==1 && i-1>0 && isValid(squeeze(f(i-1,j,:))+[1;0])
                        fijLeft = squeeze(f(i-1,j,:))+[1;0];
                        spcLeft = calcSpatialCost(i,j, fijLeft);
                        apcLeft = calcAppearanceCost(i,j, fijLeft);
                        costLeft = spcLeft*alphaSp + apcLeft*alphaAp;
                    end
                    if M(i, j-1)==1 && j-1>0 && isValid(squeeze(f(i,j-1,:))+[0;1])
                        fijUp = squeeze(f(i,j-1,:))+[0;1];
                        spcUp = calcSpatialCost(i,j, fijUp);
                        apcUp = calcAppearanceCost(i,j, fijUp);
                        costUp = spcUp*alphaSp + apcUp*alphaAp;
                    end
                    
                    if costLeft < min(costUp, costOld)
                        f(i,j,:) = fijLeft;
                        costOld = costLeft;
                    elseif costUp < min(costLeft, costOld)
                        f(i,j,:) = fijUp;
                        costOld = costUp;
                    end
                    
                    %% random search
                    drs = wrs*rrs;
                    while(drs>D(i,j)+3)
                        searchvec = round((rand(1,2)-0.5)*drs) + [i,j];
                        num = 3;
                        while(~isValid(searchvec) && num>0)
                            searchvec = round((rand(1,2)-0.5)*drs) + [i,j];
                            num = num-1;
                        end
                        fijRS = searchvec';
                        spcRS = calcSpatialCost(i,j, fijRS);
                        apcRS = calcAppearanceCost(i,j, fijRS);
                        costRS = spcRS*alphaSp + apcRS*alphaAp;
                        
                        if costRS<costOld
                            f(i,j,:) = fijRS;
                            costOld = costRS;
                        end
                        
                        drs = drs*wrs;
                    end
                end
                
            end % end for j
        end % end for i
        
    else % iter is even
        for i=R:-1:1
            for j=C:-1:1
                
                if M(i,j)==1
                    fijOld = squeeze(f(i,j,:));
                    costRight = Inf;
                    costDown = Inf;
                    spcOld = calcSpatialCost(i,j, fijOld);
                    apcOld = calcAppearanceCost(i,j, fijOld);
                    costOld = spcOld*alphaSp + apcOld*alphaAp;
                    if M(i+1, j)==1 && i+1<=R && isValid(squeeze(f(i+1,j,:))+[-1;0])
                        fijRight = squeeze(f(i+1,j,:))+[-1;0];
                        spcRight = calcSpatialCost(i,j, fijRight);
                        apcRight = calcAppearanceCost(i,j, fijRight);
                        costRight = spcRight*alphaSp + apcRight*alphaAp;
                    end
                    if M(i, j+1)==1 && j+1>0 && isValid(squeeze(f(i,j+1,:))+[0;-1])
                        fijDown = squeeze(f(i,j+1,:))+[0;-1];
                        spcDown = calcSpatialCost(i,j, fijDown);
                        apcDown = calcAppearanceCost(i,j, fijDown);
                        costDown = spcDown*alphaSp + apcDown*alphaAp;
                    end
                    
                    if costRight < min(costDown, costOld)
                        f(i,j,:) = fijRight;
                    elseif costDown < min(costRight, costOld)
                        f(i,j,:) = fijDown;
                    end
                    %% random search
                    drs = wrs*rrs;
                    while(drs>D(i,j))
                        searchvec = round((rand(1,2)-0.5)*drs) + [i,j];
                        num = 3;
                        while(~isValid(searchvec) && num>0)
                            searchvec = round((rand(1,2)-0.5)*drs) + [i,j];
                            num = num-1;
                        end
                        fijRS = searchvec';
                        spcRS = calcSpatialCost(i,j, fijRS);
                        apcRS = calcAppearanceCost(i,j, fijRS);
                        costRS = spcRS*alphaSp + apcRS*alphaAp;
                        
                        if costRS<costOld
                            f(i,j,:) = fijRS;
                            costOld = costRS;
                        end
                        
                        drs = drs*wrs;
                    end
                end
                
            end % end for j
        end % end for i
        
    end % end if
    
end
retf = f;
filledI = im;

end

% fij may change due to propagation or random search.
function spc = calcSpatialCost(i,j, fij)
global f;
global R;
global C;

spc = 0;
w = 1/8;
tao = 200;
for r = max(i-1,1):min(i+1,R)
    for c = max(j-1,1):min(j+1,C)
        % f(p) = fij, f(p+v) = f(r,c,:), v = [r-i], c-j], f(p)+v = fij+v
        v = [r,c] - [i,j];
        fv = squeeze(f(r, c, :));
        diff = fv - (fij+v');
        spc = spc + min(norm(diff)^2, tao);
    end
end
spc = spc*w;

end

% fij may change due to propagation or random search.
function apc = calcAppearanceCost(i,j, fij)
global f;
global R;
global C;
global im;
apc = 0;
w = 1/25;
for r = max(i-2,1):min(i+2,R)
    for c = max(j-2,1):min(j+2,C)
        % f(p) = fij, f(p+v) = f(r,c,:), v = [r-i, c-j], f(p)+v = fij+v
        v = [r,c] - [i,j];
        i1 = im(r,c,:);
        if isValid([fij(1)+v(1), fij(2)+v(2)])
            i2 = im(fij(1)+v(1), fij(2)+v(2), :);
        else 
            i2 = zeros(size(i1));
        end
        
        apc = apc + norm(squeeze(i1-i2))^2;
    end
end
apc = apc*w;

end

function v = isValid(vec)
global R;
global C;
global mask;

i = vec(1);
j = vec(2);

v = i>0&&i<=R && j>0&&j<=C && mask(i,j)==0;

end
