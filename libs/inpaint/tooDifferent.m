function [isdiff] = tooDifferent(im1, im2)
im1 = im2uint8(im1);
im2 = im2uint8(im2);

h1 = histImage(im1);
h2 = histImage(im2);

isdiff = norm(h1-h2) > 0.005

end

function [h] = histImage(im)
h = zeros(255,1);
tot = nnz(im);
for i=1:255
   h(i) = nnz(im==i);
end
h = h./tot;

end