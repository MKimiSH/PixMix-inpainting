function I = maskImage(I, M)

if (ismatrix(I))
    I(M) = 0;
    return;
end
i1 = I(:,:,1);
i2 = I(:,:,2);
i3 = I(:,:,3);

i1(M)=0; i2(M)=0; i3(M)=0;

I(:,:,1) = i1;
I(:,:,2) = i2;
I(:,:,3) = i3;

end
