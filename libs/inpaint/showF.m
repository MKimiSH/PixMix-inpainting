function showF(I, M, F)

figure, imshow(I), hold on

[R,C] = size(M);
for r = 1:R
    for c = 1:C
        if M(r,c) == 1 && random('unif', 0,1)>0.96
            f = squeeze(F(r,c,:));
            rf = [f(2), f(1)];
            plot([c,f(2)], [r,f(1)], 'LineWidth',1,'Color','cyan');
            plot(c,r,'x','LineWidth',2,'Color','yellow');
            plot(f(2),f(1),'x','LineWidth',2,'Color','red');
        end
    end
end
end