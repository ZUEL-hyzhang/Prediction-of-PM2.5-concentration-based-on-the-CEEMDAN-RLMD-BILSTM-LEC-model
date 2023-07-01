w1 = 0;
w2 = 1;
w = 0
ww = 1
rmse = inf;
while w1<=1
    w1 = w1+0.01;
    w2 = w2-0.01;
    s = w1*s1+w2*s2;
    rmse1 = (sum((s-or).*(s-or))/length(or))^1/2;
    if rmse1<rmse
        rmse = rmse1;
        w = w1;
        ww = w2;
    end
end
print(rmse,w,ww)

