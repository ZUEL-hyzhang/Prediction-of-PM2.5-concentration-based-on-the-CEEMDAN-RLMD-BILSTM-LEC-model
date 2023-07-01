rms = inf;
rms1 = inf;
i1 = 0;
for i = 20:20:1000
    data1 = data;
    for j = 1:length(data)-1
        if abs(data(j)-data(j+1))>=i
            data1(j+1) = data(j+1)+err(j+1);
        end
    end
    rmse1 = (sum((data1-or).*(data1-or))/length(or))^1/2;
    if rms1<rms
        rms = rms1;
        i1 = i;
    end
end