data=xlsread('Cdata.xlsx');
data=data(:,3);
Nstd = 0.2;
NR = 500;
MaxIter =2000;


[modes its]=ceemdan(data,0.2,500,2000);
t=1:length(data);

[a b]=size(modes);

figure;
subplot(a+1,1,1);
plot(t,data);
ylabel('data')
set(gca,'xtick',[])
axis tight;
title('CEEMDAN')

for i=2:a
    subplot(a+1,1,i);
    plot(t,modes(i-1,:));
    ylabel (['IMF ' num2str(i-1)]);
    set(gca,'xtick',[])
    xlim([1 length(data)])
end;

subplot(a+1,1,a+1)
plot(t,modes(a,:))
ylabel(['IMF ' num2str(a)])
xlim([1 length(data)])

