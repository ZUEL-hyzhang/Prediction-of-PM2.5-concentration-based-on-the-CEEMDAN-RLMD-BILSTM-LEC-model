clc;
clear;
close all;
fs = 2923; % sampling frequency
N =2923; % data amount
t = (1:N)/fs; % time vector
x=xlsread('1010A.csv');
x=x(:,4);
options.display = 1;
options.max_iter = 30;
options.max_pfs = 6;
[pf3, ams3, fms3, ort3] = RLMD(x,options);
figure;
plot(t,x)