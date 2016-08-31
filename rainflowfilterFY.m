clear;
clc;
addpath(genpath(pwd));
load('FY_RAVG.mat');
format long
signal.data=double(signal.data);
fy= transpose(signal.data);
tp=dat2tp(fy);
Time= double((1:signal.samples)/ signal.samplerate);

% [y] = rfcfilter(x,h,def);
%Input:
%    x   = Signal.   [nx1] OR [nx2]
 %   h   = Threshold for rainflow filter.
%   def = 0: removes cycles with range < h. (default amplitude)
%             1: removes cycles with range <= h.

amp=500;
y = rfcfilter(tp,amp,0);


figure(1),
subplot(211), 
plot(tp);
xlabel('number of recorded points','fontsize',20);
ylabel('loads(N)','fontsize',20);
title('Fy in the Left Front Weel with points recorded','fontsize',20);

subplot(212),
plot(y);
xlabel('filtered points','fontsize',20);
ylabel('loads(N)','fontsize',20);
title('Fy in the Left Front Weel with filtered points where amplitude<500N ','fontsize',20);