clear all; 
close all;

w = 23;
amplitude = 4;
intercept = 4000;
len = 1;

w_line = 60;
line_noise = 100;

FS = 128;
lo = 4;
hi = 45;

[b, a] = butter(1,[lo hi]/(FS/2), 'bandpass');
fvtool(b,a);

x = 0:(1/FS):len;

y = intercept + amplitude * sin(2*pi*w*x) + line_noise * sin(2*pi*w_line*(x-rand));

y2 = filter(b,a,y);

plot(x,y);
hold on;
plot(x,y2);
legend('original','filtered');

figure;
idx = round(length(y)/4);
plot(x(idx:end),y(idx:end));
figure;
plot(x(idx:end),y2(idx:end));


%% notch 
% use on power lines of 60 Hz
wo = 60/(FS/2);  
bw = wo/35;

[b,a] = iirnotch(wo,bw);
fvtool(b,a);

y3 = filter(b,a,y2);

plot(x,y2);
hold on;
plot(x,y3);
legend('original','filtered');

figure; 
plot(x(idx:end),y2(idx:end));
figure;
plot(x(idx:end),y3(idx:end));
