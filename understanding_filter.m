%% generate signal similar to Beta waves  

clear all; 
close all;

w = 23;
amplitude = 4;
intercept = 4000;
% set length to 2s, but only interested in 1s of data
pad_len = 1;
sig_len = 1;

len = pad_len + sig_len;

w_line = 60;
line_noise = 10;

FS = 128;
lo = 4;
hi = 45;

[b, a] = butter(1,[lo hi]/(FS/2), 'bandpass');
fvtool(b,a);

x = 0:(1/FS):len;

y = intercept + amplitude * sin(2*pi*w*x) + line_noise * sin(2*pi*w_line*(x-rand));


%% get power of initial signal

PSD_FREQS = 8:1:30; 
[Pxx1, Fxx1] = pwelch(y, [], [], PSD_FREQS, FS, 'power');

plot(Fxx1, Pxx1)
 

%% butterworth filter 
y2 = filter(b,a,y);

plot(x,y);
hold on;
plot(x,y2);
title('with padding')
legend('original','filtered');

figure;
idx = length(0:(1/FS):pad_len);
plot(x(idx:end),y(idx:end));
title('original')

figure;
plot(x(idx:end),y2(idx:end));
title('filtered')


%% notch 
% use on power lines of 60 Hz
wo = 60/(FS/2);  
bw = wo/35;

[b2,a2] = iirnotch(wo,bw);
fvtool(b2,a2);

y3 = filter(b2,a2,y2);

plot(x,y2);
hold on;
plot(x,y3);
legend('original','filtered');

figure; 
plot(x(idx:end),y2(idx:end));
title('original')

figure;
plot(x(idx:end),y3(idx:end));
title('filtered')



%% Get power of final signal
[Pxx, Fxx] = pwelch(y3, [], [], PSD_FREQS, FS, 'power');
plot(Fxx, Pxx)