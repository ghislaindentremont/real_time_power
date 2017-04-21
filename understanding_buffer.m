%% generate signal similar to Beta waves  

clear all; 
close all;

w = 23;
amplitude_rest = 10;
amplitude_NF = 2;
intercept = 4000;
noise = 2;

% set length to 1.5s, but only interested in 1s of data
pad_len = 0.5;
sig_len = 1;

len = pad_len + sig_len;

% length of half trial 
time = 5;

w_line = 60;
line_noise = 10;

FS = 128;

lo = 13;
hi = 30;

% filter for getting continous estimate of power
[b, a] = butter(3,[lo hi]/(FS/2), 'bandpass');

x = 0:(1/FS):time;

y1 = intercept + amplitude_rest * sin(2*pi*w*x) + line_noise * sin(2*pi*w_line*(x-rand)); % + noise * randn([1 (time*FS+1)]);

% add negative deflection for last 5 seconds
y2 = intercept + amplitude_NF * sin(2*pi*w*x) + line_noise * sin(2*pi*w_line*(x-rand)); % + noise * randn([1 (time*FS+1)]);

y = [y1 y2];
plot(y)

%% buffer

power_list = [];

for i=(len*FS+1):length(y)
    
    data_buffer = y(i-(len*FS-1):i);
    
    ye = filter(b,a,data_buffer);
    
    ye2 = ye.^2;
    
    power_list = [power_list, mean(ye2(pad_len*FS+1:len*FS))];
    
end


%% baseline power 

% we'll take a couple second period

disp('mean of power list')
mean(power_list(1:FS*2))

% just filter once and compute 
ye = filter(b,a,y);
ye2 = ye.^2;

disp('mean power')
mean(ye2((len*FS+1):(FS*2+len*FS)))
% the power values are similar but not exactly the same because of buffer


%% log ERS

ERS_list = power_list((1282/2-len*FS+1):end)/mean(power_list(1:FS*2));
log_ERS_list = log(ERS_list);
disp('mean log ERS');
mean(log_ERS_list)

mean_ERS = mean(ye2((time*FS+2):end))/mean(power_list(1:FS*2));
disp('log mean ERS');
log(mean_ERS)

plot(ye2((len*FS+1):end))

