%% real-time power spectrum: http://cognionics.com/wiki/pmwiki.php/Main/Real-timeDataStreamingWithLSL
close all;
clear all;

% Interval between feedback
% interval of 0.01 with ME of 0.001 is good for raw signals 
% interval of 0.02 with ME of 0.001 is good for power
% signal
FEEDBACK_INTERVAL = .02; 
FEEDBACK_ME = 0.001;

% rest aquisition time
REST_TIME = 10;

% should be equal to FS defined in next block
sampling_rate = 128;

MV_AVG_LENGTH = 1; % in seconds
DATA_POINTS = MV_AVG_LENGTH * sampling_rate;

PSD_FREQS = 13:1:30; 

% FIGURE OUT MAPPINGS
CHANNELS_OF_INTEREST = [8, 12]; % C3 = Ch8; C4 = Ch12

data_buffer = zeros(length(CHANNELS_OF_INTEREST), DATA_POINTS); %pre-allocate data

% what display?
power = true;

% what hand
hand_right = true;

%% instantiate the LSL library
% make sure that everything is on the path and LSL is loaded
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an EEG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'name','openvibeSignal'); end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

% resolve stream sample rate
FS = inlet.info.nominal_srate(); % apperently the sampling rate is 128 Hz

% create a figure
disp('Now receiving data...');
plot_fig = true;

close all;
figure('MenuBar', 'none', 'WindowStyle', 'modal');
if power
    b = barh(1e1); 
    axis([-1e2 1e2 .9 1.1]);
    set(gca,'position',[0 0 1 1],'units','normalized');
    axis off; 
else
    i = 0;
    for ch = CHANNELS_OF_INTEREST;
        i = i+1;
        ax(i) = subplot(length(CHANNELS_OF_INTEREST),1,i);

        set(ax(i), 'xlim', [0 DATA_POINTS/FS], 'ylim', [-500 500]);
        ax(i).YLabel.String =  sprintf('Ch%i', ch);
        ax(i).XLabel.String =  'Time (s)';
        lines(i) = line('Parent', ax(i));
    end;

    % Create xlabel
    time = (0:DATA_POINTS-1)/FS;

    % call plots early
    i = 0;
    for ch = CHANNELS_OF_INTEREST;
        i = i+1;
        set(lines(i), 'XData', time, 'YData', zeros(1,DATA_POINTS));
    end;
end

% Run
tic;
k = 0;
j = 0;

inlet_times = 0;
buffer_times = 0;
detrend_times = 0;
plot_times = 0;
since_last_plot_times = 0;
between_plot_times = 0;

plot_stop = 0;
while true
    k = k + 1;
    
    if k ==1
        disp('Resting EEG Aquisition...')
        start_rest = toc; 
        power_contra_rest_list = [];
        power_ipsi_rest_list = [];
        while toc - start_rest < REST_TIME
            [temp_data, ts] = inlet.pull_chunk();

            new_points = temp_data(CHANNELS_OF_INTEREST, :);
            new_length = size(new_points,2);

            data_buffer(:,1:DATA_POINTS-new_length) = data_buffer(:,new_length+1:end);
            data_buffer(:,DATA_POINTS-new_length+1:end) = new_points;

            display_buffer = detrend(data_buffer.');

            [Pxx, Fxx] = pwelch(display_buffer, [], [], PSD_FREQS, FS, 'power');
            if hand_right;
                power_contra_rest_list = [power_contra_rest_list, mean(Pxx(:,1)) ]; % Ch8 ==> C3
                power_ipsi_rest_list = [power_ipsi_rest_list, mean(Pxx(:,2)) ]; % Ch12 ==> C4
            end;
            
            if ( mod(toc - start_rest, 1) > 1-FEEDBACK_INTERVAL )
                disp(floor(toc - start_rest) )
                pause(FEEDBACK_INTERVAL)
            end
        end
        power_contra_rest = mean(power_contra_rest_list);
        power_ipsi_rest = mean(power_ipsi_rest_list);
    end
    
    % get data from the inlet
    % this keeps pulling 'chunks' as though inlet is just a place to look for them
    start_inlet = toc;
    [temp_data, ts] = inlet.pull_chunk();
    inlet_time = toc - start_inlet;
    inlet_times(k) = inlet_time;
    
    % temp data has data from multiple time points (rows) for each
    % channel (columns)
    start_buffer = toc;
    new_points = temp_data(CHANNELS_OF_INTEREST, :);
    % resolve length of data points
    new_length = size(new_points,2);
    % shift into buffer
    % this temporarily stores recent data, making room for the newest data
    data_buffer(:,1:DATA_POINTS-new_length) = data_buffer(:,new_length+1:end);
    data_buffer(:,DATA_POINTS-new_length+1:end) = new_points;
    buffer_time = toc - start_buffer;
    buffer_times(k) = buffer_time;

    % remove DC offset
    % detrend removes the mean value or linear trend from a vector or matrix, usually for FFT processing
    % it basically just transposes to zero - you can see this by comparing
    % graphs
    start_detrend = toc;
    display_buffer = detrend(data_buffer.');
    detrend_time = toc - start_detrend;
    detrend_times(k) = detrend_time;
   

    % timing
    if plot_fig == false;
        since_last_plot_time = toc - plot_stop;
    else 
        since_last_plot_time = FEEDBACK_INTERVAL;
    end
    since_last_plot_times(k) = since_last_plot_time;

    % plot
    if (FEEDBACK_INTERVAL - FEEDBACK_ME < since_last_plot_time); 
        j = j + 1;
        plot_fig = false;
       
        between_plot_time = toc - plot_stop;
        between_plot_times(j) = between_plot_time;
        
        if power
            plot_start = toc;
    %         welch power takes ~ 4 ms
            [Pxx, Fxx] = pwelch(display_buffer, [], [], PSD_FREQS, FS, 'power');
            if hand_right;
                power_contra = mean(Pxx(:,1)); % Ch8 ==> C3
                log2_ERS_contra = log2(power_contra/power_contra_rest); % where positive means syncronisation (increased power relative to baseline)
                
                power_ipsi = mean(Pxx(:,2)); % Ch12 ==> C4
                log2_ERS_ipsi = log2(power_ipsi/power_ipsi_rest);
                
                LI = log2_ERS_contra - log2_ERS_ipsi;
            end;
                
    %         Mean Theta Power
            set(b, 'YData', LI)
            plot_stop = toc;
            plot_time = plot_stop - plot_start;
            plot_times(j) = plot_time;
        else
            plot_start = toc;
            i = 0;
            for ch = CHANNELS_OF_INTEREST;
                i = i+1;
                set(lines(i), 'YData', display_buffer(:,i).')
            end
            plot_stop = toc;
            plot_time = plot_stop - plot_start;
            plot_times(j) = plot_time;
        end
      
        drawnow;
    end;
end

%% Look at timing after
close all;

% figure
% histogram(inlet_times)
% title('inlet times')

figure
histogram(plot_times)
title('plot times')

% figure
% histogram(detrend_times)
% title('detrend times')
% 
% figure
% histogram(buffer_times)
% title('buffer times')

figure
histogram(since_last_plot_times)
title('since last plot times')

figure
histogram(between_plot_times)
title('between plot times')