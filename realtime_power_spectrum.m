%% real-time power spectrum: http://cognionics.com/wiki/pmwiki.php/Main/Real-timeDataStreamingWithLSL
close all;
clear all;

% Interval between feedback
FEEDBACK_INTERVAL = .1;
FEEDBACK_ME = 0.01;

% should be equal to FS defined in next block
sampling_rate = 128;

MV_AVG_LENGTH = 5; % in seconds
DATA_POINTS = MV_AVG_LENGTH * sampling_rate;

PSD_FREQS = 13:1:30; 

% FIGURE OUT MAPPINGS
CHANNELS_OF_INTEREST = [8, 12]; % C3 = Ch8; C4 = Ch12

data_buffer = zeros(length(CHANNELS_OF_INTEREST), DATA_POINTS); %pre-allocate data

err = 0;



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
FS = inlet.info.nominal_srate();
% apperently the sampling rate is 128 Hz

% create a figure
disp('Now receiving data...');
plot_fig = true;

%% Set up plot
% Figure
close all;
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
%% Run
while true
    % get data from the inlet
    % this keeps pulling 'samples' as though inlet is just a place to look for them
    [temp_data,ts] = inlet.pull_sample();
    
    
    % temp data has data from multiple time points (rows) for each
    % channel (columns)
    new_points2 = temp_data(:, CHANNELS_OF_INTEREST);
    new_points = new_points2.';

    % resolve length of data points
    new_length = size(new_points,2);

    % shift into buffer
    % this temporarily stores recent data, making room for the newest data
    data_buffer(:,1:DATA_POINTS-new_length) = data_buffer(:,new_length+1:end);
    data_buffer(:,DATA_POINTS-new_length+1:end) = new_points;

    % remove DC offset
    % detrend removes the mean value or linear trend from a vector or matrix, usually for FFT processing
    % it basically just transposes to zero - you can see this by comparing
    % graphs
    % takes less than 0.4 ms to run
    display_buffer = detrend(data_buffer.');

    % timing
    % NOTE: time between plots takes less than 2 ms to run
    if plot_fig == false;
        time_elapsed = toc;
        if (time_elapsed > FEEDBACK_INTERVAL + FEEDBACK_ME);
            err = err + 1;
        end;
    else 
        time_elapsed = FEEDBACK_INTERVAL;
    end

    % plot
    if (FEEDBACK_INTERVAL - FEEDBACK_ME < time_elapsed); 
        plot_fig = false;
        
%         % this 'tri-plot' takes ~70 ms to execute
%         %plot raw data
%         subplot(3,1,1);
%         time = (0:DATA_POINTS-1)/FS;
%         plot(time, display_buffer, 'k');
%         axis([0 DATA_POINTS/FS -500 500]);
%         xlabel('Time (s)');
%         ylabel('Signal (\mu V)'); % DO WE KNOW IT'S IN MICROVOLTS?
% 
%         %plot PSD
%         subplot(3,1,2);
%         [Pxx, Fxx] = pwelch(display_buffer, [], [], PSD_FREQS, FS); % There are options to play around with for this function
%         semilogy(Fxx, sqrt(Pxx), 'k', 'LineWidth', 2);
%         xlabel('Frequency (Hz)');
%         ylabel('Power (\mu V/Hz ^{1/2})');
%         axis([PSD_FREQS(1) PSD_FREQS(end) 1e-1 1e2]);
%         grid on;
% 
%         % Mean Theta Power
%         subplot(3,1,3);
%         barh(mean(Pxx));
%         set(gca, 'YTick', [])
%         ylabel('Mean Theta Power');
%         xlabel('Power (\mu V/Hz ^{1/2})');
%         axis([1e-1 1e2 .9 1.1]);
        
        % plot raw data
        % takes near 200ms to plot! 
        % CHANNEL LABELS ARE GOOD IF ORDER IS KEPT
%         i = 0;
%         for ch = CHANNELS_OF_INTEREST;
%             i = i+1;
%             subplot(length(CHANNELS_OF_INTEREST),1,i);
%             time = (0:DATA_POINTS-1)/FS;
%             plot(time, display_buffer(:,i).', 'k');
%             axis([0 DATA_POINTS/FS -500 500]);
%             ylabel(sprintf('Ch%i', ch));
%             xlabel('Time (s)');
%         end
        t1 = toc;
        i = 0;
        for ch = CHANNELS_OF_INTEREST;
            i = i+1;
            set(lines(i), 'YData', display_buffer(:,i).')
        end
        toc - t1 

% %         welch power takes ~ 4 ms
%         [Pxx, Fxx] = pwelch(display_buffer, [], [], PSD_FREQS, FS); % There are options to play around with for this function
%         
% %         Mean Theta Power
% %         plot takes ~30 ms 
%         barh(mean(mean(Pxx))); 
%         axis([1e-1 1e2 .9 1.1]);
%         set(gca,'position',[0 0 1 1],'units','normalized');
%         axis off;
        
        drawnow;
        tic;
    end;
end