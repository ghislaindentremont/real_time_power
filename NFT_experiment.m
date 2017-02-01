% Clear the workspace and the screen
sca;
close all;
clearvars;

try 

    % Here we call some default settings for setting up Psychtoolbox
    PsychDefaultSetup(2);

    % Get the screen numbers
    screens = Screen('Screens');

    % Draw to the external screen if avaliable
    screen_number = max(screens);

    % Define black and white
    white = WhiteIndex(screen_number);
    black = BlackIndex(screen_number);

    % Open an on screen window
    [window, window_rect] = PsychImaging('OpenWindow', screen_number, black);
    % HideCursor(); 

    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);

    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);

    % Query the maximum priority level
    top_priority_level = MaxPriority(window);  % I may not make use of this below.

    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(window_rect);

    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');



    %----------------------------------------------------------------------
    %                       Fixation Cross 
    %----------------------------------------------------------------------

    % Here we set the size of the arms of our fixation cross
    FIX_CROSS_DIM_PIX = 40*2;

    % Now we set the coordinates (these are all relative to zero we will let
    % the drawing routine center the cross in the center of our monitor for us)
    x_fix_coords = [-FIX_CROSS_DIM_PIX FIX_CROSS_DIM_PIX 0 0];
    y_fix_coords = [0 0 -FIX_CROSS_DIM_PIX FIX_CROSS_DIM_PIX];
    all_fix_coords = [x_fix_coords; y_fix_coords];

    % Set the line width for our fixation cross
    LINE_WIDTH_PIX = 4;

    % make it green 
    FIX_COLOR = [0 1 0];

    % just vertical line
    all_cue_coords = [[0 0 0 0]; y_fix_coords];



    %----------------------------------------------------------------------
    %                         Arrow Cue 
    %----------------------------------------------------------------------

    % arrow color - red 
    ARR_COLOR = [1 0 0];

    % draw right arrow 
    rarrow_base = [xCenter + LINE_WIDTH_PIX/2
        , yCenter - FIX_CROSS_DIM_PIX/2
        , xCenter + LINE_WIDTH_PIX/2 + FIX_CROSS_DIM_PIX
        , yCenter + FIX_CROSS_DIM_PIX/2];

    rarrow_spear = [[xCenter + LINE_WIDTH_PIX/2 + FIX_CROSS_DIM_PIX, yCenter + FIX_CROSS_DIM_PIX];
        [xCenter + LINE_WIDTH_PIX/2 + FIX_CROSS_DIM_PIX*2, yCenter];
        [xCenter + LINE_WIDTH_PIX/2 + FIX_CROSS_DIM_PIX, yCenter - FIX_CROSS_DIM_PIX]];

    % draw right arrow 
    larrow_base = [xCenter - LINE_WIDTH_PIX/2 - FIX_CROSS_DIM_PIX
        , yCenter - FIX_CROSS_DIM_PIX/2
        , xCenter - LINE_WIDTH_PIX/2 
        , yCenter + FIX_CROSS_DIM_PIX/2];

    larrow_spear = [[xCenter - LINE_WIDTH_PIX/2 - FIX_CROSS_DIM_PIX, yCenter + FIX_CROSS_DIM_PIX];
        [xCenter - LINE_WIDTH_PIX/2 - FIX_CROSS_DIM_PIX*2, yCenter];
        [xCenter - LINE_WIDTH_PIX/2 - FIX_CROSS_DIM_PIX, yCenter - FIX_CROSS_DIM_PIX]];



    %----------------------------------------------------------------------
    %                         NF Bar 
    %----------------------------------------------------------------------

    NF_COLOR = [0 0 1];
    
    LI_SCALE = 20;



    %----------------------------------------------------------------------
    %                       Timing Information
    %----------------------------------------------------------------------

    % Fixation interval time in seconds and frames
    FIX_TIME = 2;
    fix_time_frames = round(FIX_TIME / ifi);

    % Cue interval time in seconds and frames 
    CUE_TIME = 1.25;
    cue_time_frames = round(CUE_TIME / ifi);

    % Numer of frames to wait before re-drawing
    WAIT_FRAMES = 1;

    % NF time frames 
    NF_TIME = 5;
    NF_WAIT_FRAMES = 3;
    NF_time_frames = round(NF_TIME / ifi / NF_WAIT_FRAMES);

    % Baseline duration
    BASELINE_TIME = 5;
    baseline_time_frames = round(BASELINE_TIME / ifi / NF_WAIT_FRAMES);



    %----------------------------------------------------------------------
    %                       Conditions
    %----------------------------------------------------------------------

    % cue location
    cue_locs_list = {'left', 'right'};
    cue_locs = [1, 2];

    TRIALS_PER_CONDITION = 1;
    cond_matrix = repmat(cue_locs, 1, TRIALS_PER_CONDITION);

    % Get the size of the matrix
    [~, num_trials] = size(cond_matrix);

    % Randomise the conditions
    shuffler = Shuffle(1:num_trials);
    cond_matrix_shuffled2 = cond_matrix(:, shuffler);

    % get ITIs 
    itis = rand(1,num_trials)+1;

    % final condition matrix
    cond_matrix_shuffled = [cond_matrix_shuffled2; itis];



    %----------------------------------------------------------------------
    %                       EEG Aquisition Setup
    %----------------------------------------------------------------------

    % should be equal to FS defined in next block
    sampling_rate = 128;

    AVG_LENGTH = 1; % in seconds
    PAD_LENGTH = 0.5; 
    
    pad_points = round(PAD_LENGTH * sampling_rate);
    
    avg_points = round(AVG_LENGTH * sampling_rate);
    
    data_points = pad_points + avg_points;

    PSD_FREQS = 13:1:30; 

    CHANNELS_OF_INTEREST = [8, 12]; % C3 = Ch8; C4 = Ch12

    data_buffer = zeros(length(CHANNELS_OF_INTEREST), data_points); %pre-allocate data

    % make sure that everything is on the path and LSL is loaded
    addpath(genpath('liblsl-Matlab'))
    
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

    
    
    %----------------------------------------------------------------------
    %                            Filters
    %----------------------------------------------------------------------
    
    LO = 4;  % keeps 13 Hz ~0 dB
    HI = 60;
    
    % 2nd order (2n, where n is first argument)
    % second argument must be between 0 and 1 (Nyquist - FS/2)
    [b, a] = butter(1,[LO HI]/(FS/2), 'bandpass');

    % use on power lines of 60 Hz
    wo = 60/(FS/2);  
    bw = wo/35;

    [b2,a2] = iirnotch(wo,bw);
            
    

    %----------------------------------------------------------------------
    %                       Experimental Loop
    %----------------------------------------------------------------------

    for trial = 1:num_trials

        %---------------- Baseline EEG Aquisition -------------------------
        if trial == 1
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'Press Any Key To Begin The Experiment',...
            'center', 'center', white );
            Screen('Flip', window);
            KbStrokeWait; 

            % baseline EEG aquisition 
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'Please Relax With Your Eyes Open',...
            'center', 'center', white );

            power_rest_list = [];

            disp('Resting EEG Aquisition...')
            vbl = Screen('Flip', window);
            
            % Flip to the screen
            for frame = 1:baseline_time_frames-1
                
                % baseline EEG aquisition 
                Screen('TextSize', window, 36); 
                DrawFormattedText(window, 'Please Relax With Your Eyes Open',...
                'center', 'center', white );
                
                [temp_data, ts] = inlet.pull_chunk();
            
%                 [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);
                
                if size(temp_data, 2) > data_points - 1
                    new_points = temp_data(CHANNELS_OF_INTEREST, size(temp_data, 2)-data_points+1:end);
                else 
                    new_points = temp_data(CHANNELS_OF_INTEREST, :);
                end

                new_length = size(new_points,2);

                data_buffer(:,1:data_points-new_length) = data_buffer(:,new_length+1:end);
                data_buffer(:,data_points-new_length+1:end) = new_points;
              
                data_buffer2 = data_buffer.';
                
                display_buffer3 = filter(b,a, data_buffer2);  % butterworth
                display_buffer4 = filter(b2,a2, display_buffer3);  % notch filter 

                display_buffer = display_buffer4((pad_points+1):data_points, :);

%                 if frame > baseline_time_frames - 2
%                     
%                     
%                     [Pxx, Fxx] = pwelch(data_buffer2((pad_points+1):data_points, :), [], [], 1:70, FS, 'power');
%                     figure
%                     plot(Fxx, Pxx)
%                     title('before filter')
%                     
% %                     figure
% %                     plot(data_buffer2)
% %                     title('pre filter')
% %                     figure
% %                     plot(data_buffer2((pad_points+1):data_points, :))
% %                     title('pre filter padding')
% %                     
% %                     figure
% %                     plot(display_buffer3)
% %                     title('post butter')
% %                     figure
% %                     plot(display_buffer3((pad_points+1):data_points, :))
% %                     title('post butter padding')
% %                     
% %                     figure
% %                     plot(display_buffer4)
% %                     title('post notch')
% %                     figure
% %                     plot(display_buffer4((pad_points+1):data_points, :))
% %                     title('post notch padding')
%                     
%                     [Pxx, Fxx] = pwelch(display_buffer, [], [], 1:70, FS, 'power');
%                     figure
%                     plot(Fxx, Pxx)
%                     title('after filter')
%                 end

                [Pxx, Fxx] = pwelch(display_buffer, [], [], PSD_FREQS, FS, 'power');

                power_rest_list = [power_rest_list, mean(mean(Pxx))]; 
                
                % Flip to the screen
                vbl = Screen('Flip', window, vbl + (NF_WAIT_FRAMES - 0.5) * ifi);

            end

            power_rest = mean(power_rest_list); 
           
        end
        %----------------------------------------------------------------------


    %     %-------------------Trial Initiation message --------------------------
    %     Screen('TextSize', window, 36; 
    %     DrawFormattedText(window, 'Press Any Key To Begin the Trial',...
    %         'center', 'center', white );
    %     Screen('Flip', window);
    %     KbStrokeWait; 
    %     %----------------------------------------------------------------------


        %----------------------------ITI --------------------------------------
        Screen('FillRect', window, black, window_rect );

        iti = cond_matrix_shuffled(2,trial);
        iti_time_frames = round(iti / ifi / WAIT_FRAMES);
        
        vbl = Screen('Flip', window);
                
        for frame = 1:iti_time_frames-1
            Screen('FillRect', window, black, window_rect );
            vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);
        end
        %----------------------------------------------------------------------


        %--------------------- Draw the fixation cross ------------------------
        Screen('DrawLines', window, all_fix_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
        vbl = Screen('Flip', window); 

        % Flip to the screen
        for frame = 1:fix_time_frames-1
            % Draw the fixation point
            Screen('DrawLines', window, all_fix_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);

            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);
        end
        %----------------------------------------------------------------------


        %--------------------- Draw the cue -----------------------------------
        cue_loc_idx = cond_matrix_shuffled(1,trial);
        cue_loc = cue_locs_list(cue_loc_idx);

        if strcmp(cue_loc, 'right')
            arrow_base = rarrow_base;
            arrow_spear = rarrow_spear;
        elseif strcmp(cue_loc, 'left')
            arrow_base = larrow_base;
            arrow_spear = larrow_spear;
        else
            disp('Error: Cue location not specified')
            sca;
        end

        % Draw cue 
        Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
        Screen('FillRect', window, ARR_COLOR, arrow_base);
        Screen('FillPoly', window, ARR_COLOR, arrow_spear, 1);

        % Flip to the screen
        vbl = Screen('Flip', window);

        % Flip to the screen
        for frame = 1:cue_time_frames-1
            Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
            Screen('FillRect', window, ARR_COLOR, arrow_base);
            Screen('FillPoly', window, ARR_COLOR, arrow_spear, 1);

            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);
        end
        %----------------------------------------------------------------------


        %--------------------- Draw NF bar -----------------------------------    
            
        % get data from the inlet
        % this keeps pulling 'chunks' as though inlet is just a place to look for them
        [temp_data, ts] = inlet.pull_chunk();
         
        [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);
        
        LI = get_LI(cue_loc, Pxx, power_rest);
        
        NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX);
        
        % Draw NF bar 
        Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
        Screen('FillRect', window, NF_COLOR, NF_bar); 
        
        % Flip to the screen
        vbl = Screen('Flip', window);

        % Flip to the screen
        for frame = 1:NF_time_frames-1

            % get data from the inlet
            % this keeps pulling 'chunks' as though inlet is just a place to look for them
            [temp_data, ts] = inlet.pull_chunk();

            [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);
            
            LI = get_LI(cue_loc, Pxx, power_rest);
        
            NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX);
        
            % Draw NF bar
            Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
            Screen('FillRect', window, NF_COLOR, NF_bar);

            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (NF_WAIT_FRAMES - 0.5) * ifi);
            
        end
        %----------------------------------------------------------------------


        %------------------- End of Experiment --------------------------------
        if trial == 4

            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'You Have Succesfully Completed The Experiment' ,...
            'center', 'center', white );
            Screen('Flip', window);
            tic;
            while toc < 3 end

            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'The Experimenter Should Be With You Shortly',...
            'center', 'center', white );
            Screen('Flip', window);
            KbStrokeWait; 

        end
        %----------------------------------------------------------------------

    end

    % Clear the screen
    sca;
catch
    sca;
    psychrethrow(psychlasterror);
end


