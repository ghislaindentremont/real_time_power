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
    %                        ID/Demographics
    %----------------------------------------------------------------------
    
    Screen('TextSize', window, 36); 
    DrawFormattedText(window, 'But really, enough about me, let''s talk about you...',...
    'center', 'center', white );
    Screen('Flip', window);

    prompt = {'Enter participant ID:'
        , 'Enter participant age:'
        ,'Enter participant sex (''m'' or ''f''):'
        , 'Enter participant handedness (''r'' or ''l''):'
        };
    dlg_title = 'Demographics';
    num_lines = 1;
    dems = inputdlg(prompt,dlg_title,num_lines);
    id = char(dems(1));
    age = char(dems(2));
    sex = char(dems(3));
    hand = char(dems(4));

    while (strcmp(hand, 'r') ~= 1 && strcmp(hand, 'l') ~= 1) ||  (strcmp(sex, 'm') ~= 1 && strcmp(sex, 'f') ~= 1)
        dems = inputdlg(prompt,dlg_title,num_lines);
        id = char(dems(1));
        age = char(dems(2));
        sex = char(dems(3));
        hand = char(dems(4));
    end
    
    % get numbers
    id = str2num(id);
    age = str2num(age);
       
    if (strcmp(sex, 'm') == 1)
        sex = 1;
    elseif (strcmp(sex, 'f') == 1)
        sex = 2;
    else
        disp('ERROR: sex');  % should also get and error when putting in response matrix
    end
    
    if (strcmp(hand, 'l') == 1)
        hand = 1;
    elseif (strcmp(hand, 'r') == 1)
        hand = 2;
    else
        disp('ERROR: handedness');  % should also get and error when putting in response matrix
    end  
    
    % time and date 
    format shortg
    c = clock;
    year = c(1);
    month = c(2);
    day = c(3);
    hour = c(4);
    minute = c(5);
    seconds = c(6);
    


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

    % Numer of frames to wait before re-drawing
    WAIT_FRAMES = 1;
    
    % NF time frames 
    NF_TIME = 5;
    NF_WAIT_FRAMES = 3;
    NF_time_frames = round(NF_TIME / ifi / NF_WAIT_FRAMES);
    
    % Fixation interval time in seconds and frames
    FIX_TIME = 2;
    fix_time_frames = round(FIX_TIME / ifi / NF_WAIT_FRAMES);

    % Cue interval time in seconds and frames 
    CUE_TIME = 1.25;
    cue_time_frames = round(CUE_TIME / ifi/ WAIT_FRAMES);



    %----------------------------------------------------------------------
    %                       Conditions
    %----------------------------------------------------------------------

    % block condition
    blk_conds_list = {'MI', 'ME'};
    blk_conds = [1, 2];
    
    BLOCKS_PER_CONDITION = 1;
    blk_matrix = repmat(blk_conds, 1, BLOCKS_PER_CONDITION);
    
    % Get the size of the matrix
    [~, num_blocks] = size(blk_matrix);

    % Randomise the conditions
    shuffler = Shuffle(1:num_blocks);
    blk_matrix_shuffled = blk_matrix(:, shuffler);
    
    % cue location
    cue_locs_list = {'left', 'right'};
    cue_locs = [1, 2];

    TRIALS_PER_CONDITION = 1;
    cond_matrix = repmat(cue_locs, 1, TRIALS_PER_CONDITION);

    % Get the size of the matrix
    [~, num_trials] = size(cond_matrix);
    
    

    %----------------------------------------------------------------------
    %                       EEG Aquisition Setup
    %----------------------------------------------------------------------

    power_rest_window = [];
    
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

    for block = 1:num_blocks
        
        block_cond = blk_matrix_shuffled(block);
        
        % Randomize the conditions
        shuffler = Shuffle(1:num_trials);
        cond_matrix_shuffled2 = cond_matrix(:, shuffler);

        % get ITIs   
        itis = rand(1,num_trials)+1;

        % final condition matrix
        cond_matrix_shuffled = [cond_matrix_shuffled2; itis];
        
        if block_cond == 1
            
            tic;
            while toc < 5
                %------------------- Block Instruction Message ------------------------
                Screen('TextSize', window, 36); 
                DrawFormattedText(window, 'Motor Imagery: You will imagine moving\nthe hand indicated by the arrow in each trial\nwithout actually moving that hand.',...
                    'center', 'center', white );
                Screen('Flip', window);
                %----------------------------------------------------------------------
            end
            
            %------------------- Block Instruction Message ------------------------
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'Motor Imagery: You will imagine moving\nthe hand indicated by the arrow in each trial\nwithout actually moving that hand.\n\n\nPress Any Key To Begin the Block',...
                'center', 'center', white );
            Screen('Flip', window);
            KbStrokeWait; 
            %----------------------------------------------------------------------
            
        elseif block_cond == 2
            
            tic;
            while toc < 5
                %------------------- Block Instruction Message ------------------------
                Screen('TextSize', window, 36); 
                DrawFormattedText(window, 'Motor Execution: You will move\nthe hand indicated by the arrow in each trial.',...
                    'center', 'center', white );
                Screen('Flip', window);
                %----------------------------------------------------------------------
            end
            
            %------------------- Block Instruction Message ------------------------
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'Motor Execution: You will move\nthe hand indicated by the arrow in each trial.\n\n\nPress Any Key To Begin the Block',...
                'center', 'center', white );
            Screen('Flip', window);
            KbStrokeWait; 
            %----------------------------------------------------------------------
            
        else
            disp('ERROR: block condition not defined')       
        end
        
        
        for trial = 1:num_trials
          

            %-------------------Trial Initiation message --------------------------
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'Press Any Key To Begin the Trial',...
                'center', 'center', white );
            Screen('Flip', window);
            KbStrokeWait; 
            %----------------------------------------------------------------------


            %----------------------------ITI --------------------------------------
            iti = cond_matrix_shuffled(2,trial);
            iti_time_frames = round(iti / ifi / WAIT_FRAMES);

            % warm up data buffer 
            [temp_data, ts] = inlet.pull_chunk();


            %------------ Save Data -----------%
            trial_stage = 1;
            iteration = 1;
            cue_loc_idx = NaN;

            nrow = size(temp_data,2);
            raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            if trial == 1
                raw_mat = raw_mat_temp;
            else
                raw_mat = [raw_mat; raw_mat_temp];
            end
            %----------------------------------%


            [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


            %------------ Save Data -----------%
            % Pxx(:,1) is C3 whereas P(:,2) is C4 (i.e. left to right)
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN];

            if trial == 1
                pwr_mat = pwr_mat_temp;
            else
                pwr_mat = [pwr_mat; pwr_mat_temp];
            end
            %----------------------------------%


            Screen('FillRect', window, black, window_rect );
            vbl = Screen('Flip', window);

            for frame = 1:iti_time_frames-1
                % warm up data buffer 
                [temp_data, ts] = inlet.pull_chunk();

                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


                 %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


                Screen('FillRect', window, black, window_rect );
                vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);
            end
            %----------------------------------------------------------------------


            %--------------------- Draw the fixation cross ------------------------
            % set power to zero
            power_rest_list = [];

            [temp_data, ts] = inlet.pull_chunk();


            %------------ Save Data -----------%
            trial_stage = 2;
            iteration = 1;

            nrow = size(temp_data,2);
            raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            raw_mat = [raw_mat; raw_mat_temp];
            %----------------------------------%


            [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


            %------------ Save Data -----------%
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN];
            pwr_mat = [pwr_mat; pwr_mat_temp];
            %----------------------------------%


            power_rest_list = [power_rest_list, mean(mean(Pxx))]; 

            Screen('DrawLines', window, all_fix_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
            vbl = Screen('Flip', window); 

            % Flip to the screen
            for frame = 1:fix_time_frames-1

                [temp_data, ts] = inlet.pull_chunk();


                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


                power_rest_list = [power_rest_list, mean(mean(Pxx))]; 

                Screen('DrawLines', window, all_fix_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                vbl = Screen('Flip', window, vbl + (NF_WAIT_FRAMES - 0.5) * ifi);
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

            [temp_data, ts] = inlet.pull_chunk();


            %------------ Save Data -----------%
            trial_stage = 3;
            iteration = 1;

            nrow = size(temp_data,2);
            raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            raw_mat = [raw_mat; raw_mat_temp];
            %----------------------------------%


            % keep updating data buffer, but not baseline/rest power...
            [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


            %------------ Save Data -----------%
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN];
            pwr_mat = [pwr_mat; pwr_mat_temp];
            %----------------------------------%


            Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
            Screen('FillRect', window, ARR_COLOR, arrow_base);
            Screen('FillPoly', window, ARR_COLOR, arrow_spear, 1);
            vbl = Screen('Flip', window);

            % Flip to the screen
            for frame = 1:cue_time_frames-1

                [temp_data, ts] = inlet.pull_chunk();


                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];
                
                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                % keep updating data buffer, but not baseline/rest power...
                [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


                Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                Screen('FillRect', window, ARR_COLOR, arrow_base);
                Screen('FillPoly', window, ARR_COLOR, arrow_spear, 1);
                vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);
            end
            %----------------------------------------------------------------------


            %--------------------- Draw NF bar -----------------------------------        
            [temp_data, ts] = inlet.pull_chunk();


            %------------ Save Data -----------%
            trial_stage = 4;
            iteration = 1;

            nrow = size(temp_data,2);
            raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            raw_mat = [raw_mat; raw_mat_temp];
            %----------------------------------%


            [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);
            
            power_rest = mean(power_rest_list);

            power_rest_window = [power_rest_window power_rest];

            if length(power_rest_window) < 40
                power_rest_mavg = mean(power_rest_window);
            else 
                power_rest_mavg = mean(power_rest_window(end-40+1:end));
            end

            LI = get_LI(cue_loc, Pxx, power_rest_mavg);


            %------------ Save Data -----------%
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration power_rest_mavg mean(Pxx(:,1)) mean(Pxx(:,2)) LI];
            pwr_mat = [pwr_mat; pwr_mat_temp];
            %----------------------------------%


            NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX);

            % Draw NF bar 
            Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
            Screen('FillRect', window, NF_COLOR, NF_bar); 

            % Flip to the screen
            vbl = Screen('Flip', window);

            % Flip to the screen
            for frame = 1:NF_time_frames-1

                [temp_data, ts] = inlet.pull_chunk();


                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);

                LI = get_LI(cue_loc, Pxx, power_rest);


                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block block_cond trial cue_loc_idx iti trial_stage iteration power_rest mean(Pxx(:,1)) mean(Pxx(:,2)) LI];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


                NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX);

                % Draw NF bar
                Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                Screen('FillRect', window, NF_COLOR, NF_bar);

                % Flip to the screen
                vbl = Screen('Flip', window, vbl + (NF_WAIT_FRAMES - 0.5) * ifi);

            end
            %----------------------------------------------------------------------

        end
        
    end
    
    
    %------------------- End of Experiment --------------------------------
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
    %----------------------------------------------------------------------
    
    
    % write response matrix to csv
    csvwrite(sprintf('C:/Users/Kine Research/Documents/MATLAB/ghis_data/%i_raw.csv', id), raw_mat);
    csvwrite(sprintf('C:/Users/Kine Research/Documents/MATLAB/ghis_data/%i_pwr.csv', id), pwr_mat);
  
    % Clear the screen
    sca;
catch
    sca;
    psychrethrow(psychlasterror);
end


