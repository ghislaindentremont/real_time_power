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
        , 'Enter task condition (''MI'' or ''ME''):'
        , 'Enter feedback condition (''C'', ''I'')'
        };
    dlg_title = 'Demographics';
    num_lines = 1;
    dems = inputdlg(prompt,dlg_title,num_lines);
    
    if strcmp(char(dems(1)), 'test')
            
        id = 99;
        age = 99;
        sex = 99;
        hand = 99; 
        
        task_str = 'MI';
        task = 1;
        feedback_str = 'C';
        feedback = 1;
        
    else
        
        id = char(dems(1));
        age = char(dems(2));
        sex = char(dems(3));
        hand = char(dems(4));
        task_str = char(dems(5));
        feedback_str = char(dems(6));

        while (strcmp(hand, 'r') ~= 1 && strcmp(hand, 'l') ~= 1) ||  (strcmp(sex, 'm') ~= 1 && strcmp(sex, 'f') ~= 1) ||  (strcmp(task_str, 'MI') ~= 1 && strcmp(task_str, 'ME') ~= 1) ||  (strcmp(feedback_str, 'C') ~= 1 && strcmp(feedback_str, 'I') ~= 1)
            dems = inputdlg(prompt,dlg_title,num_lines);
            id = char(dems(1));
            age = char(dems(2));
            sex = char(dems(3));
            hand = char(dems(4));
            task_str = char(dems(5));
            feedback_str = char(dems(6));
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

        if (strcmp(task_str, 'MI') == 1)
            task = 1;
        elseif (strcmp(task_str, 'ME') == 1)
            task = 2;
        else
            disp('ERROR: block condition');  % should also get and error when putting in response matrix
        end  
        
        if (strcmp(task_str, 'C') == 1)
            feedback = 1;
        elseif (strcmp(task_str, 'I') == 1)
            feedback = 2;
        else
            disp('ERROR: block condition');  % should also get and error when putting in response matrix
        end  
        
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
    
    
    % set escape
    KbName('UnifyKeyNames');
    escapeKey = KbName('ESCAPE');

    keysOfInterest=zeros(1,256);
    keysOfInterest(KbName({'ESCAPE'}))=1;
    KbQueueCreate(-1, keysOfInterest);
    % start 'escape' Queue
    KbQueueStart;


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
    
    % just horizontal lines
    all_NF_coords = [[-FIX_CROSS_DIM_PIX FIX_CROSS_DIM_PIX 0 0]; [0 0 0 0]];


    %----------------------------------------------------------------------
    %                         Arrow Cue 
    %----------------------------------------------------------------------

    ARR_COLOR = [0 0 1];

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
    NF_COLOR_CONTRA = [0 1 0]; % want contralateral ERD to go up
    NF_COLOR_IPSI = [1 0 0];  % want ipsilateral ERD to go up less than contra
    
    LI_SCALE = 20;
    LOG_ERS_SCALE = 20;
    
    % for Intermittent condition
    holder_bar_right = [xCenter + xCenter/4 - FIX_CROSS_DIM_PIX/2
        , yCenter + LINE_WIDTH_PIX/2 
        , xCenter + xCenter/4 + FIX_CROSS_DIM_PIX/2
        , yCenter +  LINE_WIDTH_PIX/2 + yCenter/4];

    holder_bar_left = [xCenter - xCenter/4 - FIX_CROSS_DIM_PIX/2
        , yCenter + LINE_WIDTH_PIX/2 
        , xCenter - xCenter/4 + FIX_CROSS_DIM_PIX/2;
        , yCenter + LINE_WIDTH_PIX/2 + yCenter/4];


    
    %----------------------------------------------------------------------
    %                       Timing Information
    %----------------------------------------------------------------------

    % Numer of frames to wait before re-drawing
    WAIT_FRAMES = 3;
    
    % NF time frames 
    NF_TIME = 5;
    NF_time_frames = round(NF_TIME / ifi / WAIT_FRAMES);
    
    % Fixation interval time in seconds and frames
    FIX_TIME = 2;
    fix_time_frames = round(FIX_TIME / ifi / WAIT_FRAMES);

    % Cue interval time in seconds and frames 
    CUE_TIME = 1.25;
    cue_time_frames = round(CUE_TIME / ifi/ WAIT_FRAMES);



    %----------------------------------------------------------------------
    %                       Conditions
    %----------------------------------------------------------------------
    
    NUM_EXPERIMENTAL_BLOCKS = 4;
    
    num_blocks = NUM_EXPERIMENTAL_BLOCKS + 1;
    
    % cue location
    cue_locs_list = {'left', 'right'};
    cue_locs = [1, 2];
    
    NUM_PRACTICE_TRIALS = 10;

    TRIALS_PER_CONDITION = 30;
    cond_matrix = repmat(cue_locs, 1, TRIALS_PER_CONDITION);

    % Get the size of the matrix
    [~, num_trials] = size(cond_matrix);
    
    

    %----------------------------------------------------------------------
    %                       EEG Aquisition Setup
    %----------------------------------------------------------------------

    power_rest_window = [];
    WINDOW_LENGTH = 1;
    
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
    
    LO = 13;  
    HI = 30;
    
    % 6th order (2n, where n is first argument)
    % second argument must be between 0 and 1 (Nyquist - FS/2)
    [b, a] = butter(3,[LO HI]/(FS/2), 'bandpass');

    % use on power lines of 60 Hz
    wo = 60/(FS/2);  
    bw = wo/35;

    [b2,a2] = iirnotch(wo,bw);
            
    

    %----------------------------------------------------------------------
    %                       Experimental Loop
    %----------------------------------------------------------------------

    [pressed, firstPress]=KbQueueCheck; 
    
    for block = 1:num_blocks
        
        if pressed
            if firstPress(KbName('ESCAPE'))
                break
            end
        end
        
        num_trials_in_this_block = num_trials;
        
        % Randomize the conditions
        shuffler = Shuffle(1:num_trials);
        cond_matrix_shuffled2 = cond_matrix(:, shuffler);

        % get ITIs   
        itis = rand(1,num_trials)+1;

        % final condition matrix
        cond_matrix_shuffled = [cond_matrix_shuffled2; itis];
        
        if block == 1
            
            num_trials_in_this_block = NUM_PRACTICE_TRIALS;
           
            if task == 1

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
                DrawFormattedText(window, 'Motor Imagery: You will imagine moving\nthe hand indicated by the arrow in each trial\nwithout actually moving that hand.\n\n\nPress Any Key To Move On',...
                    'center', 'center', white );
                Screen('Flip', window);
%                 KbStrokeWait; 
                %----------------------------------------------------------------------
                
                if feedback == 1
                    
                    %------------------- Block Instruction Message ------------------------
                    Screen('TextSize', window, 26); 
                    DrawFormattedText(window, 'Continuous Feedback: During the motor imagery portion of each trial\nyou will be presented with continous feedback\nregarding your performance. The feedback will be displayed in the form of two continously changing vertical bars.\nYour goal will be to increase the bar on the side that the arrow pointed\nwhile keeping the other bar near zero.\n\n\nPress Any Key To Begin the Block',...
                        'center', 'center', white );
                    Screen('Flip', window);
    %                 KbStrokeWait; 
                    %----------------------------------------------------------------------
                    
                elseif feedback == 2
                    
                    %------------------- Block Instruction Message ------------------------
                    Screen('TextSize', window, 26); 
                    DrawFormattedText(window, 'Intermittent Feedback: During the motor imagery portion of each trial\nyou will be presented with intermittent feedback at the end of each trial\nregarding your performance. The feedback will be displayed in the form of two vertical bars.\nYour goal, across trials, will be to increase the bar on the side that the arrow pointed\nwhile keeping the other bar near zero.\n\n\nPress Any Key To Begin the Block',...
                        'center', 'center', white );
                    Screen('Flip', window);
    %                 KbStrokeWait; 
                    %----------------------------------------------------------------------

                else
                    disp('ERROR: feedback condition not defined') 
                end

            elseif task == 2

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
%                 KbStrokeWait; 
                %----------------------------------------------------------------------
                
                if feedback == 1
                    
                    %------------------- Block Instruction Message ------------------------
                    Screen('TextSize', window, 26); 
                    DrawFormattedText(window, 'Continuous Feedback: During the motor execution portion of each trial\nyou will be presented with continous feedback\nregarding your performance. The feedback will be displayed in the form of two continously changing vertical bars.\nYour goal will be to increase the bar on the side that the arrow pointed\nwhile keeping the other bar near zero.\n\n\nPress Any Key To Begin the Block',...
                        'center', 'center', white );
                    Screen('Flip', window);
    %                 KbStrokeWait; 
                    %----------------------------------------------------------------------
                    
                elseif feedback == 2
                    
                    %------------------- Block Instruction Message ------------------------
                    Screen('TextSize', window, 26); 
                    DrawFormattedText(window, 'Intermittent Feedback: During the motor ececution portion of each trial\nyou will be presented with intermittent feedback at the end of each trial\nregarding your performance. The feedback will be displayed in the form of two vertical bars.\nYour goal, across trials, will be to increase the bar on the side that the arrow pointed\nwhile keeping the other bar near zero.\n\n\nPress Any Key To Begin the Block',...
                        'center', 'center', white );
                    Screen('Flip', window);
    %                 KbStrokeWait; 
                    %----------------------------------------------------------------------

                else
                    disp('ERROR: feedback condition not defined') 
                end

            else
                disp('ERROR: task condition not defined')       
            end
            
            tic;
            while toc < 5
                %------------------- Block Instruction Message ------------------------
                Screen('TextSize', window, 36); 
                DrawFormattedText(window, 'This is a practice block.',...
                    'center', 'center', white );
                Screen('Flip', window);
                %----------------------------------------------------------------------
            end
            
            %------------------- Block Instruction Message ------------------------
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'This is a practice block.\n\n\nPress Any Key To Begin the Block',...
                'center', 'center', white );
            Screen('Flip', window);
%             KbStrokeWait; 
            %----------------------------------------------------------------------
            
        elseif block == 2
            %------------------- Block Instruction Message ------------------------
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'This is an experimental block.\n\n\nPress Any Key To Begin the Block',...
                'center', 'center', white );
            Screen('Flip', window);
%             KbStrokeWait; 
            %----------------------------------------------------------------------
        else
            %------------------- Block Instruction Message ------------------------
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'Take a break!\n\n\nPress Any Key To Begin the Block',...
                'center', 'center', white );
            Screen('Flip', window);
%             KbStrokeWait; 
            %----------------------------------------------------------------------
        end
        
        
        
        for trial = 1:num_trials_in_this_block
            
            % check if escape key has been pressed. If so, exit experiment 
            [pressed, firstPress]=KbQueueCheck; 
            if pressed
                if firstPress(KbName('ESCAPE'))
                    break
                end
            end

%             -------------------Trial Initiation message --------------------------
%             Screen('TextSize', window, 36); 
%             DrawFormattedText(window, 'Press Any Key To Begin the
%             Trial',...
%                 'center', 'center', white );
%             Screen('Flip', window);
%             KbStrokeWait; 
%             ----------------------------------------------------------------------

            cue_loc_idx = cond_matrix_shuffled(1,trial);  % define cue location id ahead of time


            %----------------------------ITI --------------------------------------
            iti = cond_matrix_shuffled(2,trial);
            iti_time_frames = round(iti / ifi / WAIT_FRAMES);

            % warm up data buffer 
            [temp_data, ts] = inlet.pull_chunk();


            %------------ Save Data -----------%
            trial_stage = 1;
            iteration = 1;

            nrow = size(temp_data,2);
            raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            if trial == 1 
                raw_mat = raw_mat_temp;
            else
                raw_mat = [raw_mat; raw_mat_temp];
            end
            %----------------------------------%


            [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


            %------------ Save Data -----------%
            % Pxx(:,1) is C3 whereas P(:,2) is C4 (i.e. left to right)
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN NaN NaN];

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
                raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


                 %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN NaN NaN];
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
            raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            raw_mat = [raw_mat; raw_mat_temp];
            %----------------------------------%


            [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


            %------------ Save Data -----------%
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN NaN NaN];
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
                raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN NaN NaN];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


                power_rest_list = [power_rest_list, mean(mean(Pxx))]; 

                Screen('DrawLines', window, all_fix_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);
            end
            %----------------------------------------------------------------------


            %--------------------- Draw the cue -----------------------------------
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
            raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            raw_mat = [raw_mat; raw_mat_temp];
            %----------------------------------%


            % keep updating data buffer, but not baseline/rest power...
            [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


            %------------ Save Data -----------%
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN NaN NaN];
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
                raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];
                
                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                % keep updating data buffer, but not baseline/rest power...
                [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);


                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration NaN mean(Pxx(:,1)) mean(Pxx(:,2)) NaN NaN NaN];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


                Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                Screen('FillRect', window, ARR_COLOR, arrow_base);
                Screen('FillPoly', window, ARR_COLOR, arrow_spear, 1);
                vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);
            end
            %----------------------------------------------------------------------


            %--------------------- Draw NF bar ----------------------------------- 
            if feedback == 2
               log_ERS_ipsi_list = [];
               log_ERS_contra_list = [];
            end
            
            [temp_data, ts] = inlet.pull_chunk();


            %------------ Save Data -----------%
            trial_stage = 4;
            iteration = 1;

            nrow = size(temp_data,2);
            raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
            raw_info_mat = repmat(raw_info, nrow, 1);

            raw_eeg_mat = [ts.' temp_data.'];

            raw_mat_temp = [raw_info_mat raw_eeg_mat];

            raw_mat = [raw_mat; raw_mat_temp];
            %----------------------------------%


            [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);
            
            power_rest = mean(power_rest_list);

            power_rest_window = [power_rest_window power_rest];

            if length(power_rest_window) < WINDOW_LENGTH
                power_rest_mavg = mean(power_rest_window);
            else 
                power_rest_mavg = mean(power_rest_window(end-WINDOW_LENGTH+1:end));
            end

            [LI, log_ERS_ipsi, log_ERS_contra] = get_LI(cue_loc, Pxx, power_rest_mavg);
            
                   
            log_ERS_ipsi_list = [log_ERS_ipsi_list log_ERS_ipsi];
            log_ERS_contra_list = [log_ERS_contra_list log_ERS_contra];
            

            %------------ Save Data -----------%
            pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration power_rest_mavg mean(Pxx(:,1)) mean(Pxx(:,2)) LI log_ERS_ipsi log_ERS_contra];
            pwr_mat = [pwr_mat; pwr_mat_temp];
            %----------------------------------%


            if strcmp(cue_loc, 'right')
                NF_color_right = NF_COLOR_CONTRA;
                NF_color_left = NF_COLOR_IPSI;
            elseif strcmp(cue_loc, 'left')
                NF_color_right = NF_COLOR_IPSI;
                NF_color_left = NF_COLOR_CONTRA;
            end
            
%             NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc);
            [NF_bar_left, NF_bar_right] = get_NF_bars(log_ERS_ipsi, log_ERS_contra, LOG_ERS_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc);
           
            
            if feedback == 1
                
    %             Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter+xCenter/4 yCenter], 2);
                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter-xCenter/4 yCenter], 2);
    %             Screen('FillRect', window, NF_COLOR, NF_bar); 
                Screen('FillRect', window, NF_color_left, NF_bar_left); 
                Screen('FillRect', window, NF_color_right, NF_bar_right);
                
            elseif feedback == 2
                
                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter+xCenter/4 yCenter], 2);
                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter-xCenter/4 yCenter], 2);
                Screen('FillRect', window, NF_color_left, holder_bar_left); 
                Screen('FillRect', window, NF_color_right, holder_bar_right);
                
            else 
                disp('ERROR: feedback condition undefined')
            end

            
            % Flip to the screen
            vbl = Screen('Flip', window);

            % Flip to the screen
            for frame = 1:NF_time_frames-1

                [temp_data, ts] = inlet.pull_chunk();


                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%


                [Pxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2);

                [LI, log_ERS_ipsi, log_ERS_contra] = get_LI(cue_loc, Pxx, power_rest);
                
                
                log_ERS_ipsi_list = [log_ERS_ipsi_list log_ERS_ipsi];
                log_ERS_contra_list = [log_ERS_contra_list log_ERS_contra];


                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration power_rest mean(Pxx(:,1)) mean(Pxx(:,2)) LI log_ERS_ipsi log_ERS_contra];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


%                 NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc);
                [NF_bar_left, NF_bar_right] = get_NF_bars(log_ERS_ipsi, log_ERS_contra, LOG_ERS_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc);
                
                   
                if feedback == 1

        %             Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                    Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter+xCenter/4 yCenter], 2);
                    Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter-xCenter/4 yCenter], 2);
        %             Screen('FillRect', window, NF_COLOR, NF_bar); 
                    Screen('FillRect', window, NF_color_left, NF_bar_left); 
                    Screen('FillRect', window, NF_color_right, NF_bar_right);

                elseif feedback == 2

                    Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter+xCenter/4 yCenter], 2);
                    Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter-xCenter/4 yCenter], 2);
                    Screen('FillRect', window, NF_color_left, holder_bar_left); 
                    Screen('FillRect', window, NF_color_right, holder_bar_right);

                else 
                    disp('ERROR: feedback condition undefined')
                end
                

                % Flip to the screen
                vbl = Screen('Flip', window, vbl + (WAIT_FRAMES - 0.5) * ifi);

            end
            
            % trial info
            iteration = 0;
            trial_stage = 5;
            
            % keep the feedback on for another 5 seconds
            if feedback == 1
                
                [temp_data, ts] = inlet.pull_chunk();


                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%
                
                
                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration power_rest mean(Pxx(:,1)) mean(Pxx(:,2)) LI log_ERS_ipsi log_ERS_contra];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%
                
                
    %             Screen('DrawLines', window, all_cue_coords, LINE_WIDTH_PIX, FIX_COLOR, [xCenter yCenter], 2);
                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter+xCenter/4 yCenter], 2);
                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter-xCenter/4 yCenter], 2);
    %             Screen('FillRect', window, NF_COLOR, NF_bar); 
                Screen('FillRect', window, NF_color_left, NF_bar_left); 
                Screen('FillRect', window, NF_color_right, NF_bar_right);
                
                % Flip to the screen
                vbl = Screen('Flip', window);

            elseif feedback == 2
                 
                [temp_data, ts] = inlet.pull_chunk();


                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%
                
                
                log_ERS_ipsi = mean(log_ERS_ipsi_list);
                log_ERS_contra = mean(log_ERS_contra_list);
                
               
                %------------ Save Data -----------%
                pwr_mat_temp = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration power_rest mean(Pxx(:,1)) mean(Pxx(:,2)) NaN log_ERS_ipsi log_ERS_contra];
                pwr_mat = [pwr_mat; pwr_mat_temp];
                %----------------------------------%


%                 NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc);
                [NF_bar_left, NF_bar_right] = get_NF_bars(log_ERS_ipsi, log_ERS_contra, LOG_ERS_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc); 

                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter+xCenter/4 yCenter], 2);
                Screen('DrawLines', window, all_NF_coords, LINE_WIDTH_PIX, NF_COLOR, [xCenter-xCenter/4 yCenter], 2);
                Screen('FillRect', window, NF_color_left, NF_bar_left); 
                Screen('FillRect', window, NF_color_right, NF_bar_right);
                
                % Flip to the screen
                vbl = Screen('Flip', window);

            else 
                disp('ERROR: feedback condition undefined')
            end
            
            tic;
            while toc<5
                [temp_data, ts] = inlet.pull_chunk();

                %------------ Save Data -----------%
                iteration = iteration + 1;

                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds block task trial cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%
            end
            %----------------------------------------------------------------------

        end
        
        % collect raw data for an added amount of time to avoid
        % edge-effects with offline filtering
        if block == 4 | firstPress(KbName('ESCAPE'))
            Screen('TextSize', window, 36); 
            DrawFormattedText(window, 'The experiment is over\n\n\nThe experimenter should be with you shortly' ,...
            'center', 'center', white );
            Screen('Flip', window);
            tic;
            while toc < 20
                [temp_data, ts] = inlet.pull_chunk();

                %------------ Save Data -----------%
                nrow = size(temp_data,2);
                raw_info = [id age sex hand year month day hour minute seconds 99 task 99 cue_loc_idx iti trial_stage iteration];
                raw_info_mat = repmat(raw_info, nrow, 1);

                raw_eeg_mat = [ts.' temp_data.'];

                raw_mat_temp = [raw_info_mat raw_eeg_mat];

                raw_mat = [raw_mat; raw_mat_temp];
                %----------------------------------%
            end
        end
        
        % write response matrix to csv
        csvwrite(sprintf('C:/Users/Kine Research/Documents/MATLAB/ghis_data/%s_p%i_raw_%s_block_%i.csv', feedback_str, id, task_str, block), raw_mat);
        csvwrite(sprintf('C:/Users/Kine Research/Documents/MATLAB/ghis_data/%s_p%i_pwr_%s_block_%i.csv', feedback_str, id, task_str, block), pwr_mat);

    end
    
    % turn off screen
    sca;
    
catch
    sca;
    psychrethrow(psychlasterror);
end


