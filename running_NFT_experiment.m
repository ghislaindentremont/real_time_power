% Clear the workspace and the screen
sca;
close all;
clearvars;

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
HideCursor(); 

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
fix_cross_dim_pix = 40*2;

% Now we set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
x_fix_coords = [-fix_cross_dim_pix fix_cross_dim_pix 0 0];
y_fix_coords = [0 0 -fix_cross_dim_pix fix_cross_dim_pix];
all_fix_coords = [x_fix_coords; y_fix_coords];

% Set the line width for our fixation cross
line_width_pix = 4;

% make it green 
fix_color = [0 1 0];

% just vertical line
all_cue_coords = [[0 0 0 0]; y_fix_coords];



%----------------------------------------------------------------------
%                         Arrow Cue 
%----------------------------------------------------------------------

% arrow color - red 
arr_color = [1 0 0];

% draw right arrow 
rarrow_base = [xCenter + line_width_pix/2
    , yCenter - fix_cross_dim_pix/2
    , xCenter + line_width_pix/2 + fix_cross_dim_pix
    , yCenter + fix_cross_dim_pix/2];

rarrow_spear = [[xCenter + line_width_pix/2 + fix_cross_dim_pix, yCenter + fix_cross_dim_pix];
    [xCenter + line_width_pix/2 + fix_cross_dim_pix*2, yCenter];
    [xCenter + line_width_pix/2 + fix_cross_dim_pix, yCenter - fix_cross_dim_pix]];

% draw right arrow 
larrow_base = [xCenter - line_width_pix/2 - fix_cross_dim_pix
    , yCenter - fix_cross_dim_pix/2
    , xCenter - line_width_pix/2 
    , yCenter + fix_cross_dim_pix/2];

larrow_spear = [[xCenter - line_width_pix/2 - fix_cross_dim_pix, yCenter + fix_cross_dim_pix];
    [xCenter - line_width_pix/2 - fix_cross_dim_pix*2, yCenter];
    [xCenter - line_width_pix/2 - fix_cross_dim_pix, yCenter - fix_cross_dim_pix]];



%----------------------------------------------------------------------
%                         NF Bar 
%----------------------------------------------------------------------

NF_color = [0 0 1];



%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------

% Baseline duration
baseline_time = 2;

% Fixation interval time in seconds and frames
fix_time = 2;
fix_time_frames = round(fix_time / ifi);

% Cue interval time in seconds and frames 
cue_time = 1.25;
cue_time_frames = round(cue_time / ifi);

% Numer of frames to wait before re-drawing
wait_frames = 1;

% NF time frames 
NF_time = 5;
NF_wait_frames = 3;
NF_time_frames = round(NF_time / ifi / NF_wait_frames);




%----------------------------------------------------------------------
%                       Conditions
%----------------------------------------------------------------------

% cue location
cue_locs_list = {'left', 'right'};
cue_locs = [1, 2];

trials_per_condition = 2;
cond_matrix = repmat(cue_locs, 1, trials_per_condition);

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
%                       Experimental Loop
%----------------------------------------------------------------------


for trial = 1:num_trials
    
    %---------------- Baseline EEG Aquisition -------------------------
    if trial == 1
        Screen('TextSize', window, 60); 
        DrawFormattedText(window, 'Press Any Key To Begin The Experiment',...
        'center', 'center', white );
        Screen('Flip', window);
        KbStrokeWait; 
        
        % baseline EEG aquisition 
        Screen('TextSize', window, 60); 
        DrawFormattedText(window, 'Please Relax With Your Eyes Open',...
        'center', 'center', white );
        Screen('Flip', window);
        
        tic;
        while toc < baseline_time end 
        
    end
    %----------------------------------------------------------------------

    
%     %-------------------Trial Initiation message --------------------------
%     Screen('TextSize', window, 60); 
%     DrawFormattedText(window, 'Press Any Key To Begin the Trial',...
%         'center', 'center', white );
%     Screen('Flip', window);
%     KbStrokeWait; 
%     %----------------------------------------------------------------------


    %----------------------------ITI --------------------------------------
    Screen('FillRect', window, black, window_rect );
    Screen('Flip', window);
    
    iti = cond_matrix_shuffled(2,trial);
    tic;
    while toc < iti end
    %----------------------------------------------------------------------


    %--------------------- Draw the fixation cross ------------------------
    Screen('DrawLines', window, all_fix_coords, line_width_pix, fix_color, [xCenter yCenter], 2);
    vbl = Screen('Flip', window); 

    % Flip to the screen
    for frame = 1:fix_time_frames-1
        % Draw the fixation point
        Screen('DrawLines', window, all_fix_coords, line_width_pix, fix_color, [xCenter yCenter], 2);

        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (wait_frames - 0.5) * ifi);
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
    Screen('DrawLines', window, all_cue_coords, line_width_pix, fix_color, [xCenter yCenter], 2);
    Screen('FillRect', window, arr_color, arrow_base);
    Screen('FillPoly', window, arr_color, arrow_spear, 1);

    % Flip to the screen
    vbl = Screen('Flip', window);

    % Flip to the screen
    for frame = 1:cue_time_frames-1
        Screen('DrawLines', window, all_cue_coords, line_width_pix, fix_color, [xCenter yCenter], 2);
        Screen('FillRect', window, arr_color, arrow_base);
        Screen('FillPoly', window, arr_color, arrow_spear, 1);

        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (wait_frames - 0.5) * ifi);
    end
    %----------------------------------------------------------------------

    
    %--------------------- Draw NF bar -----------------------------------    
    log_power_ratio_in_pixels = round(randn * fix_cross_dim_pix + 1);
    
    % draw right arrow 
    NF_bar_right = [xCenter + line_width_pix/2
        , yCenter - fix_cross_dim_pix/2
        , xCenter + line_width_pix/2 + log_power_ratio_in_pixels
        , yCenter + fix_cross_dim_pix/2];

    NF_bar_left = [xCenter - line_width_pix/2 + log_power_ratio_in_pixels
        , yCenter - fix_cross_dim_pix/2
        , xCenter - line_width_pix/2 
        , yCenter + fix_cross_dim_pix/2];
    
    if log_power_ratio_in_pixels == 0
        log_power_ratio_in_pixels = 1;
    end 

    if log_power_ratio_in_pixels > 0
        NF_bar = NF_bar_right;
    else
        NF_bar = NF_bar_left;
    end
    
    % Draw NF bar 
    Screen('DrawLines', window, all_cue_coords, line_width_pix, fix_color, [xCenter yCenter], 2);
    Screen('FillRect', window, NF_color, NF_bar);

    % Flip to the screen
    vbl = Screen('Flip', window);

    % Flip to the screen
    for frame = 1:NF_time_frames-1
        
        change_pixels = randn;
        
        log_power_ratio_in_pixels = round(log_power_ratio_in_pixels + change_pixels);
        
        if log_power_ratio_in_pixels == 0
            log_power_ratio_in_pixels = 1;
        end
    
        if log_power_ratio_in_pixels > 0
            NF_bar = [xCenter + line_width_pix/2
            , yCenter - fix_cross_dim_pix/2
            , xCenter + line_width_pix/2 + log_power_ratio_in_pixels
            , yCenter + fix_cross_dim_pix/2];
        else
            NF_bar = [xCenter - line_width_pix/2 + log_power_ratio_in_pixels
            , yCenter - fix_cross_dim_pix/2
            , xCenter - line_width_pix/2 
            , yCenter + fix_cross_dim_pix/2];
        end
        
        % Draw NF bar
        Screen('DrawLines', window, all_cue_coords, line_width_pix, fix_color, [xCenter yCenter], 2);
        Screen('FillRect', window, NF_color, NF_bar);
       
        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (NF_wait_frames - 0.5) * ifi);

    end
    %----------------------------------------------------------------------
    
    
    %------------------- End of Experiment --------------------------------
    if trial == 4
        
        Screen('TextSize', window, 60); 
        DrawFormattedText(window, 'You Have Succesfully Completed The Experiment' ,...
        'center', 'center', white );
        Screen('Flip', window);
        tic;
        while toc < 3 end
         
        Screen('TextSize', window, 60); 
        DrawFormattedText(window, 'The Experimenter Should Be With You Shortly',...
        'center', 'center', white );
        Screen('Flip', window);
        KbStrokeWait; 
        
    end
    %----------------------------------------------------------------------
    
end
    
% Clear the screen
sca;


