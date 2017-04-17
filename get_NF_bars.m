function [NF_bar_left, NF_bar_right] = get_NF_bars(log_ERS_ipsi, log_ERS_contra, LOG_ERS_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc)
        
    % since y-coordinates are flipped such that zero is at the top of the
    % screen, we treat ERS as positive so that when flipped, ERS is
    % represented as a downward bar whereas ERD is represented as a upward
    % bar
    contra_pix = log_ERS_contra * LOG_ERS_SCALE;
    ipsi_pix = log_ERS_ipsi * LOG_ERS_SCALE;

    if contra_pix > yCenter - LINE_WIDTH_PIX/2 
        contra_pix = yCenter - LINE_WIDTH_PIX/2;
    elseif contra_pix < -(yCenter - LINE_WIDTH_PIX/2) 
        contra_pix = -(yCenter - LINE_WIDTH_PIX/2);
    end

    if ipsi_pix > yCenter - LINE_WIDTH_PIX/2 
        ipsi_pix = yCenter - LINE_WIDTH_PIX/2;
    elseif ipsi_pix < -(yCenter - LINE_WIDTH_PIX/2) 
        ipsi_pix = -(yCenter - LINE_WIDTH_PIX/2);
    end
    
    if strcmp(cue_loc, 'right')
        
        if contra_pix < 0
        
            NF_bar_right = [xCenter + xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter - LINE_WIDTH_PIX/2 + contra_pix
                , xCenter + xCenter/4 + FIX_CROSS_DIM_PIX/2
                , yCenter - LINE_WIDTH_PIX/2];
            
        else
            
            NF_bar_right = [xCenter + xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter + LINE_WIDTH_PIX/2 
                , xCenter + xCenter/4 + FIX_CROSS_DIM_PIX/2
                , yCenter +  LINE_WIDTH_PIX/2 + contra_pix];
            
        end
        
        if ipsi_pix < 0 
            
            NF_bar_left = [xCenter - xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter - LINE_WIDTH_PIX/2 + ipsi_pix
                , xCenter - xCenter/4 + FIX_CROSS_DIM_PIX/2;
                , yCenter - LINE_WIDTH_PIX/2];
            
        else
            
            NF_bar_left = [xCenter - xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter + LINE_WIDTH_PIX/2 
                , xCenter - xCenter/4 + FIX_CROSS_DIM_PIX/2;
                , yCenter + LINE_WIDTH_PIX/2 + ipsi_pix];
            
        end



        
    elseif strcmp(cue_loc, 'left')
        
        if ipsi_pix < 0
        
            NF_bar_right = [xCenter + xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter - LINE_WIDTH_PIX/2 + ipsi_pix
                , xCenter + xCenter/4 + FIX_CROSS_DIM_PIX/2
                , yCenter - LINE_WIDTH_PIX/2];
            
        else
            
            NF_bar_right = [xCenter + xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter + LINE_WIDTH_PIX/2 
                , xCenter + xCenter/4 + FIX_CROSS_DIM_PIX/2
                , yCenter +  LINE_WIDTH_PIX/2 + ipsi_pix];
            
        end
        
        if contra_pix < 0 
            
            NF_bar_left = [xCenter - xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter - LINE_WIDTH_PIX/2 + contra_pix
                , xCenter - xCenter/4 + FIX_CROSS_DIM_PIX/2;
                , yCenter - LINE_WIDTH_PIX/2];
            
        else
            
            NF_bar_left = [xCenter - xCenter/4 - FIX_CROSS_DIM_PIX/2
                , yCenter + LINE_WIDTH_PIX/2 
                , xCenter - xCenter/4 + FIX_CROSS_DIM_PIX/2;
                , yCenter + LINE_WIDTH_PIX/2 + contra_pix];
            
        end
        
    end        
    
end
