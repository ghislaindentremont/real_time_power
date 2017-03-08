function NF_bar = get_NF_bar(LI, LI_SCALE, xCenter, yCenter, LINE_WIDTH_PIX, FIX_CROSS_DIM_PIX, cue_loc)
    
    if strcmp(cue_loc, 'right')
        LI_pix = LI * LI_SCALE;
    elseif strcmp(cue_loc, 'left')
        LI_pix = - LI * LI_SCALE;
    else
        error('Cue location not properly defined')
    end
        
    if LI_pix > xCenter - LINE_WIDTH_PIX/2 
        LI_pix = xCenter - LINE_WIDTH_PIX/2;
    elseif LI_pix < -(xCenter - LINE_WIDTH_PIX/2)
        LI_pix = -(xCenter - LINE_WIDTH_PIX/2);
    end

    % draw right arrow 
    NF_bar_right = [xCenter + LINE_WIDTH_PIX/2
        , yCenter - FIX_CROSS_DIM_PIX/2
        , xCenter + LINE_WIDTH_PIX/2 + LI_pix;
        , yCenter + FIX_CROSS_DIM_PIX/2];

    NF_bar_left = [xCenter - LINE_WIDTH_PIX/2 + LI_pix;
        , yCenter - FIX_CROSS_DIM_PIX/2
        , xCenter - LINE_WIDTH_PIX/2 
        , yCenter + FIX_CROSS_DIM_PIX/2];

    if LI_pix > 0
        NF_bar = NF_bar_right;
    else
        NF_bar = NF_bar_left;
    end
    
end
