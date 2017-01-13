function NF_bar = get_NF_bar(LI)

    LI_pix = LI * LI_scale;
        
    if LI_pix > xCenter - line_width_pix/2 
        LI_pix = xCenter - line_width_pix/2;
    elseif LI_pix < -(xCenter - line_width_pix/2)
        LI_pix = -(xCenter - line_width_pix/2);
    end

    % draw right arrow 
    NF_bar_right = [xCenter + line_width_pix/2
        , yCenter - fix_cross_dim_pix/2
        , xCenter + line_width_pix/2 + LI_pix;
        , yCenter + fix_cross_dim_pix/2];

    NF_bar_left = [xCenter - line_width_pix/2 + LI_pix;
        , yCenter - fix_cross_dim_pix/2
        , xCenter - line_width_pix/2 
        , yCenter + fix_cross_dim_pix/2];

    if LI > 0
        NF_bar = NF_bar_right;
    else
        NF_bar = NF_bar_left;
    end
    
end
