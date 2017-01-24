function LI = get_LI(cue_loc, Pxx, power_rest)

    if strcmp(cue_loc, 'right')
        power_contra = mean(Pxx(:,1)); % Ch8 ==> C3
        log2_ERS_contra = log2(power_contra/power_rest); % where positive means syncronisation (increased power relative to baseline)

        power_ipsi = mean(Pxx(:,2)); % Ch12 ==> C4
        log2_ERS_ipsi = log2(power_ipsi/power_rest);

        LI = (log2_ERS_ipsi - log2_ERS_contra);

    elseif strcmp(cue_loc, 'left')
        power_contra = mean(Pxx(:,2)); % Ch12 ==> C4 
        log2_ERS_contra = log2(power_contra/power_rest); % where positive means syncronisation (increased power relative to baseline)

        power_ipsi = mean(Pxx(:,1)); % Ch8 ==> C3
        log2_ERS_ipsi = log2(power_ipsi/power_rest);

        LI =  -(log2_ERS_ipsi - log2_ERS_contra);
    else
        error('Cue location not properly defined')
    end

end