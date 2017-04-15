function [LI, log_ERS_ipsi, log_ERS_contra] = get_LI(cue_loc, Pxx, power_rest_mavg)

    if strcmp(cue_loc, 'right')
        power_contra = mean(Pxx(:,1)); % Ch8 ==> C3
        log_ERS_contra = log(power_contra/power_rest_mavg); % where positive means syncronisation (increased power relative to baseline)

        power_ipsi = mean(Pxx(:,2)); % Ch12 ==> C4
        log_ERS_ipsi = log(power_ipsi/power_rest_mavg);

    elseif strcmp(cue_loc, 'left')
        power_contra = mean(Pxx(:,2)); % Ch12 ==> C4 
        log_ERS_contra = log(power_contra/power_rest_mavg); % where positive means syncronisation (increased power relative to baseline)

        power_ipsi = mean(Pxx(:,1)); % Ch8 ==> C3
        log_ERS_ipsi = log(power_ipsi/power_rest_mavg);

    else
        error('Cue location not properly defined')
    end
    
     LI = (log_ERS_ipsi - log_ERS_contra);

end