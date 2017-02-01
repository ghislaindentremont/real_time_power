function [Pxx, Fxx, data_buffer] = get_power(temp_data, data_points, data_buffer, pad_points, CHANNELS_OF_INTEREST, PSD_FREQS, FS, a, b, a2, b2)

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

    [Pxx, Fxx] = pwelch(display_buffer, [], [], PSD_FREQS, FS, 'power');

end

