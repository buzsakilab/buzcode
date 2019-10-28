function [t2, v2] = upsample2x(t, v);

% [t2, v2] = upsample2x(t, v);
%
% Upsample an array of sampled voltages by a factor of 2.  Cublic spline
% interpolation is used to create the new samples.  The function returns
% the interpolated data along with a new time vector.
%
% Example:
% >> [t_amplifier_2x, amplifier_data_2x] = upsample2x(t_amplifier, amplifier_data);

tic;
fprintf(1, '\n');
fprintf(1, 'Upsampling waveforms by 2X...\n');

sample_rate = 1 / (t(2) - t(1));
[num_channels, num_samples] = size(v);

t2 = [t ; (t + 1/(2 * sample_rate))];
t2 = reshape(t2, 1, 2 * num_samples);

v2 = zeros(num_channels, 2 * num_samples);

print_increment = 10;
percent_done = print_increment;
for i=1:num_channels
    v2(i,:) = interp1(t, v(i,:), t2, 'spline');

    fraction_done = 100 * (i / num_channels);
    if (fraction_done >= percent_done)
        fprintf(1, '%d%% done...\n', percent_done);
        percent_done = percent_done + print_increment;
    end
end

fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', toc);
fprintf(1, '\n');

return
