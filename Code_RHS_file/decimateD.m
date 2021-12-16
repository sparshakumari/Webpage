%-------------------------------------------------------------------------------------------------------
% This code is used to downsample the EEG data
% It is downsampled by a factor of 'factor' - just divide the original sampling frequency by the 'factor'
%---------------------------------------------------------------------------------------------------------

function [down_time, down_data] = decimateD(amplitudes, timestamps, factor)
[x,y] = size(amplitudes); %#ok<ASGLU>
for i = 1:y
    if(i == 1)
    down_data = decimate(amplitudes(:,i), factor,'fir');
    down_time = downsample(timestamps(:, i), factor);
    else
    down_data = [down_data(:, 1:end) decimate(amplitudes(:,i), factor,'fir')];
    down_time = [down_time(:, 1:end) downsample(timestamps(:, i), factor)];
    end
end