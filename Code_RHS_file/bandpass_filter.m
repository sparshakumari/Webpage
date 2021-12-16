%-----------------------------------------------------------------------------------------
% Purpose:
% Function to filter the EEG signal using bandpass filter.
% It is the most preferred form of filter for EEG data.
%
% Inputs:
% inputData - raw data           ---- type : vector
% fs        - sampling frequency ---- type : scalar (downsampled frequency)
%
% Outputs:
% filteredData - output ---- type : vector
% 
% Function used: +

% bandpass(input_data, [lowpass_value highpass_value], sampling frequency)
% 
% Eg:
% input_data = [1, 2, 3, 4, 5, 6]
% lowpass_value = 0.1
% highpass_value = 100
% fs = 1000          
%
% Author: Sparsha Kumari, Manipal School of Life Sciences
%
% Started: 01-01-2020
%
% Last update: 20-02-2020
%-----------------------------------------------------------------------------------------

function [filteredData] = bandpass_filter(inputData, l, h, fs)
filteredData = bandpass(inputData, [l h], fs);
end