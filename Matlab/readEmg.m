clear all, close all, clc

% Load file
[fileName, path] = uigetfile('*.csv');
file = [path, fileName];

% Read file
data = csvread(file, 1, 0);
fid = fopen(file);
headers = textscan(fid, '%s %s %s %s %s %s %s %s %s', 1, 'delimiter', ',');
fclose(fid);

% Rewriting time
timestamp = data(:, 1) ./ 1000000;
dateTime = datestr(timestamp./(24*60*60) + datenum(1970, 1, 1), 'dd-mm-yyyy HH:MM:SS.FFF');

% Calculate sample rate
timeDiff = (data(end, 1) - data(1, 1)) / 1000000;
sampleRate = length(dateTime) / timeDiff;

% Loop through the columns
for i = 2:numel(headers);
    % Load data
    s = data(:, i);
    Fs = sampleRate;
    t = linspace(0,((length(s))/Fs), length(s));
    
    % Show signal
    figure(i-1), subplot(2, 1, 1);
    plot(t, s), title('EMG'), xlabel('t (s)'), ylabel('Amplitude');
    
    % Fourier analysis
    f = 0:1/max(t):Fs/2;
    y = abs(fft(s));
    n = length(f);
    y = y(1:n);
    
    % Show frequency spectrum
    figure(i-1), subplot(2, 1, 2);
    plot(f, y), title('FFT'), xlabel('Frequency (Hz)'), ylabel('Amplitude');
end





