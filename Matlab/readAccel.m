clear all, close all, clc

% Load file
[fileName, path] = uigetfile('*.csv');
file = [path, fileName];

% Read file
data = csvread(file, 1, 0);
fid = fopen(file);
headers = textscan(fid, '%s %s %s %s', 1, 'delimiter', ',');
fclose(fid);

% Rewriting time
timestamp = data(:, 1) ./ 1000000;
dateTime = datestr(timestamp./(24*60*60) + datenum(1970, 1, 1), 'dd-mm-yyyy HH:MM:SS.FFF');

% Calculate sample rate
timeDiff = (data(end, 1) - data(1, 1)) / 1000000;
sampleRate = length(dateTime) / timeDiff;

% Time per acceleration
t = 1/sampleRate;

% Calculate gravitational acceleration
GRAVITY = 9.81; % m/s^2
data(:, 2:4) = data(:, 2:4) * GRAVITY;

% Compensate for drift or bad IMU orientation
% DRIFT_SAMPLE_TIME = 1; % Seconds
% xDrift = mean(data(1:int64(sampleRate) * DRIFT_SAMPLE_TIME, 2));
% yDrift = mean(data(1:int64(sampleRate) * DRIFT_SAMPLE_TIME, 3));
% zDrift = mean(data(1:int64(sampleRate) * DRIFT_SAMPLE_TIME, 4));

% Compensate for gravity
xData = data(:, 2);% - xDrift;
yData = data(:, 3);% - yDrift;
zData = data(:, 4);% - zDrift;

for k = 2:4
    s = data(:, k);
    Fs = sampleRate;
    t = linspace(0,((length(s))/Fs), length(s));
    
    % Show signal
    figure(20 + k), subplot(3, 1, 1);
    plot(t, s), title('Accelerometer'), xlabel('t (s)'), ylabel('Amplitude');
    
    % Fourier analysis
    f = 0:1/max(t):Fs/2;
    y = abs(fft(s));
    n = length(f);
    y = y(1:n);
    
    % Show frequency spectrum
    figure(20 + k), subplot(3, 1, 2);
    plot(f, y), title('FFT'), xlabel('Frequency (Hz)'), ylabel('Amplitude');
    
    % Create simple filter
    BAND_A = 1;
    BAND_B = 50;
    filter = zeros(length(y), 1);
    filter(BAND_A:BAND_B) = 1;
    hold on
    plot(f, filter * max(y), '-r');
    plot(f, filter .* y, '-g');
    hold off
    
    % Filter and inverse Fourier
    y2 = y .* filter;
    s2 = ifft(y2);
    figure(20 + k), subplot(3, 1, 3);
    plot(f, s2), title('Accelerometer - filtered'), xlabel('t (s)'), ylabel('Amplitude');
end

% Filtering data
[b, a] = butter(10, 0.9, 'high');
xAccHigh = filtfilt(b, a, xData);
xAccHighFilt = sqrt(xAccHigh.^2);
[b2, a2] = butter(1, 0.1, 'low');
xAccLow = filtfilt(b2, a2, xData);
xAcc = xAccHighFilt .* xAccLow;
yAcc = filtfilt(b, a, yData);
zAcc = filtfilt(b, a, zData);

% Time per acceleration
t = 1/sampleRate;

% Calculate velocity
xVel = cumsum(xAcc) * t;
yVel = cumsum(yAcc) * t;
zVel = cumsum(zAcc) * t;

% Calculate position
xPos = cumsum(xVel) * t;
yPos = cumsum(yVel) * t;
zPos = cumsum(zVel) * t;

% Plotting data seperately
figure(1);
subplot(6, 1, 1), plot(xData);
subplot(6, 1, 2), plot(xAccHighFilt);
subplot(6, 1, 3), plot(xAccLow);
subplot(6, 1, 4), plot(xAcc);
subplot(6, 1, 5), plot(xVel);
subplot(6, 1, 6), plot(xPos);

figure(2);
subplot(4, 1, 1), plot(yData);
subplot(4, 1, 2), plot(yAcc);
subplot(4, 1, 3), plot(yVel);
subplot(4, 1, 4), plot(yPos);

figure(3);
subplot(4, 1, 1), plot(zData);
subplot(4, 1, 2), plot(zAcc);
subplot(4, 1, 3), plot(zVel);
subplot(4, 1, 4), plot(zPos);



% Plotting axes
% xmin = abs(min(xpos));
% xmax = max(xpos);
% ymin = abs(min(ypos));
% ymax = max(ypos);
% zmin = abs(min(zpos));
% zmax = max(zpos);
% 
% axSize = max([xmin, xmax, xmin, ymax, zmin, zmax]);

% hold on
% plot3([-axSize axSize], [0 0], [0 0], '-k'); % for x-axis
% plot3([0 0], [-axSize axSize], [0 0], '-k'); % for y-axis
% plot3([0 0], [0 0], [-axSize axSize], '-k'); % for z-axis
% hold off
% xlabel(headers{2}), ylabel(headers{3}), zlabel(headers{4});

% Estimate drift
% DRIFT_SAMPLE_TIME = 1; % seconds
% xDrift = mean(diff(xacc(end - DRIFT_SAMPLE_TIME * int64(sampleRate) : end)));
% yDrift = mean(diff(yacc(end - DRIFT_SAMPLE_TIME * int64(sampleRate) : end)));
% zDrift = mean(diff(zacc(end - DRIFT_SAMPLE_TIME * int64(sampleRate) : end)));

% Plotting path
for i = 1:length(data(:,1))
    
    % Calculating positions
    x = xPos(i);
    y = yPos(i);
    z = zPos(i);

    figure(4);
    hold on
    plot3(x, y, z, '.b');
    %pause(0.02);
    hold off
    
end
