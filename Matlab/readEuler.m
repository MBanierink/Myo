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

% Plotting axes
figure(1);
hold on
plot3([min([data(:, 2); 0]) max([data(:, 2); 0])], [0 0],[0 0], '-k'); % for roll-axis
plot3([0 0], [min([data(:, 3); 0]) max([data(:, 3); 0])],[0 0], '-k'); % for pitch-axis
plot3([0 0],[0 0], [min([data(:, 4); 0]) max([data(:, 4); 0])], '-k'); % for yaw-axis
xlabel(headers{2}), ylabel(headers{3}), zlabel(headers{4});
hold off

% Plotting path

for i = 1:length(data(:,1))
    yaw = data(i, 4);
    x1 = cos(yaw);
    y1 = sin(yaw);
    
    pitch = data(i, 3);
    x2 = cos(pitch);
    z1 = sin(pitch);
    
    roll = data(i, 2);
    
    
    hold on
    %plot3([0, data(i, 2)], [0, data(i, 3)], [0, data(i, 4)], '-b');
    plot3(data(i, 2), data(i, 3), data(i, 4), '.b');
    pause(0.05);
    hold off
end
