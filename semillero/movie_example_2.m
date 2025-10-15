clc; clear; close all;

% Parameters
x = linspace(0, 2*pi, 200);
y = sin(x);

% Create figure with fixed size
fig = figure('Position', [100 100 600 400]); % fixed window size
h = plot(x, y, 'LineWidth', 2);
axis([0 2*pi -1.5 1.5]);     % fix axis limits
axis manual;                 % prevents automatic rescaling
grid on;
title('Animated Moving Sine Wave');
xlabel('x');
ylabel('sin(x - t)');

% Prepare video writer
v = VideoWriter('sine_wave_stable.mp4', 'MPEG-4');
v.FrameRate = 30;
v.Quality = 100;
open(v);

% Animation loop
for t = 0:0.1:10
    y = sin(x - t);
    set(h, 'YData', y);
    drawnow;
    
    % Capture consistent frame size
    frame = getframe(fig);   % use fixed figure handle
    writeVideo(v, frame);
end

% Close video
close(v);
disp('Video saved as sine_wave_stable.mp4');
