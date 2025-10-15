clc; clear; close all;

% Parameters
x = linspace(0, 2*pi, 200);
y = sin(x);
h = plot(x, y, 'LineWidth', 2);
axis([0 2*pi -1.5 1.5]);
grid on;
title('Animated Moving Sine Wave');
xlabel('x');
ylabel('sin(x - t)');

% Create a video writer object
v = VideoWriter('sine_wave_animation.mp4', 'MPEG-4'); % Output file
v.FrameRate = 30;  % Frames per second
open(v);           % Open file for writing

% Animation loop
for t = 0:0.1:10
    y = sin(x - t);
    set(h, 'YData', y);
    drawnow;
    
    % Capture the current frame and write to video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);
disp('Video saved as sine_wave_animation.mp4');
