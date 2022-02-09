function []=plot_scope_spectrum_analyzer(fignum,t,y,zoomt, zoomf, zoomp)
% P. Axelrad 
% Generates plots in the time domain and frequency domain of a signal
% Both are generated at the full scale and at a zoomed in scale that you
% choose, for a total of 4 subplots.
% 
% Inputs 
% - fignum (which figure number to use)
% t - time vector (assumed to be evenly spaced at the sample time
% y - signal (amplitude) that you are analyzing
% zoomt - is max value of time for the zoomed in plot
%
% zoomf - [min max] value of frequency for the zoomed in plot
% zoomp - [min max] values of power for plots
% 
% Example call: 
% t = 0:.001:10;
% y = sin(2*pi*3*t); 
% plot_scope_spectrum_analyzer(2,t,y,2,[0 5],[-100 0])
% 
npts = length(t);
fs = 1/(t(2)-t(1));
f = fs*(0:npts-1)/npts;

% Time Plots - Signal Amplitude vs Time
% Plots on the left are at full scale, right is zoomed in
figure(fignum)
subplot(2,2,1),plot(t,y),xlabel('TIME (SEC)'),ylabel('AMPLITUDE (V)')
ylim([-1.1*max(y) 1.1*max(y)]), grid
subplot(2,2,2),plot(t,y),xlabel('TIME (SEC)')
ylim([-1.1*max(y) 1.1*max(y)]), grid
if nargin>3 
    xlim([0 zoomt])
end

py = 20*log10(sqrt(2)*abs(fft(y)/npts)); 
disp(['Question ' num2str(fignum) ', max P: ' num2str(max(py))])

subplot(2,2,3),plot(f,py,'b-',f,py,'b.'),xlabel('FREQ (Hz)'),ylabel('POWER (dBW)')
xlim([0 max(f)/2]), grid
if nargin>5
    % ylim(zoomp)
end
subplot(2,2,4),plot(f,py,'b-',f,py,'b.'),xlabel('FREQ (Hz)'),
if nargin>4
    xlim(zoomf)
end
if nargin>5
    ylim(zoomp)
end
grid

end
