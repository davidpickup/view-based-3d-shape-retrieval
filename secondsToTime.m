function [hours, minutes, seconds] = secondsToTime(seconds)
% [hours, minutes, seconds] = secondsToTime(seconds)
%
% David Pickup 2015

hours = floor((seconds / 60) / 60);
seconds = seconds - (hours*60*60);
minutes = floor(seconds/60);
seconds = round(seconds - (minutes*60));

return;