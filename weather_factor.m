clear all;clc;close all;

% Constants for weather factor calculation
MAX_TEMP = 30; % 
MIN_TEMP = 0; % 
MAX_PRECIP = 50; % millimeters
MIN_PRECIP = 0; % millimeters
SHADING =1; % percentage of shading (between 0 (not cloudy and 1 is cloudy)

% Temperature, precipitation, and shading data for winter
temperature_winter = [11.4, 10.9, 10.8, 10.8, 11, 10.8, 10.8, 10.7, 10.7, ...
10.2, 11.1, 11.1, 11.2, 11.3, 11.3, 11.1, 10.3, 9.8, 9.5, 8.9, 8.4, ...
8.2, 8, 7.8]; % degrees Celsius
precipitation_winter = [4.733, 0.263, 0, 0, 0, 0, 9.378, 0.061, 0, 0, ...
0, 0, 0, 0, 0.106, 0.071, 0.167, 0, 0.155, 0.666, 0.848, 0.409, ...
0.391, 0.406]; % millimeters

% Calculate winter weather factor
temp_factor_winter = (MAX_TEMP - temperature_winter) / (MAX_TEMP - MIN_TEMP);
shading_factor_winter = 1 - SHADING;
weather_factor_winter_cloudy = (temp_factor_winter * shading_factor_winter + ...
(MAX_PRECIP - precipitation_winter) / (MAX_PRECIP - MIN_PRECIP)) / 2;

% Temperature, precipitation, and shading data for summer
temperature_summer = [13.3, 12.4, 11.7, 11, 10.4, 10.8, 10.5, 12, 14.7, ...
18.6, 20.8, 21.5, 21.8, 21.9, 22.4, 22.8, 22.9, 22.6, 22.4, 22.3, ...
20.7, 18.5, 16.5, 14.8]; % degrees Celsius
precipitation_summer = zeros(size(temperature_summer)); % millimeters

% Calculate summer weather factor
temp_factor_summer = (MAX_TEMP - temperature_summer) / (MAX_TEMP - MIN_TEMP);
shading_factor_summer = 1 - SHADING;
weather_factor_summer_cloudy = (temp_factor_summer * shading_factor_summer + ...
(MAX_PRECIP - precipitation_summer) / (MAX_PRECIP - MIN_PRECIP)) / 2;
% save ('WEATHER_FACTOR.mat',  'weather_factor_summer_cloudy',  'weather_factor_winter_cloudy','-append')
