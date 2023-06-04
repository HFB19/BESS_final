 clear all
load import_tariff.mat
load("basepower.mat")


%%%%%%%%%%%%%%%%%
% Input date and usage probability  here(0 to 1)%
%%%%%%%%%%%%%%%%%%
D =' 26-June-2022';%%only choose winter
usage_probability=0.5;
[DayNumber,DayName] = weekday(D);

if DayNumber == 1 || DayNumber == 7
    weekend = 1;
    weekday = 0;
    disp('Your day is a weekend')
else
    weekend = 0;
    weekday = 1;
    disp('Your day is a weekday')
end

% Convert '1-Jan-2021' to a MATLAB serial date number
d1 = datenum('1-Jan-2022');

% Convert D to a MATLAB serial date number
d2 = datenum(D);

% Calculate the difference in days
NumDays= abs(d2 - d1);


if (NumDays >= 1 && NumDays <= 59) || (NumDays>= 244 && NumDays<= 365)
    season = 'winter';
    disp('You chose a winter day!') 
elseif NumDays >= 60 && NumDays <= 243
    season = 'summer';
    disp('You chose a summer day!') 
else
    season = '';
    error('Invalid date. Please choose another day.');
end 

% Determine the appropriate model based on weekday/weekend and season
if weekday && strcmp(season, 'winter')
    day_season = 11; % Winter weekday model
    import_tariff=import_tariff_weekday_winter;
    base=basepower_weekday_winter*1000;
elseif weekday && strcmp(season, 'summer')
    day_season = 12; % Summer weekday model
    import_tariff=import_tariff_weekday_summer;
    base=basepower_weekday_summer*1000;
elseif weekend && strcmp(season, 'winter')
    day_season = 21; % Winter weekend mode
    import_tariff=import_tariff_weekend_winter;
    base=basepower_weekend_winter*1000;

elseif weekend && strcmp(season, 'summer')
    day_season = 22; % Summer weekend model
     import_tariff=import_tariff_weekend_summer;
     base=basepower_weekend_summer*1000;

else
    day_season = 0; % Invalid model
    error('Invalid date. Please choose another day.');
end

%======================================================================
%                           Categorizing Tariffs
%======================================================================
% Description:
% This section categorizes tariffs based on their pricing structure and
% usage patterns.


hour=transpose(1:48);
for i = 1:48
   half_hour_average(i) =mean(import_tariff(i,1));
end

% Sort the data into quartiles
quartiles = prctile(half_hour_average, [0 25 50 75 100]);

% Classify the half_hour_average into 4 import_categories based on quartiles
import_categories = zeros(size(half_hour_average));
import_categories(half_hour_average <= quartiles(2)) = 1; % lower quartile
import_categories(half_hour_average > quartiles(2) & half_hour_average <= quartiles(3)) = 2; % middle quartiles
import_categories(half_hour_average > quartiles(3) & half_hour_average <= quartiles(4)) = 3;
import_categories(half_hour_average > quartiles(4)) = 4; % upper quartile
import_categories=import_categories';




% Get the list of files in the probability directory
files = dir('probability/*.mat');

% Loop through each file and load its contents
for i = 1:length(files)
    load(fullfile(files(i).folder, files(i).name));
end
% initialize demand
no_demand =zeros(48,1);
demand=zeros(48,1);


t_1440 = linspace(1, 24, 48)';
cook_consumption=zeros(48,1);
tv_consumption=zeros(48,1);
%%activate all threshold to be zero everyday at the start
threshold_reached = false;


for y = 1:48
    % Fridge and heater consumption is made from a different file based on researched data
    demand(y) = no_demand(y) + fridge_consumption(y);


    % Additional consumption sources: %https://data.open-power-system-data.org/household_data/
    demand(y) = demand(y) + 10; % Modern
    demand(y) = demand(y) + 366.28; % Freezer


    if TV_Probability(y) > usage_probability
        demand(y) = demand(y) + 100;
        tv_consumption(y) = 100;
    end
    
    if Cook_Probability(y) > usage_probability
        demand(y) = demand(y) + 2000;
        cook_consumption(y) = 2000;
        threshold_reached = true;
        
        if threshold_reached
            Cook_Probability(y+1:end) = 0; % Set the probability to zero for all further y values
        else
            Cook_Probability(y) = Cook_Probability(y);
        end
    end

%======================================================================
%                           Cooling & Heating appliances
%                           (temperature dependant)
%======================================================================

    if season=='winter'
        demand(y)=demand(y)+heater(y);
    elseif season=='summer'
      air_cond=repelem(ac_power_consumption', 2)';
      demand(y)=demand(y)+air_cond(y)*1000;
    end

%======================================================================
%                           Iron and Vacuum appliances
%                           (only done during weekends)
%======================================================================
    if weekend==1 
        if iron_probability(y) > usage_probability
            demand(y) = demand(y) + 1100;
            threshold_reached = true;
                if threshold_reached
                    iron_probability(y+1:end) = 0; % Set the probability to zero for all further y values
                else
                   iron_probability(y) = iron_probability(y);
                end
        end
    
            if vacuum_probability(y) > usage_probability
            demand(y) = demand(y) + 1400;
            threshold_reached = true;
                if threshold_reached
                    vacuum_probability(y+1:end) = 0; % Set the probability to zero for all further y values
                else
                   vacuum_probability(y) = vacuum_probability(y);
                end
            end
    end


end

%======================================================================
%                           Plotting long duration appliances
%======================================================================
t_48=linspace(1,24,48);
plot(t_48,demand, 'LineWidth',2)
xlim([1 24])
title('Huge appliances demand', 'FontSize', 20)
ylabel('Power(W)', 'FontSize', 14); xlabel('Day Hour(hr)', 'FontSize', 20)



%======================================================================
%                           Plotting short duration appliances
%======================================================================
minute_demand = demand;


% Define a new time vector with 30 points between each data point
t_48 = linspace(1, 48,48);
t_new= linspace(1, 48,1440);

% Interpolate the data using the new time vector
long_durationappliances_demand = transpose(interp1(t_48, minute_demand, t_new));

% Set the total number of 30-minute intervals in 24 hours
n_intervals = 48; % 24 hours x 2 (30-minute intervals per hour)

% Create an array to store the power consumption for each interval
kettle_consumption = zeros(n_intervals, 1);
hairdryer_consumption = zeros(n_intervals, 1);
waterheater_consumption = zeros(n_intervals, 1);
microwave_consumption=zeros(n_intervals, 1);
% Simulate the power consumption for each 30-minute interval
for i = 1:n_intervals


if (import_categories(i) == 1)
   if (laundry_probability(i) > usage_probability)
        laundry_probability(i+1:end) = 0; % set all future laundry values to 0
        laundry_pattern(1:225) = average_laundry_data; % set pattern to average laundry data
        idx = (i-1)*30 + 1;% starting index for this hour of data
        laundry_consumption=zeros(1440,1);
        laundry_consumption(idx:idx+224) = laundry_pattern'; % write the first half hour of pattern
        dishwasher_probability(i+1:end) = 0; % set all future  dishwasher values to 0
        dishwasher_pattern(1:77) = average_dishwasher_data; % set pattern to average laundry data
        dishwasher_consumption=zeros(1440,1);
        dishwasher_consumption(idx:idx+76) =  dishwasher_pattern'; 
    end
    end
    
    
    
%======================================================================
%                           Kettle pattern
%======================================================================
        % Determine if the kettle is on for this interval
        p1 = 2000; % watts (power consumption during the kettle-on period)
        p2 = 0; % watts (power consumption during the kettle-off period)
        if kettle(i) > 0.6% probability threshold
            kettle_pattern(1:5) = 2e03; % kettle on for 2 minutes
        else
            kettle_pattern(1:30) = p2; % kettle off for the entire interval
        end
        % Store the 30-minute kettle consumption pattern in the kettle_consumption array
        idx = (i-1)*30 + 1; % starting index for this 30-minute interval
        kettle_consumption(idx:idx+29) = kettle_pattern;
%======================================================================
%                           Hairdryer pattern
%======================================================================
        p1 = 2000; % watts (power consumption during the hairdryer-on period)
        p2 = 0; % watts (power consumption during the hairdryer-off period)
        hairdryer_pattern = zeros(30,1);
        if hairdryer(i) > usage_probability
            hairdryer_pattern(1:5) = p1; % hairdryer on for 5 minutes
            hairdryer(i+1:end) = 0; % set all future hairdryervalues to 0
        else
            hairdryer_pattern(1:30) = p2;   
        end
        
        % Store the 30-minute hairdryer consumption pattern in the hairdryerconsumption array
        idx = (i-1)*30 + 1; % starting index for this 30-minute interval
        hairdryer_consumption(idx:idx+29) = hairdryer_pattern;
%======================================================================
%                           Water Heater pattern
%======================================================================

        p1 = 4e03; % watts (power consumption during the waterheater-on period)
        p2 = 0; % watts (power consumption during the waterheater-off period)
        waterheater_pattern = zeros(30,1);
        if hairdryer(i) > usage_probability
            waterheater_pattern(1:10) = p1;% waterheater on for 10 minutes
            waterheater(i+1:end) = 0; % set all future hairdryer values to 0
        else
            waterheater_pattern(1:30) = p2;   
        end
        
        % Store the 30-minute waterheater consumption pattern in the waterheater_consumption array
        idx = (i-1)*30 + 1; % starting index for this 30-minute interval
        waterheater_consumption(idx:idx+29) = waterheater_pattern;

%======================================================================
%                           Microwave pattern
%======================================================================

        p1 = 1e03; % watts (power consumption during the microwave-on period)
        p2 = 0; % watts (power consumption during the microwave-off period)
        microwave_pattern = zeros(30,1);
        if microwave(i) > usage_probability 
            microwave_pattern(1:5) = p1; % microwave on for 2 minutes
            microwave_pattern(6:30)= p2;
    
        else
             microwave_pattern(1:30) = p2; % microwave off for the entire interval
        end
        % Store the 30-minute microwave consumption pattern in the microwave_consumption array
        idx = (i-1)*30 + 1; % starting index for this 30-minute interval
        microwave_consumption(idx:idx+29) =  microwave_pattern;

    
%======================================================================
%                           Breakfast appliances pattern
%======================================================================
        p1 = 3e03; % watts (power consumption during the breakfast-on period)
        p2 = 0; % watts (power consumption during the breakfast-off period)
        breakfast_pattern = zeros(30,1);
        if Breakfast_Probability(i) > usage_probability 
            breakfast_pattern(1:3) = p1; % kettle on for 3minutes
            Breakfast_Probability(i+1:end) = 0; % set all future hairdryer values to 0
        else
             breakfast_pattern(1:30) = 0;   
        end
        
        % Store the 30-minute breakfast consumption pattern in the breakfast_consumption array
        idx = (i-1)*30 + 1; % starting index for this 30-minute interval
        breakfast_consumption(idx:idx+29) =  breakfast_pattern;

 
 end

% Define the time axis
t_1440 = linspace(0, 24, 1440);

short_duration_appliances_demand = waterheater_consumption + dishwasher_consumption +microwave_consumption+...
    breakfast_consumption' + hairdryer_consumption + kettle_consumption + laundry_consumption;

% Plot individual minute consumption with legend
plot(t_1440, waterheater_consumption, 'LineWidth', 2);
hold on;
plot(t_1440, dishwasher_consumption, 'LineWidth', 2);
plot(t_1440, microwave_consumption, 'LineWidth', 2)
plot(t_1440, breakfast_consumption', 'LineWidth', 2);
plot(t_1440, hairdryer_consumption, 'LineWidth', 2);
plot(t_1440, kettle_consumption, 'LineWidth', 2);
plot(t_1440, laundry_consumption, 'LineWidth', 2);

% Plot total minute consumption
% plot(t, total_minute_consumption, 'LineWidth', 2.5);

% Set plot properties
xlim([0 24])
title('Individual minute consumption', 'FontSize', 20)
ylabel('Power (W)', 'FontSize', 15)
xlabel('Time (hr)', 'FontSize', 15)
legend('Water heater', 'Dishwasher', 'Microwave','Breakfast appliances', 'hairdryer appliances', 'Kettle', 'Laundry appliances')
set(gca,'FontSize', 20)


% Plot the daily demand  consumption
total_demand =long_durationappliances_demand+short_duration_appliances_demand;

hrs_tranformed = linspace(1,24,1440); % Convert minutes to hours

plot(hrs_tranformed,total_demand, 'LineWidth',1)
xlim([1 24])

title('Total Consumption (Minute-by-minute)', 'FontSize', 20)
xlabel('Time (hours)', 'FontSize', 20);
ylabel('Power consumption (W)', 'FontSize', 20);
% legend('Minute-by-minute consumption', 'FontSize', 20)
hold off

%======================================================================
%                           Adjustment to elexon demand profile
%======================================================================
t_48_48=linspace(1,48,48);
% define the spline interpolation
base_power=interp1(t_48_48, base, hrs_tranformed);


% Calculate the total demand as the sum of the base power and traffic flow
demand_before_adjust= total_demand'+ base_power;% Demand Predicted +Elexon base_power;
demand = base_power;% Elexon base_power;

time=1:1440;
time_intervals = reshape(time, 30, [])';
demand_intervals = reshape(demand(1:end), 30, [])';
demand_adjust_intervals = reshape(demand_before_adjust(1:end), 30, [])';
figure
plot(t_1440,  total_demand)
hold on
plot(t_1440, demand_before_adjust)
hold off
legend ('Base  demand', 'Base demand+predicted demand')


% Calculate the area under the curve for each interval using trapz
difference_areas=zeros(size(demand_intervals, 1), 1);
demand_areas = zeros(size(demand_intervals, 1), 1);
demand_adjust_areas = zeros(size(demand_adjust_intervals, 1), 1);

for i = 1:size(demand_intervals, 1)
    demand_areas(i) = trapz(time(1:30), demand_intervals(i,:));
    demand_adjust_areas(i) = trapz(time(1:30), demand_adjust_intervals(i,:));
    difference_areas(i)=demand_areas(i)-demand_adjust_areas(i) ;
end

% Check if the areas match for each interval
% set tolerance for area difference
tolerance = 0.1; 
for i = 1:48    
    if abs(demand_adjust_areas(i) - demand_areas(i)) > tolerance
        % adjust demand_adjust vector by a factor until areas match
        factor = demand_areas(i) / demand_adjust_areas(i);
        while abs(demand_adjust_areas(i)*factor - demand_areas(i)) > tolerance
            factor = factor + 0.001; % increase factor by small increment
        end
        demand_adjust_intervals(i,:) = demand_adjust_intervals(i,:) * factor;
        demand_adjust_areas(i) = trapz(time_intervals(i,:), demand_adjust_intervals(i,:));
        demand_after_adjust = reshape(demand_adjust_intervals', [], 1); 
    end
end



%======================================================================
%                           Save all files obtained
%======================================================================

% if weekday && strcmp(season, 'winter')
%     demand_adjust_weekday=demand_after_adjust;
%     save demand_adjust.mat demand_adjust_weekday
% 
% elseif weekday && strcmp(season, 'summer')
%     demand_adjust_weekday_summer=demand_after_adjust;
%     save demand_adjust_summer.mat demand_adjust_weekday_summer 
% 
% elseif weekend && strcmp(season, 'winter')
%     demand_adjust_weekend=demand_after_adjust;
%     save demand_adjust_weekend.mat demand_adjust_weekend 
% 
% elseif weekend && strcmp(season, 'summer')
%      demand_adjust_weekend_summer=demand_after_adjust;
%     save demand_adjust_weekend_summer.mat demand_adjust_weekend_summer
% 
% else
%     day_season = 0; % Invalid model
%     error('Invalid date. Please choose another day.');
% end

   

