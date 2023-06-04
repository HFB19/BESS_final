clear all; close all; clc;
% Load minute-by-minute solar and demand data from csv file

load supply_data.mat;load demand_adjust_summer.mat;load export_tariff.mat;load import_tariff.mat
load demand_adjust.mat;load demand_adjust_summer.mat;load demand_adjust_weekend.mat;load demand_adjust_weekend_summer.mat

%%%%%%%%%%%%%%%%%
% Input date%
%%%%%%%%%%%%%%%%%%
D =' 26-June-2022';%%only choose winter

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
% Calculate the difference in days
NumDays= abs(d2 - d1);

if (NumDays >= 1 && NumDays <= 59) || (NumDays>= 244 && NumDays<= 365)
    season = 'winter';
    disp('You chose a winter week!') 
elseif NumDays >= 60 && NumDays <= 243
    season = 'summer';
    disp('You chose a summer week!') 
else
    season = '';
    error('Invalid date. Please choose another day.');
end 


% Determine if the input date is a weekday or weekend day
% Determine the appropriate model based on weekday/weekend and season
if weekday && strcmp(season, 'winter')
    day_season = 11; % Winter weekday model
    import_tariff=import_tariff_weekday_winter;
    export_tariff=export_tariff_weekday_winter;

elseif weekday && strcmp(season, 'summer')
    day_season = 12; % Summer weekday model
    import_tariff=import_tariff_weekday_summer;
    export_tariff=export_tariff_weekday_summer;

elseif weekend && strcmp(season, 'winter')
    day_season = 21; % Winter weekend mode
    import_tariff=import_tariff_weekend_winter;
    export_tariff=export_tariff_weekend_winter;
  
elseif weekend && strcmp(season, 'summer')
    day_season = 22; % Summer weekend model
     import_tariff=import_tariff_weekend_summer;
     export_tariff=export_tariff_weekend_summer;

else
    day_season = 0; % Invalid model
    error('Invalid date. Please choose another day.');
end


%======================================================================
%                   Concatenating demand and supply data
%======================================================================
week_in_test = 1;
random_integer = randi([1,7]);

% Determine the appropriate model based on weekday/weekend and season
if strcmp(season, 'winter')
    solar_data = horzcat(repmat(supply_winter, 1, random_integer),repmat(supply_winter_cloudy, 1, 7-random_integer ));
    demand_data = horzcat(repmat(demand_adjust_weekday'/1000, 1, 5), repmat(demand_adjust_weekend'/1000, 1, 2));
elseif strcmp(season, 'summer')
    solar_data = horzcat(repmat(supply_summer, 1, random_integer),repmat(supply_summer_cloudy, 1, 7-random_integer ));
  demand_data = horzcat(repmat(demand_adjust_weekday_summer'/1000, 1, 5), repmat(demand_adjust_weekend_summer'/1000, 1, 2));
else
    error('Invalid date. Please choose another day.');
end

% Randomly select the order of sunny and cloudy days in the week

day_in_test=week_in_test*7;



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
import_categories=repelem(import_categories,30);

import_categories_days= repmat(import_categories', 1,day_in_test);


for i = 1:48
   export_half_hour_average(i) =mean(export_tariff(i,1));
end

% Sort the data into quartiles
quartiles = prctile(export_half_hour_average, [0 25 50 75 100]);

% Classify the half_hour_average into 4 import_categories based on quartiles
export_categories = zeros(size(export_half_hour_average));
export_categories(export_half_hour_average <= quartiles(2)) = 1; % lower quartile
export_categories(export_half_hour_average > quartiles(2) & export_half_hour_average <= quartiles(3)) = 2; % middle quartiles
export_categories(export_half_hour_average> quartiles(3) & export_half_hour_average<= quartiles(4)) = 3;
export_categories(export_half_hour_average > quartiles(4)) = 4; % upper quartile
export_categories=export_categories';
export_categories=repelem(export_categories,30);
export_categories_days= repmat(export_categories', 1,day_in_test);


%======================================================================
%                           Simulation Parameters
%======================================================================
start_time = 0;
end_time = 1440*day_in_test; % 24 hours in minutes
time_step = 1; % in minutes
% Initialize variables
n_steps = (end_time - start_time)/time_step;
time = linspace(start_time, end_time, n_steps);
% Initialize demand and solar power vectors
demand_power = zeros(1, n_steps);
solar_power = zeros(1, n_steps);

% Fill demand and solar power vectors
for i = 1:n_steps
    demand_power(i) =demand_data(i);
    solar_power(i) = solar_data(i) ;
end
plot(solar_data)
plot(demand_data)

%======================================================================
%                           Battery Parameters
%======================================================================
battery_capacity=13;
max_charge_rate = 3; % kW
max_discharge_rate = 3; % kW
charge_efficiency =0.9;
discharge_efficiency = 0.9;
max_import_rate = 3; % kW


% Initialize battery variables
battery_energy = zeros(1, n_steps);
battery_power = zeros(1, n_steps);
charge_power = zeros(1, n_steps);
discharge_power = zeros(1, n_steps);
charge_energy = zeros(1, n_steps);
discharge_energy = zeros(1, n_steps);
net_power = zeros(1, n_steps);
battery_energy(1)=0;

%======================================================================
%                           Grid Parameters 
%======================================================================
grid_import_energy = zeros(1, n_steps);
grid_export_energy = zeros(1, n_steps);
grid_import_rate = import_tariff'/100; % $/kWh
grid_export_rate =export_tariff'/100; % $/kWh




%======================================================================
%                           Prediction Parameters
%======================================================================

level_charged_at_night=zeros(1,day_in_test);
level_charged_at_night(1,1) =1;  % Define the battery capacity values you want to test
export_threshold=zeros(1,day_in_test);
export_threshold(1) =1;  % Define the battery capacity values you want to test
profit_daily = zeros(1, day_in_test);
export_daily = zeros(1, day_in_test);
import_daily = zeros(1, day_in_test);
import_daily_after7 = zeros(1, day_in_test);
net_power_daily= zeros(1, day_in_test);
import_dailyafter_6 = zeros(1, day_in_test);
%======================================================================
%                           Troubleshooting Arrays
%======================================================================
caseArray = [];
timeArray = [];
day = 1; % Initialize day counter
total_profit_value=zeros(11,11);
level_charged_at_night(1,1)=0.0;
export_threshold(1,1)=0;



for t = 2:n_steps
    net_power(t) = demand_power(t)-solar_power(t); %+ battery_power(t);
    buy = import_categories_days(t) == 1;
    sell = export_categories_days(t) >= 4;
    % Determine battery charge for discharge during high tariff period
        if net_power( t) > 0


            if battery_energy(t-1) <=level_charged_at_night(1,day)*battery_capacity
                if buy 

                    charge_power(t) = max_charge_rate; % this is linear
                    charge_energy(t) = charge_power(t) / 60;
                    battery_energy(t) = min(battery_energy(t-1) + charge_energy(t), battery_capacity);
                    battery_power(t) = min(charge_power(t) * charge_efficiency, max_charge_rate);
                    grid_import_energy(t) = charge_power(t) + net_power( t);
                    caseArray = [caseArray, 1];
                    %%Charge battery and net_power from grid


            elseif sell
                if mod(t-1, 1440) > (18 * 60)
                    %%% discharge everything to net_power at the end of the day
                    if battery_energy(t-1)~=0
                    discharge_power(t) = min(net_power( t), max_discharge_rate);
                    discharge_energy(t) = discharge_power(t) / 60;
                    battery_energy(t) = max(battery_energy(t-1) - discharge_energy(t) * discharge_efficiency, 0);
                    battery_power(t) = min(discharge_power(t) * discharge_efficiency, max_discharge_rate);
                    grid_export_energy(t) = discharge_power(t) - net_power( t);
                     grid_import_energy(t)=0;
                    else
                        grid_import_energy(t)=net_power(t);
                         grid_export_energy(t)=0;
                    end
                     caseArray = [caseArray, 2];
                elseif battery_energy(t-1) > export_threshold(1,day)*battery_capacity 
                    discharge_power(t) = min(net_power( t), max_discharge_rate);
                    discharge_energy(t) = discharge_power(t) / 60;
                    battery_energy(t) = max(battery_energy(t-1) - discharge_energy(t) * discharge_efficiency, 0);
                    battery_power(t) = min(discharge_power(t) * discharge_efficiency, max_discharge_rate);
                    grid_export_energy(t) = discharge_power(t) - net_power( t);
                     grid_import_energy(t)=0;
                     caseArray = [caseArray, 3];
                else
                    % Battery reaches threshold during discharging, stop discharging and maintain energy
                    battery_energy(t) = battery_energy(t-1);
                    battery_power(t) = 0;
                    grid_export_energy(t) = 0;
                    grid_import_energy(t) = net_power( t);
                    caseArray = [caseArray, 4];
                end



          else 

               if mod(t-1, 1440) > (18 * 60)
                    %%% discharge everything to net_power at the end of the day
                    if battery_energy(t-1)~=0
                    discharge_power(t) = net_power( t);
                    discharge_energy(t) = discharge_power(t) / 60;
                    battery_energy(t) = max(battery_energy(t-1) - discharge_energy(t) * discharge_efficiency, 0);
                    battery_power(t) = min(discharge_power(t) * discharge_efficiency, max_discharge_rate);
                    grid_export_energy(t) =0;
                    grid_import_energy(t)=0;
                    else
                     grid_import_energy(t)=net_power(t);
                     grid_export_energy(t)=0;
                    end
                    caseArray = [caseArray, 5];
                elseif battery_energy(t-1) > 0
                    discharge_power(t) = min(net_power( t), max_discharge_rate);
                    discharge_energy(t) = discharge_power(t) / 60;
                    battery_energy(t) = max(battery_energy(t-1) - discharge_energy(t) * discharge_efficiency, 0);
                    battery_power(t) = min(discharge_power(t) * discharge_efficiency, max_discharge_rate);
                    grid_import_energy(t) = 0;
                    grid_export_energy(t) = 0;
                    caseArray = [caseArray, 6];
                else
                    battery_energy(t) = battery_energy(t-1);
                     battery_power(t)=0;
                    discharge_power(t) = 0;
                    grid_import_energy(t) = net_power( t);
                    caseArray = [caseArray, 7];
                    
                    
                 end
                %%%selfsupply from battery->house
            end

            elseif battery_energy(t-1)>level_charged_at_night(1,day)*battery_capacity
                if buy 
                     grid_import_energy(t)=net_power( t);
                      battery_energy(t) = battery_energy(t-1);
                    battery_power(t) = 0;
                    caseArray = [caseArray, 8];

                elseif sell
                     discharge_power(t) = min(net_power( t),max_discharge_rate)  ;
                     discharge_energy(t) = discharge_power(t) / 60;
                     battery_energy(t) = max(battery_energy(t-1) - discharge_energy(t) * discharge_efficiency, 0);
                     battery_power(t) = min(discharge_power(t) * discharge_efficiency, max_discharge_rate);
%                      grid_import_energy(t)=net_power( t);
                     grid_import_energy(t)=net_power( t);
                     grid_export_energy(t)=discharge_power(t) ;
                     caseArray = [caseArray, 9];
                else
                    if battery_energy(t-1)~=0
                     discharge_power(t) = net_power( t)  ;
                     discharge_energy(t) = discharge_power(t) / 60;
                     battery_energy(t) = max(battery_energy(t-1) - discharge_energy(t) * discharge_efficiency, 0);
                     battery_power(t) = min(discharge_power(t) * discharge_efficiency, max_discharge_rate);
                     grid_import_energy(t)=0;
                     grid_export_energy(t)=0;
                    else
                      grid_import_energy(t)=net_power( t)  ;
                     battery_energy(t) = battery_energy(t-1);
                    battery_power(t) = 0;
                    grid_export_energy(t)=0;
                    
                      
                    end
                    caseArray = [caseArray, 10];
                end

            end

            %%%netpower<0
        elseif net_power( t-1) <= 0
            if battery_energy(t-1) >export_threshold(1,day) * battery_capacity
                if sell
                    if battery_energy(t-1)~=0
                        discharge_power(t) = max_discharge_rate;
                    else
                        discharge_power(t) = 0;
                    end
                    discharge_energy(t) = discharge_power(t) / 60;
                    battery_energy(t) = max(battery_energy(t-1) - discharge_energy(t) * discharge_efficiency, 0);
                    battery_power(t) = min(discharge_power(t) * discharge_efficiency, max_discharge_rate);
                    grid_export_energy(t) = -net_power( t) + discharge_power(t);
                else
                    grid_export_energy(t) = -net_power( t);
                    battery_energy(t) = battery_energy(t-1);
                    battery_power(t) = 0;
                end
                caseArray = [caseArray, 11];
            elseif battery_energy(t-1) <= export_threshold(1,day) * battery_capacity
                if sell
                    battery_energy(t) = battery_energy(t-1);
                    battery_power(t) = 0;
                    grid_export_energy(t) = -net_power( t);
                elseif buy
                    charge_power(t) = max_charge_rate-net_power(t);
                    charge_energy(t) = charge_power(t) / 60;
                    battery_energy(t) = min(battery_energy(t-1) + charge_energy(t), battery_capacity);
                    battery_power(t) = min(charge_power(t) * charge_efficiency, max_charge_rate);
                    grid_import_energy(t) = charge_power(t);
%                     grid_export_energy(t) = -net_power( t);
                else
                    charge_power(t) = -net_power(t);
                    charge_energy(t) = charge_power(t) / 60;
                    battery_energy(t) = min(battery_energy(t-1) + charge_energy(t), battery_capacity);
                    battery_power(t) = min(charge_power(t) * charge_efficiency, max_charge_rate);
                    grid_import_energy(t) = 0;
                     grid_export_energy(t) = 0;
                end
                caseArray = [caseArray, 12];
            end
        else
            battery_energy(t) = battery_energy(t-1);
            battery_power(t) = 0;
            caseArray = [caseArray, 13];
        end
    
    % Store the time value
    timeArray = [timeArray, t];

 


 

    if mod(t, 1440) == 0

        day = ceil(t / 1440); % Update day counter and round up to the nearest integer
        
        % Reshape the data to have 60-minute intervals
        hourly_grid_export_energy = reshape(grid_export_energy(t-1439:t), 60, [])';
        hourly_grid_import_energy = reshape(grid_import_energy(t-1439:t), 60, [])';
        hourly_net_demand = reshape(net_power(t-1439:t), 60, [])';

        % Calculate energy consumed/exported per hour
        hourly_energy_export = zeros(24, 1);
        hourly_energy_import = zeros(24, 1);
        hourly_net_power = zeros(24, 1);
        
        for i = 1:24
            % Calculate the time vector (in minutes) for the current hour
            time_vec = ((i-1)*60 + 1):(i*60);
            
            % Calculate the time interval (in hours) between consecutive points
            dt = diff(time_vec) / 60;  % divide by 60 to convert to hours
            
            % Integrate the half-hourly energy values, multiplied by the time interval, to get the total energy consumed/exported in kWh for the current hour
            hourly_energy_export(i) = trapz(time_vec, hourly_grid_export_energy(i,:) .* [dt, dt(end)]/2); % multiply by dt and use trapezoidal rule
            hourly_energy_import(i) = trapz(time_vec, hourly_grid_import_energy(i,:) .* [dt, dt(end)]/2);
            hourly_net_power(i) = trapz(time_vec,  hourly_net_demand(i,:) .* [dt, dt(end)]);
        end
        
        % Reshape the data to have 30-minute intervals
        halfhourly_grid_export_energy = reshape(grid_export_energy(t-1439:t), 30, [])';
        halfhourly_grid_import_energy = reshape(grid_import_energy(t-1439:t), 30, [])';
        
        % Calculate energy consumed/exported per half hour
        half_hourly_energy_export = zeros(48, 1);
        half_hourly_energy_import = zeros(48, 1);
        
        for i = 1:48
            % Calculate the time vector (in minutes) for the current half hour
            time_vec = ((i-1)*30 + 1):(i*30);
            
            % Calculate the time interval (in hours) between consecutive points
            dt = diff(time_vec) / 60;  % divide by 60 to convert to hours
            
            % Integrate the half-hourly energy values, multiplied by the time interval, to get the total energy consumed/exported in kWh for the current hour
            half_hourly_energy_export(i) = trapz(time_vec, halfhourly_grid_export_energy(i,:) .* [dt, dt(end)]) / 2;  % multiply by dt and use trapezoidal rule
            half_hourly_energy_import(i) = trapz(time_vec, halfhourly_grid_import_energy(i,:) .* [dt, dt(end)]) / 2;
            total_import_cost(i) = half_hourly_energy_import(i) * grid_import_rate(i);
            total_export_profit(i) = half_hourly_energy_export(i) * grid_export_rate(i);
            total_profit(i) = total_export_profit(i) - total_import_cost(i);
        end


      end



     if mod(t, 1440) == 0
        % Update daily variables
        profit_daily(day) = sum(total_profit);
        export_daily(day) = sum(hourly_energy_export);
        import_daily(day)= sum(hourly_energy_import);
        net_power_daily(day)=sum( hourly_net_power);
        import_daily_after7(day) = sum(hourly_net_power(7:24));
        import_dailyafter_6(day) = sum(hourly_net_power(20:24));
        export_threshold(day+1) = min(1,import_dailyafter_6(day) / 13);
        level_charged_at_night(day+1) = max(min(net_power_daily(day) / 13, 1),0);

        % Accumulate import_daily from previous days
        import_daily_cumulative(day) = import_daily_after7(day);
        net_power_daily_cumulative(day)=net_power_daily(day);
        import_dailyafter_6_cumulative(day)=import_dailyafter_6(day);
        n=zeros(1,day_in_test);
        if day > 1
            import_daily_cumulative(day) = import_daily_cumulative(day) + import_daily_cumulative(day-1);
            net_power_daily_cumulative(day) = net_power_daily_cumulative(day) + net_power_daily_cumulative(day-1);
            import_dailyafter_6_cumulative(day) = import_dailyafter_6_cumulative(day) + import_dailyafter_6_cumulative(day-1);
            n(1, day) = import_dailyafter_6_cumulative(day) / day;
            export_threshold(day+1) = min(1, n(1, day) / 13);
            level_charged_at_night(day+1) = max(min(import_daily_cumulative(day) / (day * 13), 1), 0);
        
        else
            import_daily_cumulative(day) = import_daily_cumulative(day);
            net_power_daily_cumulative(day) = net_power_daily_cumulative(day);
            import_dailyafter_6_cumulative(day) = import_dailyafter_6_cumulative(day);
            export_threshold(day+1) = min(1, mean(import_dailyafter_6_cumulative(1:day)) / 13);
            level_charged_at_night(day+1) = max(min(mean(import_daily_cumulative(1:day)) / (day * 13), 1), 0);
        end
        level_charged_at_night(day+1);
         day=day+1;

     end



end
