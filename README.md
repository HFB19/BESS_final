# BESS_final
Please ensure that the required datasets (dependencies) are available in the working directory before running the script.

## demand_model.m

This MATLAB script calculates the electricity demand for a given date based on various appliances' usage probabilities and power consumption patterns. The script categorizes tariffs based on their pricing structure and usage patterns. It then simulates the power consumption for different appliances and plots the demand profiles.


The implemented code allows users to simulate and analyze patterns over weeks and months. It prompts the user to input a specific date and the corresponding usage probability. By incorporating this information, the code generates a simulated pattern that reflects the expected usage over the specified time period

### Input Parameters
D: The chosen date in the format 'dd-Month-yyyy'.
usage_probability: The probability of appliance usage (0 to 1).

### Dependencies
import_tariff.mat: A dataset containing import tariff information.
basepower.mat: A dataset containing Elexon base power information.


### Code Description
1. Load the necessary datasets (import_tariff.mat and basepower.mat).
2. Check if the chosen date is a weekday or a weekend and determine the season (winter or summer).
3. Determine the appropriate (taiff and Elexon base model) model based on the day type and season.
4. Categorize the import tariffs based on quartiles of their half-hour average.
5. Load the probability datasets for different appliances.
6. Calculate the demand for each 30-minute interval based on the usage probabilities and power consumption patterns.
7. Plot the demand profile for long-duration appliances.
8. Interpolate the data to obtain minute-by-minute consumption profiles for short-duration appliances.
9. Simulate the power consumption patterns for different appliances.
10. Interpolate the total demand profile on a minute-by-minute basis.
11. Adjust the demand model area to match the Elexon demand profile on a half-hourly basis.
12. Plot the total consumption profile.



## supply.m

This MATLAB script simulates the energy generation of solar panels based on various parameters and weather conditions. It calculates the total energy generated throughout the day and plots the results.

## Prerequisites

- MATLAB software installed
- Input file 'WEATHER_FACTOR.mat' containing weather factor data for the simulation


### Input Parameters

- Modify the following parameters in the code according to your requirements:

  - `D`: Input date in the format 'DD-MMM-YYYY'.
  - `condition`: Weather condition, either 'sunny' or 'cloudy'.
 

### Dependencies
The simulation requires weather factor data for different hours of the day. The data should be provided in the 'WEATHER_FACTOR.mat' file.


### Simulation Parameters

- The following simulation parameters can be adjusted:

  - `h_start`: Start of the day (in hours).
  - `h_end`: End of the day (in hours).
  - `Tnom`: Nominal operating temperature (in degrees Celsius).

### Solar Panel Parameters

- Modify the following parameters to match your solar panel setup:

  - `pollution`: Pollution factor.
  - `inverter_efficiency`: Inverter efficiency.
  - `I_o`: Solar constant.
  - `theta_lat`: Latitude.
  - `theta_tilt`: Tilt angle.
  - `panel_number`: Number of panels.
  - `panel_area`: Panel area.
  - `panel_tilt`: Panel tilt angle.
  - `panel_azimuth`: Panel azimuth angle.
  - `panel_efficiency`: Panel efficiency.
  - `house_idle`: House idle power consumption.

- 

### Code Description

1. Ensure that the 'WEATHER_FACTOR.mat' file is in the same directory as the MATLAB script.
2. Open the MATLAB script in MATLAB software.
3. Modify the input and simulation parameters if desired.
4. Run the script.
5. The script will generate a plot showing the total energy generated by the solar panels throughout the day.



## battery_optimised.m

### Input Parameters

- Modify the following parameters in the code according to your requirements:

  - `D`: Input date in the format 'DD-MMM-YYYY'.

 
### Dependencies


supply_data.mat: Contains minute-by-minute solar supply data.
export_tariff.mat: Contains the export tariff data.
import_tariff.mat: Contains the import tariff data.
combined_demand_adjust.mat: Contains the combined demand adjustment data.

### Code Description

1. Loads the required data files.
2. The user  input a specific date.
3. Determines if the input date is a weekday or weekend and identifies the season (winter or summer).
4. Based on the day type and season, selects the appropriate model and corresponding tariff rates.
5. Categorizes the import and export tariffs based on their pricing structure and usage patterns.
6. Sets simulation parameters such as start and end times, time step, and initializes variables.
7. Simulates the battery behavior and energy exchange with the grid for each time step.
8. Calculates the profit, import, export, and net power values for each day of the simulation.
9. Outputs the final results, including the total profit value and various arrays for analysis.
  
## battery_weekly_pattern.m

This code is similar to battery_optimized.m but has additional functionality to simulate battery operation over a weekly and monthly pattern. It performs similar tasks of optimizing battery usage and energy exchange in a solar-powered system, but extends the simulation to cover longer time periods.
