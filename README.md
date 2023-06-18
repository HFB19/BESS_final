# BESS_final
Please ensure that the required datasets (dependencies) are available in the working directory before running the script.

## demand_model.m

This MATLAB script calculates the electricity demand for a given date based on various appliances' usage probabilities and power consumption patterns. The script categorizes tariffs based on their pricing structure and usage patterns. It then simulates the power consumption for different appliances and plots the demand profiles.


## supply.m

## battery_optimised.m


## battery_weekly_pattern.m

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



