#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>

/**
 * @brief Class for logging simulation progress and configuration information.
 * 
 * This class provides utilities for logging time steps, user configuration parameters,
 * and messages to the standard output with formatted output and automatic line breaking.
 */
class Logger {

public:
    /**
     * @brief Constructs a Logger object.
     */
    Logger();

    /**
     * @brief Logs a time step with dynamic interval calculation.
     * 
     * The logging interval is automatically adjusted based on the total number of steps
     * to prevent console overload while maintaining progress visibility.
     * 
     * @param t Current time step
     * @param tTotal Total number of time steps in the simulation
     */
    void logStep(int t, int tTotal);

    /**
     * @brief Logs a general message to the output.
     * 
     * @param msg Message string to be logged
     */
    void logMessage(const std::string& msg);

    /**
     * @brief Ends the current line if incomplete and resets the line counter.
     * 
     * This method ensures proper formatting by breaking incomplete lines and
     * preparing the logger for a new sequence of log entries.
     */
    void endLine();
    
    /**
     * @brief Logs user configuration parameters in a formatted table.
     * 
     * Displays the simulation parameters in a centered, formatted table with
     * automatic differentiation between integer and floating-point values.
     * 
     * @param title Title of the configuration section
     * @param params Vector of key-value pairs containing parameter names and values
     */
    void logUserConfig(const std::string& title,
                      const std::vector<std::pair<std::string, double>>& params);
                      
    void logUserConfig(const std::string& title,
                      const std::vector<std::pair<std::string, std::string>>& params);
private:
    int valuesPerLine_;       ///< Number of values to print per line
    int count_;               ///< Counter for values in the current line
    
    /**
     * @brief Calculates the appropriate logging interval based on total steps.
     * 
     * Determines the step interval for logging to balance information density
     * with console output performance.
     * 
     * @param tTotal Total number of time steps
     * @return Calculated logging interval
     */
    int calculateLogInterval(int tTotal) const;

    void printUserConfigTable(const std::string& title,
                            const std::vector<std::pair<std::string, std::string>>& params,
                            bool isNumeric = false);
};