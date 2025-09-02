#pragma once

#include <chrono>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <sstream>

#include "Logger.hpp"

/**
 * @brief Simple timer class for measuring execution time of code blocks.
 *
 * This class allows timing multiple labeled code sections using high-resolution
 * clocks. It supports starting, stopping, querying, resetting, and reporting
 * the elapsed times for each label with formatted output using Logger.
 */
class Timer {

public:
    
    using Clock = std::chrono::high_resolution_clock; ///< High-resolution clock type
    using Microseconds = std::chrono::microseconds;   ///< Duration type in microseconds

    /**
     * @brief Constructs a Timer with a reference to a Logger.
     * @param logger Reference to Logger for formatted output
     */
    Timer(Logger& logger) : logger_(logger) {}

    /**
     * @brief Starts the timer for a given label.
     * @param label Identifier for the timed section
     */
    void start(const std::string& label);

    /**
     * @brief Stops the timer for a given label and accumulates elapsed time.
     * @param label Identifier for the timed section
     */
    void stop(const std::string& label);

    /**
     * @brief Prints a formatted summary report using the Logger.
     */
    void report() const;

    /**
     * @brief Returns the accumulated time for a given label.
     * @param label Identifier for the timed section
     * @return Elapsed time in microseconds
     */
    Microseconds get(const std::string& label) const;

    /**
     * @brief Resets all timers and accumulated times.
     */
    void reset();

private:
    Logger& logger_; ///< Reference to Logger for output
    std::unordered_map<std::string, Clock::time_point> starts; ///< Start times for each label
    std::unordered_map<std::string, Microseconds> times;       ///< Accumulated times for each label
    
    /**
     * @brief Formats a duration into a human-readable string.
     * @param duration Duration to format
     * @return Formatted string with appropriate units
     */
    std::string formatDuration(Microseconds duration) const;
};