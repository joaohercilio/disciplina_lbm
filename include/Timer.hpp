#pragma once

#include <chrono>
#include <string>
#include <unordered_map>
#include <iostream>

class Timer {
public:
    using Clock = std::chrono::high_resolution_clock;
    using Microseconds = std::chrono::microseconds;

    void start(const std::string& label);
    void stop(const std::string& label);
    void report() const;
    Microseconds get(const std::string& label) const;
    void reset();

private:
    std::unordered_map<std::string, Clock::time_point> starts;
    std::unordered_map<std::string, Microseconds> times;
};
