#include "Timer.hpp"

void Timer::start(const std::string& label) {
    starts[label] = Clock::now();
}

void Timer::stop(const std::string& label) {
    auto end = Clock::now();
    auto elapsed = std::chrono::duration_cast<Microseconds>(end - starts[label]);
    times[label] += elapsed;
}

std::string Timer::formatDuration(Microseconds duration) const {
    double time_ms = duration.count() / 1000.0;
    double time_s = time_ms / 1000.0;
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3);
    
    if (time_s >= 1.0) {
        oss << time_s << " s";
    } else {
        oss << static_cast<int>(time_ms) << " ms";
    }
    return oss.str();
}

void Timer::report() const {
    if (times.empty()) {
        return;
    }
    
    std::vector<std::pair<std::string, std::string>> timingData;
    
    if (times.find("Initialization") != times.end()) {
        timingData.emplace_back("Initialization", formatDuration(times.at("Initialization")));
    }
    
    for (const auto& [label, time] : times) {
        if (label != "Initialization") {
            timingData.emplace_back(label, formatDuration(time));
        }
    }
    
    Microseconds totalTime(0);
    for (const auto& [label, time] : times) {
        totalTime += time;
    }
    timingData.emplace_back("TOTAL", formatDuration(totalTime));
    
    logger_.logUserConfig("TIMING SUMMARY", timingData);
}

Timer::Microseconds Timer::get(const std::string& label) const {
    auto it = times.find(label);
    return (it != times.end()) ? it->second : Microseconds::zero();
}

void Timer::reset() {
    starts.clear();
    times.clear();
}