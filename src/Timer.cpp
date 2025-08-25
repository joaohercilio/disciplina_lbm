#include "Timer.hpp"

void Timer::start(const std::string& label) {
    starts[label] = Clock::now();
}

void Timer::stop(const std::string& label) {
    auto end = Clock::now();
    auto elapsed = std::chrono::duration_cast<Microseconds>(end - starts[label]);
    times[label] += elapsed;
}

void Timer::report() const {
    for (const auto& [label, time] : times) {
        std::cout << label << " duration: " << time.count() / 1000 << " ms\n";
    }
}

Timer::Microseconds Timer::get(const std::string& label) const {
    auto it = times.find(label);
    return (it != times.end()) ? it->second : Microseconds::zero();
}

void Timer::reset() {
    starts.clear();
    times.clear();
}