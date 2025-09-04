#include "Logger.hpp"

Logger::Logger() : valuesPerLine_(10), count_(0) {}

void Logger::logStep(int t, int tTotal) {
    int logInterval = std::max(100, tTotal / 100);
    
    if(t % logInterval != 0) return;
    
    std::cout << std::setw(8) << t;
    count_++;
    
    if(count_ % valuesPerLine_ == 0) {
        std::cout << std::endl;
    } else {
        std::cout << " ";
    }
    
    std::cout.flush();
}

void Logger::logMessage(const std::string& msg) {
    std::cout << msg << std::endl;
}

std::string Logger::to_string(double v) {
    std::ostringstream oss;
    oss.precision(17);   
    oss << v;
    return oss.str();
}

void Logger::endLine() {
    if(count_ % valuesPerLine_ != 0) {
        std::cout << std::endl;
    }
    count_ = 0;
    std::cout << std::endl;
}

void Logger::printUserConfigTable(const std::string& title,
                                 const std::vector<std::pair<std::string, std::string>>& params,
                                 bool isNumeric) {
    if(params.empty()) {
        return;
    }

    size_t maxNameLength = 0;
    size_t maxValueLength = 0;
    
    for(const auto& [name, value] : params) {
        maxNameLength = std::max(maxNameLength, name.size());
        maxValueLength = std::max(maxValueLength, value.size());
    }

    const int nameWidth = static_cast<int>(maxNameLength) + 2;
    const int valueWidth = static_cast<int>(maxValueLength) + 2;
    const int totalWidth = nameWidth + valueWidth;
    
    const int separatorWidth = totalWidth + 10;
    const std::string separator(separatorWidth, '=');
    
    const int timeStepWidth = 7 * valuesPerLine_;
    const int offset = (timeStepWidth - separatorWidth) / 2;
    const std::string padding(offset > 0 ? offset : 0, ' ');

    std::cout << std::endl << padding << separator << std::endl;
    
    int titlePadding = (separatorWidth - title.size()) / 2;
    if(titlePadding < 0) titlePadding = 0;
    std::cout << padding << std::string(titlePadding, ' ') << title << std::endl;
    
    std::cout << padding << separator << std::endl;

    for(const auto& [name, value] : params) {
        int linePadding = (separatorWidth - totalWidth) / 2;
        std::string linePaddingStr(linePadding > 0 ? linePadding : 0, ' ');
        
        std::cout << padding << linePaddingStr 
                  << std::left << std::setw(nameWidth) << name
                  << std::right << std::setw(valueWidth) << value
                  << std::endl;
    }

    std::cout << padding << separator << std::endl << std::endl;
}

void Logger::logUserConfig(const std::string& title,
                          const std::vector<std::pair<std::string, double>>& params) {
    std::vector<std::pair<std::string, std::string>> stringParams;
    
    for(const auto& [name, value] : params) {
        std::string valueStr;
        if (std::abs(value - std::round(value)) < 1e-10) {
            valueStr = std::to_string(static_cast<int>(value));
        } else {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(6) << value;
            valueStr = oss.str();
        }
        stringParams.emplace_back(name, valueStr);
    }
    
    printUserConfigTable(title, stringParams, true);
}

void Logger::logUserConfig(const std::string& title,
                          const std::vector<std::pair<std::string, std::string>>& params) {
    printUserConfigTable(title, params, false);
}