#include <chrono>

class Timer {
private:
    double elapsed_microsec;
    std::chrono::system_clock::time_point start_time;
    int count;

public:
    void reset() {
        elapsed_microsec = 0;
        count            = 0;
    }

    Timer() {
        reset();
    }

    void start() {
        start_time = std::chrono::system_clock::now();
    }

    void stop() {
        std::chrono::system_clock::time_point stop_time = std::chrono::system_clock::now();
        elapsed_microsec += std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count();
        count++;
    }

    double calc_ave_microsec() {
        return elapsed_microsec / (double)count;
    }
};