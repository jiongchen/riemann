#ifndef HIGH_RESOLUTION_TIMER_H
#define HIGH_RESOLUTION_TIMER_H

#include <chrono>

namespace riemann {

class high_resolution_timer
{
public:
  typedef std::chrono::high_resolution_clock clock_type;
  void start() {
    t0_ = clk_.now();
  }
  void stop() {
    t1_ = clk_.now();
  }
  void log() const {
    uint64_t period = std::chrono::duration_cast<std::chrono::milliseconds>(t1_-t0_).count();
    if ( period/1000 > 0 )
      printf("[info] time duration: %.3lf sec\n", period/1000.0);
    else
      printf("[info] time duration: %.3lf msec\n", period/1.0);
  }
private:
  std::chrono::high_resolution_clock clk_;
  clock_type::time_point t0_, t1_;
};

}
#endif
