#include <iostream>
#include <boost/filesystem.hpp>

#include "src/wave_constructor.h"

using namespace std;
using namespace riemann;

int main(int argc, char *argv[])
{
  wave_constructor wc;
  wc.test_wave_conditions();
  cout << "[info] done\n";
  return 0;
}
