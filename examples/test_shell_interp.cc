#include <iostream>
#include <fstream>

#include "src/shell.h"

using namespace std;
using namespace zjucad::matrix;
using namespace riemann;

int main(int argc, char *argv[])
{
  shell_solver sol;
  sol.temp_test();
  cout << "done\n";
  return 0;
}
