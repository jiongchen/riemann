#include <iostream>
#include <fstream>

#include "src/json.h"

using namespace std;

int main(int argc, char *argv[])
{
  ifstream ifs("../../dat/karate.json");
  Json::Reader reader;
  Json::Value json;
  reader.parse(ifs, json);
  ifs.close();
  cout << json["cloths"][1]["materials"][0]["strain_limits"][0].isString() << endl;
  cout << "done\n";
  return 0;
}
