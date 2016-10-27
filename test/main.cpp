#include <iostream>
#include "AceEvalCpp.hpp"

using std::cout;
using std::endl;

int main(int argc, char **argv) {
  cout << "starting ..." << endl;
  OnlineEngine e (argv[1], argv[2]);
  Evidence ev(e);
  cout << "engine & evidence created!" << endl;
  
  e.assertEvidence(ev, false);
  cout << e.probOfEvidence() << endl;
  Variable x1 = e.varForName("x1");
  ev.varCommit(x1, 0);
  e.assertEvidence(ev, false);
  cout << e.probOfEvidence() << endl;
  return 0;
}
