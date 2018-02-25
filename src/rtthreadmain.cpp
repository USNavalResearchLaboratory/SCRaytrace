// $Id: rtthreadmain.cpp,v 1.2 2010-08-16 20:44:20 colaninn Exp $


#include <iostream>
#include "rtthread.h"
#include <cstdlib>

int main(int argc, char *argv[])
{
  std::cout << "In main of rtthreadmain.cpp" << std::endl;
  
  int out=rtthreadtest();

  std::cout << "out : " << out << std::endl;
  
  
  return EXIT_SUCCESS;
}
