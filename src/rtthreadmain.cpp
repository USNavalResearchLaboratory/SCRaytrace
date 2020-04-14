/** \file rtthreadmain.cpp
 * \brief Provides terminal command line executable for testing the multi-threading.
 */
 

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
