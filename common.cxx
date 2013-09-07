#include "common.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
//
//  command line option processing
//
int FindOption(int argc, char **argv, const char *option)
{
  for (int i = 1; i<argc; i++)
    if(std::strcmp(argv[i],option)==0)
      return i;
  return -1;
}

int ReadIntOption(int argc, char **argv, const char *option, int default_value)
{
  int iplace = FindOption(argc,argv,option);
  if (iplace>=0 && iplace<argc-1)
    return std::atoi(argv[iplace+1]);
  return default_value;
}

char *ReadStringOption(int argc, char **argv, const char *option, char *default_value)
{
  int iplace = FindOption(argc,argv,option);
  if (iplace>=0 && iplace<argc-1)
    return argv[iplace+1];
  return default_value;
}


