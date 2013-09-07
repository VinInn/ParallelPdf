#ifndef COMMON_H_
#define COMMON_H_

#include <iomanip>
#include <iostream>
#include <cstring>

//
//  argument processing routines
//

int FindOption(int argc, char **argv, const char *option);
int ReadIntOption(int argc, char **argv, const char *option, int default_value = 1000);
char *ReadStringOption(int argc, char **argv, const char *option, char* default_value = 0);

#endif
