// compile with c++ -std=c++11 -O3 -fopenmp -lnuma testNuma.cpp
/*
setenv GOMP_CPU_AFFINITY '0-27'
c++ -std=c++11 -O3 -fopenmp -lnuma testNuma.cpp
foreach i ( 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27)
foreach j ( 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27)
numactl -C 0-27 ./a.out $i $j 8 | grep time
numactl -C 0-27 ./a.out $i $j 120 | grep time
end
end

*/

#include <numa.h>
#include <omp.h>
#include<iostream>
#include<cstdio>
#include <cstdlib>
#include <atomic>
#include <cassert>
#include <x86intrin.h>
unsigned int taux=0;
inline unsigned long long rdtscp() {
 return __rdtscp(&taux);
}


int main(int na,  char * arg[]) {

  if (na<3) std::cout << "please give two thread numbers and buffer size" << std::endl; 

  const int one = ::atoi(arg[1]);
  const int two = ::atoi(arg[2]);
  const unsigned int bs = ::atoi(arg[3]);

  std::cout << "writing " << bs << "MB " << "in thread " << one << " read in thread " << two << std::endl;

  int * buffer = nullptr;
  unsigned int size = bs*1000000/4;
  bool first=true;

  long long sum=0;
  long long t1=0;

  for (int kk=0; kk<50; ++kk) {

#pragma omp parallel
  {
    // assume each thread will allocate in its own NUMA side
    if (one==omp_get_thread_num()) {
      if (first) {
        std::cout << "in thread " << one << " of " << omp_get_num_threads() << std::endl;
      }
      buffer = (int*)numa_alloc_local(size*4);
      // force the OS to allocate physical memory for the region
      for (auto i=0U; i<size;++i) buffer[i]=i;
    } else {

    }
  }

  assert(buffer);
  assert(10==buffer[10]); 

#pragma omp parallel
  {
    // assume each thread will allocate in its own NUMA side
    if (two==omp_get_thread_num()) {
      if (first) {
        std::cout << "in thread " << two << " of " << omp_get_num_threads() << std::endl;
      }
      t1 -= rdtscp();
      for (auto i=0U; i<size;++i) sum +=buffer[i];
      t1 += rdtscp();
    } else {
  
    }
  }


  numa_free(buffer, size*4);
  buffer = nullptr;
  if (first) {
     std::cout << "sum " << sum << std::endl;
     t1=0;
  }

  first = false;
  sum=0;

  }

  std::cout << "time " << bs << "MB "<< one << ' ' << two << ' ' << double(t1)*1.e-6 << std::endl;
  return 0;
}
