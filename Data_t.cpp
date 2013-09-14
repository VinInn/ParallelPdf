#define private public
#include "Data.h"
#undef private

#include<cassert>
#include<iostream>


void one(int np, size_t nv, size_t size) {


  Data::nPartions=np;

  std::cout << "Partitions " << Data::nPartions
	    << " size " << size
	    << " nvar " << nv
	    << std::endl;

  Data data("bha","bho",size,nv);

  for ( auto i=0U; i!=Data::nPartions; ++i)
    std::cout  << i 
	      << " " << data.m_stride[i]
	      << " " << data.m_size[i]
	      << " " << data.m_capacity[i]
	      << " " << data.m_start[i]
	       << '\n'
      ;
    std::cout << std::endl;

    for (auto i=0U; i!=size; ++i) {
      auto k = data.partition(i);
      assert(k<Data::nPartions);
      assert(i>=data.m_start[k]);
      assert(i-data.m_start[k] < data.m_stride[k]);
      for (auto iv=0U; iv!=nv; ++iv) {
	auto p = data.GetData(iv,i);
	auto d = p-data.m_data[k];
	assert(d>=0);
	size_t ud = d;
	assert(ud<data.m_capacity[k]);
	assert(ud<(iv+1)*data.m_stride[k]);
	assert(ud>=iv*data.m_stride[k]);
      }
    }
}


int main() {

  std::cout<< "test of Data Class" << std::endl;

  one(1,3,4096);
  one(2,3,4096*2);

  one(1,3,123*17*53);
  one(2,3,123*17*53);


  return 0;
}
