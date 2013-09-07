
#ifndef LIST
#define LIST

#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm>

#include "Named.h"
#include "TMath.h"

class Variable;
class AbsPdf;

template<class T> class List : public Named {
 public:
 
 List() : Named("list","list") { ResetIterator(); }
 explicit List(T &element1) : Named("list","list") { AddElement(element1); ResetIterator(); }
 List(T &element1, T &element2) : Named("list","list") { AddElement(element1); AddElement(element2); ResetIterator(); }
 List(T &element1, T &element2, T &element3) : Named("list","list") {
    AddElement(element1); AddElement(element2); AddElement(element3);
    ResetIterator();
  }
  virtual ~List() { }

  Bool_t AddElement(T &element) {
    m_list.push_back(&element);
    m_resetIter = false;

    return true;
  }
  
  Bool_t add(T &element) {
    return AddElement(element);
  }

 private:
  typedef typename std::vector<T*>::iterator tIterator;

 public:

  Bool_t AddElement(List<T> &elements) {
    tIterator iter = elements.m_list.begin();
    while (iter!=elements.m_list.end()) {
      m_list.push_back(*iter);
      iter++;
    }
    m_resetIter = false;
    return true;
  }

  Bool_t add(List<T> &elements) {
    return AddElement(elements);
  }
  
  // no thread-safe
  inline void ResetIterator() { m_iter = m_list.begin(); m_resetIter = true; }
  
  // no thread-safe
  inline T *Next() {
    T *element(0);
    
    if (m_resetIter && m_iter!=m_list.end()) {
      element = *m_iter;
      m_iter++;
    }

    return element;
  }
  
  // thread-safe iterator
  // Always start for the initial element
  class Iterator {
  public: 
    Iterator(std::vector<T*>& elements) : m_iterator(elements.begin()), m_end(elements.end()) { }
    Iterator(const Iterator& iterator) : m_iterator(iterator.m_iterator), m_end(iterator.m_end) { }

    inline T *Next() {
      T *element(0);
    
      if (m_iterator!=m_end) {
	element = *m_iterator;
	m_iterator++;
      }

      return element;
    }

  private:
    tIterator m_iterator;
    tIterator m_end;
  };

  inline Iterator GetIterator() {
    return Iterator(m_list);
  }

  UInt_t GetSize() const { return m_list.size(); }
  
  Int_t Index(const T &element) const {
    Int_t size(m_list.size());
    for (Int_t pos(0); pos<size; pos++) {
      if (0==std::strcmp(element.GetName(),m_list[pos]->GetName()))
	return pos;
    }
    
    return -1; // not found
  }
  
  
  T* GetElement(const UInt_t index) const {
    if (index >= m_list.size())
      return 0;
    return m_list[index];
  }

  void Print(const Bool_t onlyFloat = kFALSE) {
    tIterator iter = m_list.begin();
    std::cout << "# Total Parameters = " << GetSize() << std::endl;
    Int_t nSize(0);
    while (iter!=m_list.end()) {
      if (!onlyFloat || !((*iter)->IsConstant())) {
	(*iter)->Print();
	nSize++;
      }
      iter++;
    }
    std::cout << "# Float Parameters = " << nSize << std::endl;

  }
  
  void Sort() {
    std::sort(m_list.begin(),m_list.end(),sortfunc);
  }

 private:

  static Bool_t sortfunc(T* elA, T* elB) {
    return std::strcmp(elA->GetName(),elB->GetName())<0;
  }

  std::vector<T*> m_list;
  tIterator m_iter; //!
  Bool_t m_resetIter;
  
};

#endif
