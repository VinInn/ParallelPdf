#ifndef NAMED
#define NAMED

#include <string>

// Basic definitions by ROOT
#define kTRUE true
#define kFALSE false

class Named {
 public:

  enum rtti { pdf, var, list, unknown};

  Named() {}
  Named(const char* name, const char* title);

  Named(const char* name, const char* title, rtti w) :
    m_name(name), m_title(title), m_who(w){}

  virtual ~Named() { }

  rtti who() const { return m_who;}
  void setNum(int i) { m_num=i;}
  int num() const { return m_num;}


  inline const char* GetName() const { return m_name.c_str(); }
  inline const char* GetTitle() const { return m_title.c_str(); }
  inline const char* name() const { return m_name.c_str(); }
  inline const char* title() const { return m_title.c_str(); }

  
 private:
  std::string m_name;
  std::string m_title;

  rtti m_who=unknown;
  int m_num=-1;

};

#endif
