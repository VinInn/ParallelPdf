// Thread-safe ostream for parallel applications
// OpenMP
// MPI
// By default only the master thread (id=0) of the master process (id=0) 
// can write on the screen, all other are discarded
// NOTE: CERR is not supported

#ifndef MSG_SERVICE
#define MSG_SERVICE

#include <streambuf>
#include <string>
#include <vector>
#include <sstream>

// Shortcut definitions
#define gcoutMPI parallel_ostream::init(parallel_ostream::kCOUT,0,-1,true).log()
#define pendl    std::endl; parallel_ostream::setmask_previous();

class parallel_ostream : public std::streambuf {
 public:

  enum ostream_type {
    kCOUT=0,
    kCLOG,
    kALL
  }; // kALL must be the last element  

  // the first number is the OpenMP mask: >=0 represents the id of the single thread that can output,
  //                                       <0 sets the number of OpenMP threads and all threads can output
  // the second number is the MPI mask:   >=0 represents the id of the single process that can output,
  //                                       <0 (any value) all processes can output
  struct ostream_mask {
    ostream_mask(int aomp, int ampi) : omp(aomp), mpi(ampi) { }
    int omp;
    int mpi;
  };

  static parallel_ostream& init(ostream_type type, ostream_mask mask, bool silent = false);
  static parallel_ostream& init(ostream_type type = kALL, int omp_mask = 0, int mpi_mask = 0, bool silent = false);
  static void cleanup(ostream_type type = kALL);

  static void setmask(ostream_mask mask, bool silent = false);
  static inline void setmask(int omp_mask = 0, int mpi_mask = 0, bool silent = false) {
    setmask(ostream_mask(omp_mask,mpi_mask),silent);
  }
  static inline ostream_mask getmask() { return m_mask; }

  static inline bool is_initialized(ostream_type type = kALL) { 
    if (type==kALL) { // check if all streams are initialized
      int itype(0);
      for (; itype<type && is_initialized(static_cast<ostream_type>(itype)); itype++);
      return (itype==type);
    }
      
    return Instance(type)!=0;
  }

  static inline void setmask_previous() {
    setmask(m_mask_previous,true);
  }

  inline std::ostream& log() {
    return m_stream;
  }

 protected:
  parallel_ostream(std::ostream& stream);
  virtual ~parallel_ostream();

  void SetNStreams(unsigned int nstreams);

  virtual int_type overflow(int_type nChar = traits_type::eof());
  virtual int sync();

  static parallel_ostream *m_instance[kALL]; // cout, clog
  static inline parallel_ostream*& Instance(ostream_type type) {
    return m_instance[type];
  }

  static ostream_mask m_mask;
  static ostream_mask m_mask_previous;

 private:
  std::ostream& m_stream;
  std::vector<std::string> m_buffer;
  std::streambuf *m_stdBuffer;
  std::stringstream m_label;

};

#endif
