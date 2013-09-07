/*

Turbo testing/setting tool for select Intel CPUs
(C) Andrzej Nowak, CERN openlab 2011 <andrzej.nowak@cern.ch>

*/

#include <stdio.h>
#include <fcntl.h>
#include <assert.h>
#include <stdint.h>  
#include <stdlib.h>
#include <getopt.h>

#define VERSION "1.1"

#define MSR_PLATFORM_INFO		0xCE
#define IA32_MISC_ENABLE		0x1A0
#define MSR_TURBO_RATIO_LIMIT 		0x1AD

#define TURBO_ENABLE_BIT 38

#define CPUID_AVX_CAPABILITY		28
#define CPUID_SSE42_CAPABILITY		20
#define CPUID_SSE41_CAPABILITY		19
#define CPUID_AESNI_CAPABILITY		25

#define single_bit(arg) (1ULL<<arg)

#define set_bit(val, bit) 	((val) |= (single_bit(bit)))
#define clear_bit(val, bit)	((val) &= ~(single_bit(bit)))

#define bit_isset(val, bit)	((val >> bit) & 1)

// ARCH capabilities
#define ARCH_CAP_NHM 0x1
#define ARCH_CAP_WSM 0x2
#define ARCH_CAP_SNB 0x4	// Vol 3C 34-1, Oct 2011

//#define VARNAME(arg) #arg

enum arch_enum {
  ARCH_NHM,
  ARCH_WSM,
  ARCH_SNB
};

char *arch_strings[] = {
  "Nehalem",
  "Westmere",
  "Sandy Bridge"
};

double fmul[] = {
  133.33,	// NHM
  133.33,	// WSM
  100.0		// SNB
};

struct option opt[] = {
  { "help", 0, 0, 'h' },
  { "set-turbo", 1, 0, 't' },
  { "core", 1, 0, 'c' },
  { 0, 0, 0, 0 }
};

uint64_t bits(uint64_t val, int hbit, int lbit) {
  int sum = hbit - lbit + 1;
  val >>= lbit;
  val &= (1ULL << sum) - 1;
  return val;
}

#define read_msr(msrcode) read_msr_ext(msrcode, core);

uint64_t read_msr_ext(int msrcode, int cpu) {
  char msr_dev[50];
  int msr_fd;
  uint64_t val;
  
  sprintf(msr_dev, "/dev/cpu/%d/msr", cpu);
  msr_fd = open(msr_dev, O_RDONLY);
  
  assert(msr_fd >= 0);

  pread(msr_fd, &val, sizeof(val), msrcode);
  
  close(msr_fd);
  
  return val;
}

#define write_msr(msrcode, val) write_msr_ext(msrcode, val, core);

void write_msr_ext(int msrcode, uint64_t val, int cpu) {
  char msr_dev[50];
  int msr_fd;
  int retval = -1;
    
  sprintf(msr_dev, "/dev/cpu/%d/msr", cpu);
  msr_fd = open(msr_dev, O_WRONLY);
      
  assert(msr_fd >= 0);
        
  retval = pwrite(msr_fd, &val, sizeof(val), msrcode);
  if(retval != sizeof(val)) {
    printf("Writing value %p to register %d failed\n", val, msrcode);
  }
          
  close(msr_fd);           
}

int get_cpu_arch() {
  int cpuid_eax01_ecx_value = -1;

  __asm__ __volatile__ ("cpuid": "=c" (cpuid_eax01_ecx_value): "a" (1));
  
//  printf("%p\n", cpuid_eax01_ecx_value);    
  if(bit_isset(cpuid_eax01_ecx_value, CPUID_AVX_CAPABILITY))
    return ARCH_SNB;
    
  if(bit_isset(cpuid_eax01_ecx_value, CPUID_AESNI_CAPABILITY))
    return ARCH_WSM;
    
  if(bit_isset(cpuid_eax01_ecx_value, CPUID_SSE41_CAPABILITY))
    return ARCH_NHM;
    
  return -1;
}

void print_help(char *argv[]) {
        printf("Usage: %s [options]\n", argv[0]);
        printf(
        "	--help			-h	Print help (this message)\n"
        "	--set-turbo FLAG	-t	Set turbo (1/on or 0/off)\n"
        "	--core #		-c	Select core (default: 0)\n"
        "	--version		-V	Prints version\n"
        );
}

int main(int argc, char *argv[]) {
  uint64_t mask = 0ULL;
  int opt_index = 0;
  int retval = 0;
  
  int set_turbo = -1;
  int core = 0;
  int arch = get_cpu_arch();
  uint64_t val;
  
  while(1) {
    retval = getopt_long(argc, argv, "ht:c:V", opt, &opt_index);
    if (retval == -1)
      break;
      
    switch(retval) {
      case 'h':
        print_help(argv);
        exit(0);
      case 'V':
        printf("Turbo manipulation tool for Intel CPUs, version %s\n", VERSION);
        printf("(C) 2011 CERN\n\n");
        printf("Written by Andrzej Nowak\n");
        exit(0);
      case 't':
        if (!strncmp(optarg, "on", 2)) {
           set_turbo = 1;
        } else if(!strncmp(optarg, "off", 3)) {
            set_turbo = 0;
        } else {
          set_turbo = strtoul(optarg, NULL, 0);
        }
        break;
      case 'c':
        core = strtoul(optarg, NULL, 0);
        printf("Active core set to core #%d\n", core, core);
        break; 
      default:
        print_help(argv);
        exit(0);
    }
  }  
    
  mask = single_bit(TURBO_ENABLE_BIT);
  printf("The detected architecture is %s\n\n", arch_strings[arch]);

  val = read_msr(IA32_MISC_ENABLE);
  printf("Turbo state: %s\n", val & mask ? "OFF" : "ON");
  val = read_msr(MSR_PLATFORM_INFO);
  int nombins = bits(val, 15, 8);
  printf("Nominal frequency: %.2f GHz\n", nombins*fmul[arch]/1e3);

  val = read_msr(MSR_TURBO_RATIO_LIMIT);
  int bin1 = bits(val, 7, 0);
  int bin2 = bits(val, 15, 8);
  int bin3 = bits(val, 23, 16);
  int bin4 = bits(val, 31, 24);
  
  if (arch == ARCH_NHM) {
      printf("Turbo bins for 1/2/3/4 cores active: %d/%d/%d/%d\n", bin1-nombins, bin2-nombins, bin3-nombins, bin4-nombins);
  } else if (arch == ARCH_WSM) {
      int bin5 = bits(val, 39, 32);
      int bin6 = bits(val, 47, 40);
      printf("Turbo bins for 1/2/3/4/5/6 cores active: %d/%d/%d/%d/%d/%d\n", bin1-nombins, bin2-nombins, bin3-nombins, bin4-nombins, bin5-nombins, bin6-nombins);
  } else if (arch == ARCH_SNB) {
      int bin5 = bits(val, 39, 32);
      int bin6 = bits(val, 47, 40);
      int bin7 = bits(val, 55, 48);
      int bin8 = bits(val, 63, 56);
      printf("Turbo bins for 1/2/3/4/5/6/7/8 cores active: %d/%d/%d/%d/%d/%d/%d/%d\n", bin1-nombins, bin2-nombins, bin3-nombins, bin4-nombins, bin5-nombins, bin6-nombins, bin7-nombins, bin8-nombins);  
  } else
      printf("Turbo bins could not be detected (unknown architecture?)\n");

  printf("\n");
  
  if(set_turbo > -1) {
    val = read_msr(IA32_MISC_ENABLE);
    if (set_turbo == 0) {
      set_bit(val, 38);
    }
    if (set_turbo == 1) {
      clear_bit(val, 38);
    }
    write_msr(IA32_MISC_ENABLE, val);

    val = read_msr(IA32_MISC_ENABLE);
    printf("Turbo state is now: %s\n", val & mask ? "OFF" : "ON");
  }
}
