/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Indicates the presence of the Cluster arch plugin. */
/* #undef CLUSTER_DEV */

/* Disables warnings when using functions deprecated in OpenCL 2.0 */
/* #undef CL_USE_DEPRECATED_OPENCL_2_0_APIS */

/* Enables DLB support. */
/* #undef DLB */

/* Enables FPGA support */
/* #undef FPGA_DEV */

/* Defines which multithreading support is required to GASNet */
/* #undef GASNET_PAR */

/* Indicates the presence of the GPU arch plugin. */
/* #undef GPU_DEV */

/* Define to 1 if you have the <CL/opencl.h> header file. */
/* #undef HAVE_CL_OPENCL_H */

/* Define to 1 if you have the <dlb.h> header file. */
/* #undef HAVE_DLB_H */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <gasnet.h> header file. */
/* #undef HAVE_GASNET_H */

/* Define to 1 if you have the <hwloc.h> header file. */
/* #undef HAVE_HWLOC_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <libxdma.h> header file. */
/* #undef HAVE_LIBXDMA_H */

/* Define to 1 if you have the <memkind.h> header file. */
/* #undef HAVE_MEMKIND_H */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <mpi.h> header file. */
/* #undef HAVE_MPI_H */

/* Defined if the compiler supports the new gcc atomic builtins */
/* #undef HAVE_NEW_GCC_ATOMIC_OPS */

/* Define to 1 if you have the <OpenCL/opencl.h> header file. */
/* #undef HAVE_OPENCL_OPENCL_H */

/* Define if you have POSIX threads libraries and header files. */
#define HAVE_PTHREAD 1

/* Have PTHREAD_PRIO_INHERIT. */
#define HAVE_PTHREAD_PRIO_INHERIT 1

/* Define to 1 if you have the <spi/include/kernel/location.h> header file. */
/* #undef HAVE_SPI_INCLUDE_KERNEL_LOCATION_H */

/* Define to 1 if you have the <spi/include/kernel/process.h> header file. */
/* #undef HAVE_SPI_INCLUDE_KERNEL_PROCESS_H */

/* Define to 1 if you have the <sqlite3.h> header file. */
#define HAVE_SQLITE3_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if the compiler provides the __sync_add_and_fetch function for
   uint64 */
#define HAVE_SYNC_ADD_AND_FETCH_8 1

/* Define to 1 if the compiler provides the __sync_and_and_fetch function for
   uint64 */
#define HAVE_SYNC_AND_AND_FETCH_8 1

/* Define to 1 if the compiler provides the __sync_bool_compare_and_swap
   function for uint64 */
#define HAVE_SYNC_BOOL_COMPARE_AND_SWAP_8 1

/* Define to 1 if the compiler provides the __sync_or_and_fetch function for
   uint64 */
#define HAVE_SYNC_OR_AND_FETCH_8 1

/* Define to 1 if the compiler provides the __sync_sub_and_fetch function for
   uint64 */
#define HAVE_SYNC_SUB_AND_FETCH_8 1

/* Define to 1 if the compiler provides the __sync_xor_and_fetch function for
   uint64 */
#define HAVE_SYNC_XOR_AND_FETCH_8 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Indicates the presence of hwloc library. */
/* #undef HWLOC */

/* BlueGene/Q host compatibility is enabled. */
/* #undef IS_BGQ_MACHINE */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Enables memkind support */
/* #undef MEMKIND_SUPPORT */

/* Ignore cxx seek errors when including mpi.h in C++ */
/* #undef MPICH_IGNORE_CXX_SEEK */

/* Identifies which MPI implementation is being used. Supported values: intel,
   mpich, openmpi */
/* #undef MPI_IMPLEMENTATION */

/* API version of the CUDA package specified by the user */
/* #undef NANOS_CUDA_VERSION */

/* Version of the hwloc package specified by the user */
/* #undef NANOS_HWLOC_VERSION */

/* Specifies whether Nanos++ memtracker module has been enabled */
/* #undef NANOS_MEMTRACKER_ENABLED */

/* Indicates whether resiliency features should be used or not. */
/* #undef NANOS_RESILIENCY_ENABLED */

/* Specifies whether Nanos++ allocator has been enabled */
/* #undef NANOS_USE_ALLOCATOR */

/* Build version */
#define NANOX_BUILD_VERSION "git master b8af14c 2019-02-04 15:30:26 +0100 developer version"

/* Configure line */
#define NANOX_CONFIGURE_ARGS "./configure --prefix=/homelocal/lcostero_local/odroid/REDO_volta/ompss_20162003/nanox_nuevo --with-extrae=/usr/local/ LDFLAGS=-lcpufreq CXXFLAGS=-DVOLTA1 --disable-instrumentation --disable-debug --enable-static"

/* Enables task callback feature */
/* #undef NANOX_TASK_CALLBACK */

/* Indicates the presence of the OpenCL arch plugin. */
/* #undef OpenCL_DEV */

/* Name of package */
#define PACKAGE "nanox"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "pm-tools@bsc.es"

/* Define to the full name of this package. */
#define PACKAGE_NAME "nanox"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "nanox 0.15a"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "nanox"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.15a"

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* Indicates support for user level threads */
#define SMP_SUPPORTS_ULT /**/

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "0.15a"

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */
