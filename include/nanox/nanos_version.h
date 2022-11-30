#ifndef  _NANOS_enabled_H_
#define  _NANOS_enabled_H_
#define NANOS_API_MASTER 5041
#define NANOS_API_WORKSHARING 1000
#define NANOS_API_DEPS_API 1001
#define NANOS_API_COPIES_API 1005
#define NANOS_API_TASK_REDUCTION 1002
#define NANOS_API_OPENMP 8
#define NANOS_API_INSTRUMENTATION_API 1001
#define NANOS_API_RESILIENCY 1000
#define NANOS_API_OPENCL 1003
#ifdef _MERCURIUM
#ifndef _NANOS_INTERNAL
#pragma nanos interface family(master) version(5041)
#pragma nanos interface family(worksharing) version(1000)
#pragma nanos interface family(deps_api) version(1001)
#pragma nanos interface family(copies_api) version(1005)
#pragma nanos interface family(task_reduction) version(1002)
#pragma nanos interface family(openmp) version(8)
#pragma nanos interface family(instrumentation_api) version(1001)
#pragma nanos interface family(resiliency) version(1000)
#pragma nanos interface family(opencl) version(1003)
#endif /* _NANOS_INTERNAL */
#endif /* _MERCURIUM */
#endif /* _NANOS_enabled_H_ */
