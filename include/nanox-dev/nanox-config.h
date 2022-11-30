/*************************************************************************************/
/*      Copyright 2017 Barcelona Supercomputing Center                               */
/*                                                                                   */
/*      This file is part of the NANOS++ library.                                    */
/*                                                                                   */
/*      NANOS++ is free software: you can redistribute it and/or modify              */
/*      it under the terms of the GNU Lesser General Public License as published by  */
/*      the Free Software Foundation, either version 3 of the License, or            */
/*      (at your option) any later version.                                          */
/*                                                                                   */
/*      NANOS++ is distributed in the hope that it will be useful,                   */
/*      but WITHOUT ANY WARRANTY; without even the implied warranty of               */
/*      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                */
/*      GNU Lesser General Public License for more details.                          */
/*                                                                                   */
/*      You should have received a copy of the GNU Lesser General Public License     */
/*      along with NANOS++.  If not, see <http://www.gnu.org/licenses/>.             */
/*************************************************************************************/

/* Devices */
#ifndef NO_GPU_DEV
#define NO_GPU_DEV
#endif

#ifndef NO_OpenCL_DEV
#define NO_OpenCL_DEV
#endif

#ifndef NO_CLUSTER_DEV
#define NO_CLUSTER_DEV
#endif

#ifndef NO_FPGA_DEV
#define NO_FPGA_DEV
#endif

/* External libraries */
#ifndef NO_HWLOC
#define NO_HWLOC
#endif

#ifndef NO_DLB
#define NO_DLB
#endif

/* Internal features */
#ifndef NO_NANOS_MEMTRACKER_ENABLED
#define NO_NANOS_MEMTRACKER_ENABLED
#endif

#ifndef NO_NANOS_RESILIENCY_ENABLED
#define NO_NANOS_RESILIENCY_ENABLED
#endif

/* Architecture */
#ifndef NO_HAVE_NEW_GCC_ATOMIC_OPS
#define NO_HAVE_NEW_GCC_ATOMIC_OPS
#endif

#ifndef SMP_SUPPORTS_ULT
#define SMP_SUPPORTS_ULT
#endif
