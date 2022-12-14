/*************************************************************************************/
/*      Copyright 2009-2018 Barcelona Supercomputing Center                          */
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
/*      along with NANOS++.  If not, see <https://www.gnu.org/licenses/>.            */
/*************************************************************************************/

#ifndef _NANOS_PTHREAD
#define _NANOS_PTHREAD

/**@ @**/
//#define _GNU_SOURCE
#include <sched.h>
/**@ @**/

#include "pthread_decl.hpp"
#include "smpprocessor.hpp"



namespace nanos {

inline int PThread::getCpuId() const { return _core->getBindingId(); }

inline size_t PThread::getStackSize () { return _stackSize; }

inline void PThread::setStackSize( size_t size ) { _stackSize = size; }

} // namespace nanos

#endif
