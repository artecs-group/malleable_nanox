#ifndef _FGCS_MALLEABLE_FREQ_ARCH_SCHED_HPP_
#define _FGCS_MALLEABLE_FREQ_ARCH_SCHED_HPP_

#include "schedule.hpp"
#include "wddeque.hpp"
#include "plugin.hpp"
#include "system.hpp"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <semaphore.h>
#include <set>

#include "client.h"

/*@ @*/
#include <sys/time.h>
#include <cpufreq.h>
#include "common/dynamic_bitset.hpp"
#include "smpplugin_decl.hpp"

#ifdef NANOS_DEBUG_ENABLED
#define ONLY_DBG(f) f;				
#else
#define ONLY_DBG(f) ;				
#endif

// Nyapas para hacer pruebas rapidas. Poner esto bien cuando corresponda
#if defined(KELVIN) || defined(kelvin)
#define MAX_FREQ 4000000
#define MIN_FREQ 1000000
#elif defined(MAKALU) || defined(makalu) || defined (VOLTA1) || defined (volta1)
#define MAX_FREQ 1700000
#define MIN_FREQ 1000000
#else
#error "machine not defined"
#endif


namespace nanos {
namespace ext {    
  
  float getTimeStamp();
	
  enum threadStatus_t{
    NOT_TIED  = (1 << 0), // We don't have any core to run this thread
    RUNNING   = (1 << 1), // The thread is running on its own core
    SLEEPING  = (1 << 2), // There's a core assigned, but there is not any task to execute
    WAKING_UP = (1 << 3)  // This thread is receiving msg from the server. Probably oversuscription
  };
	
  enum appStatus_t {
    START,
    STOP
  };


  class FgcsMalleableFreqArch : public SchedulePolicy {
  private:
    struct TeamData : public ScheduleTeamData {
      WDPriorityQueue<> *_readyQueue;
      TeamData () : ScheduleTeamData(), _readyQueue(NULL) {_readyQueue = NEW WDPriorityQueue<>;}
      ~TeamData () { delete _readyQueue; }
    };
    
    struct { //-----------------------------------------------------------------
      int numCores; //numner of cores (i.e., the number of worker threads created)
		
      //---- LOCK wake up / sleep cores
      Lock  waking_lock;
      sem_t *sems;

      //---- LOCK
      Lock threadStatus_lck;		 
      threadStatus_t *threadStatus;     //NOT_TIED, RUNNING, SLEEPING, WAKING_UP
      bool           *to_untie;         //NOT_TIED, RUNNING, SLEEPING, WAKING_UP
      int            *coreOccupation;   //[coreId] = threadId
      int            *threadAssignment; //[threadId] = coreId
      Lock wakeSleeping_lck;		 
      
      appStatus_t    appStatus;
      dynamic_bitset currMask;
      int            nActiveThreads;

      /*** MALLEABLE FREQUENCY **/
      //---- LOCK:
      Lock   changeInProgress;
      bool   futureChange;
      int    runningFreq_idx;

      FILE*  logFile;
    } _expControl; //-----------------------------------------------------------


    //   1000 -> 2000
    const int _NFREQS;
    double _powers[11];
    const double _powerIdle;
    double _maxBudget;
	    
  public:
    FgcsMalleableFreqArch();
    virtual ~FgcsMalleableFreqArch ();
    
    virtual size_t getTeamDataSize () const { return sizeof(TeamData); }
    virtual size_t getThreadDataSize () const { return 0; }
    
    virtual ScheduleTeamData * createTeamData () { return NEW TeamData(); }
    virtual ScheduleThreadData * createThreadData () { return 0; }
    virtual void queue ( BaseThread *thread, WD &wd );
    virtual void queue ( BaseThread ** threads, WD ** wds, size_t numElems );
    bool isValidForBatch ( const WD * wd ) const { return true;	}

    void atSuccessor   ( DependableObject &successor, DependableObject &predecessor );
    virtual WD *atSubmit ( BaseThread *thread, WD &newWD );
    void atShutdown();
	    
    WD * atIdle ( BaseThread *thread, int numSteal );
    WD * atPrefetch ( BaseThread *thread, WD &current );
    WD * atBeforeExit ( BaseThread *thread, WD &current, bool schedule );

    //virtual WD * atWakeUp      ( BaseThread *thread, WD &wd );


  private:
    /*==========================================================================
     ** Workers/affinity/sleep/wakeUp related  **
     =========================================================================*/

    /** RECEIVES the new masks from the server, continuously. In the case that all the
     ** cores are assigned to this app, the threads stop listening and starts executing
     ** some tasks. **/
    void waitForNewMask(int threadId);

          
    /** Given a new AFFINITY (grater than the current one), assings the new cores
     ** to idle worker threads.
     ** Moves the threads to SLEEPING status, but DO NOT WAKE up any thread   **/
    void increaseAffinity(const dynamic_bitset& bt);


    /** WAKE UP as many threads as possible.
     ** To wake up a thread, it has to be in sleeping state (i.e., to have
     **+a core assigned)                                                      **/
    inline void wakeUp_sleepingWorkers(int numWorkers=1);

	    
    void initializeWaitingThread();
	    
          
    /** FREEs a CORE associated to a worker thread.
     ** If there are threads sleeping, frees one of those cores,
     **+if not, frees its own core.                                          
     ** RETURNS true if a sleeping thread was release. False in other case   **/
    void untie_workerThread(int nToUntie); //If sleeping -> freeSleepingWorker

    bool untie_sleepingWorker();
    void order_untie  ();
    void untie_myself (int threadId);

    
    //void untie_thread (BaseThread *thread);
    
    /** SLEEP a specific thread.         **/
    void sleepThread(BaseThread *thread);



    /*==========================================================================
    ** Frequency related**
     =========================================================================*/
    void adjustFrequency();
    //DO NOT CALL DIRECTLY (assumes the lock is blocked)
      void updateFreqIdx();
      void applyFreq();
    

		
    inline void changeFreq(int threadId, long long freq) const {
      ((SMPPlugin*)sys.getSMPPlugin())->setFrequency(threadId, freq);    }

    inline long long getFreqFromIdx_Hz(int idx) const { return (10+idx)*100000; }

    inline long double estimatePower_idx(unsigned nIdle,
					 unsigned nRunning, unsigned rFreqIdx) const {
      return (nIdle*_powerIdle + (nRunning*_powers[rFreqIdx]));
    }

    inline float getCurrentPowerConsumption() const {
      return (((_expControl.currMask.count() -_expControl.nActiveThreads)*_powerIdle) +
	      (_expControl.nActiveThreads*_powers[_expControl.runningFreq_idx]));
    }
    inline float getDesiredBudget() const {
      return (((_expControl.currMask.count()-_expControl.nActiveThreads)*_powerIdle) +
	      (_expControl.nActiveThreads*_powers[_NFREQS-1]));
    }

    inline float getDesiredBudget(int nThs) const {
      ONLY_DBG ( assert (nThs <= _expControl.numCores) );
      return (((_expControl.currMask.count()-nThs)*_powerIdle) +
	      (nThs*_powers[_NFREQS-1]));
    }

    inline int getMaxThreadsWithBudget(double maxBudget) const {
      return ((maxBudget - _expControl.currMask.count()*_powerIdle) / ( 2.5 - _powerIdle)); //ponemos 2.5 en vez de _powers[0]
    }
  };

	
  class FgcsMalleableFreqArchPlugin : public Plugin {
  public:
    FgcsMalleableFreqArchPlugin() : Plugin( "FGCS Malleable Plugin",1 ) {}

    virtual void config ( Config &cfg ){ }

    virtual void init() {
      sys.setDefaultSchedulePolicy(NEW FgcsMalleableFreqArch());
    }
  };

} //namespace ext
} //namespace nanos




#endif //_FGCS_MALLEABLE_FREQ_ARCH_SCHED_HPP_

