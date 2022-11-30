#ifndef _FGCS_MALLEABLE_FREQ_SCHED_HPP_
#define _FGCS_MALLEABLE_FREQ_SCHED_HPP_

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

#ifdef NANOS_DEBUG_ENABLED
#define ONLY_DBG(f) f;				
#else
#define ONLY_DBG(f) ;				
#endif


// Nyapas para hacer pruebas rapidas. Poner esto bien cuando corresponda
#if defined(KELVIN) || defined(kelvin)
#define MAX_FREQ 4000000
#define MIN_FREQ 1000000
#elif defined(MAKALU) || defined(makalu)
#define MAX_FREQ 1700000
#define MIN_FREQ 1000000
#else
#error "machine not defined"
#endif

namespace nanos {    
    namespace ext {

	float getTimeStamp();
	
	typedef enum {
	    RUNNING,
	    IDLE,
	    EXPROPRIATED,
	    UNDEF
	} ThreadStatus;

	enum AppStatus {
	    START,
	    STOP
	};

	
	class FgcsMalleableFreq : public SchedulePolicy {
	private:
	    struct TeamData : public ScheduleTeamData {
		WDPriorityQueue<> *_readyQueue;
		TeamData () : ScheduleTeamData(), _readyQueue(NULL) {_readyQueue = NEW WDPriorityQueue<>;}
		~TeamData () { delete _readyQueue; }
	    };

	    struct {
		int numCores;
		int startBinding;
		
		//---- LOCK: wake/sleep cores
		Lock         waking_lock;
		sem_t        *sems;

		//---- LOCK: threadId. AppStatus. nActiveThreads
		Lock threadStatus_lck;
		
		ThreadStatus   *threadStatus;
		AppStatus      appStatus;
		int            nActiveThreads;

		//---- LOCK:
		Lock   changeInProgress;
		//Lock   futureChange;
		bool   futureChange;
		int    runningFreq_idx;

		FILE*  logFile;
	    } _expControl;
	    
	    //   1000 -> 2000
	    const int _NFREQS;
	    double _powers[11];
	    const double _powerIdle;
	    double _maxBudget;
	    
	public:
	    FgcsMalleableFreq();
	    virtual ~FgcsMalleableFreq();
	    

	private:
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


	    virtual WD * atWakeUp      ( BaseThread *thread, WD &wd );


	    /**********
	     ** EXPS **
	     *********/
	private:
	    inline void sleepThread( int threadId );
	    inline void wakeUp_threads( int maxThreads );

	    void calculateNewFreq();
	    
	    inline void changeFreq(int idx, long long freq) const;
	    inline void increaseFreq(int idx, long long freq) const;
	    inline void decreaseFreq(int idx, long long freq) const;

	    

	    inline long long getFreqFromIdx_Hz(int idx) const { return (10+idx)*100000; }
	    inline long double estimatePower_idx(unsigned nIdle,
						 unsigned nRunning, unsigned rFreqIdx) const {
		return (nIdle*_powerIdle + (nRunning*_powers[rFreqIdx]));    }

	    inline float getCurrentPowerConsumption() {
		return (((_expControl.numCores-_expControl.nActiveThreads)*_powerIdle) +
			 (_expControl.nActiveThreads*_powers[_expControl.runningFreq_idx]));
	    }
	    inline float getDesiredBudget() {
		return (((_expControl.numCores-_expControl.nActiveThreads)*_powerIdle) +
			 (_expControl.nActiveThreads*_powers[_NFREQS-1]));
	    }

	    inline float getDesiredBudget(int nThs) {
		return (((_expControl.numCores-nThs)*_powerIdle) +
			(nThs*_powers[_NFREQS-1]));
	    }

            inline int getMaxThreadsWithBudget(double maxBudget) {
              return ((maxBudget - 10.0*_powerIdle) / ( 2.5 - _powerIdle)); //ponemos 2.5 en vez de _powers[0]
            }
          
	    
		    
	};
	
	    
	class FgcsMalleableFreqPlugin : public Plugin {
	public:
	    FgcsMalleableFreqPlugin() : Plugin( "FGCS MalleableFreq Plugin",1 ) {}

	    virtual void config ( Config &cfg ){ }
	    virtual void init() {
		sys.setDefaultSchedulePolicy(NEW FgcsMalleableFreq());
	    }
	};

    }
}



#endif //_FGCS_MALLEABLE_FREQ_SCHED_HPP_










