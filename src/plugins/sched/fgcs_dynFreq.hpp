#ifndef _FGCS_DYNFREQ_HPP_
#define _FGCS_DYNFREQ_HPP_

#include "schedule.hpp"
#include "wddeque.hpp"
#include "plugin.hpp"
#include "system.hpp"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <semaphore.h>

#include "client.h"

/*@ @*/
#include <cpufreq.h>

#if !defined(VOLTA1) && !defined(volta1) && !defined(MAKALU) && !defined(makalu)
#error "ERROR: Esta politica solo está definida para makalu"
#endif

#ifdef NANOS_DEBUG_ENABLED
#define ONLY_DBG(f) f;
#else
#define ONLY_DBG(f) ;				
#endif


#define NFREQS 11
/*@-- @*/




namespace nanos {
namespace ext {

typedef enum {
  CORE_RUNNING,
  CORE_IDLE,
  CORE_UNDEF
} coreStatus_t;

typedef enum {
  APP_RUNNING,
  APP_STOPPING,
  APP_UNDEF
} appStatus_t;
    
	 

/**-------------------**/
/** ____BotLevCfg____ **/
/**-------------------**/
/* Configuration options for the BotLev scheduling policy */
struct BotLevCfg {
  static int   updateFreq;
  static int   numSpins;
  static int   hpFrom;
  static int   hpTo;
  static int   hpSingle;
  static int   steal;
  static int   maxBL;
  static int   strict;
  static int   taskNumber;
  NANOS_INSTRUMENT( static int   numCritical; )   //! The number of critical tasks (for instrumentation)
};


/**---------------------------**/
/** ____Experiment Config____ **/
/**---------------------------**/
struct ExperimentsGlobalCfg {
  static long double MaxBudget;
};

/**----------------------**/
/** ____BotLevDOData____ **/
/**----------------------**/
class BotLevDOData : public DOSchedulerData {
 public:
  typedef std::set<BotLevDOData *> predecessors_t;
 private:
  int               _taskNumber;    //! Keep task number just for debugging, remove?
  int               _botLevel;      //! Bottom Level Value
  bool              _isReady;       //! Is the task ready
  Lock              _lock;          //! Structure lock
  WD               *_wd;            //! Related WorkDescriptor
  predecessors_t    _predecessors;  //! List of BotLevDOData predecessors 
  short             _isCritical;    //! Is the task critical -- needed for reordering the ready queues

 public:
  BotLevDOData(int tnum, int blev);
  ~BotLevDOData();
  void reset ();

  int getTaskNumber() const;

  int getBotLevel() const;
  bool setBotLevel( int blev );

  bool getReady() const;
  void setReady();
 
  void setCriticality( short c );
  short getCriticality();

  WD* getWorkDescriptor() const;
  void setWorkDescriptor(WD *wd);

  std::set<BotLevDOData *> *getPredecessors();
  void addPredecessor(BotLevDOData *predDod);
};


class expData_t {
 public:
  Lock  mtxArch;
  int nRunning;
  //int nIdle; == nWorkers - nRunning
  const int nWorkers;  

  int runningFreq;
  const int idleFreq;

  coreStatus_t *coreStatus; //CORE_RUNNING, CORE_IDLE, CORE_UNDEF
  sem_t        *sems;

  long double power[11];
  long double powerIdle;

  appStatus_t appStatus;
  expData_t(int nCores) :
      nRunning(nCores), nWorkers(nCores), runningFreq(0), idleFreq(0),
      appStatus(APP_RUNNING)
  {
    coreStatus = new coreStatus_t[nCores];
    sems       = new sem_t[nCores];
    
    for(int i=0; i<nCores; i++){
      coreStatus[i]=CORE_RUNNING;
      sem_init(&sems[i], 0, 1); //init threads unblocked.
    }
    
    power[0] = 5.1241; //1000
    power[1] = 5.5746; //1100
    power[2] = 5.9394; //1200
    power[3] = 5.9938; //1300
    power[4] = 6.3564; //1400
    power[5] = 6.7961; //1500
    power[6] = 6.9460; //1600
    power[7] = 7.4403; //1700
    power[8] = 8.2584; //1800
    power[9] = 8.5159; //1900
    power[10] = 9.4979; //2000

    powerIdle = 1.474423832;
  }

    ~expData_t(){
	for(int i=0; i<nWorkers; i++) sem_destroy(&sems[i]);
	delete[] sems;
	delete[] coreStatus;
    }
};


/**----------------**/
/** ____BotLev____ **/
/**----------------**/
class BotLev : public SchedulePolicy
{
 public:
  using SchedulePolicy::queue;
  typedef std::stack<BotLevDOData *>   bot_lev_dos_t;
  typedef std::set<std::pair< unsigned int, DependableObject * > > DepObjVector; /**< Type vector of successors  */
	    
 private:
  bot_lev_dos_t     _blStack;       //! tasks added, pending having their bottom level updated
  Lock              _stackLock;
  DepObjVector      _topSuccesors;  //! Successors of the last maxPriority task
  Lock              _botLevLock;    //! Lock used for topSuccessors and currMax
  int               _currMax;       //! The priority of the last critical task
  int               _maxBotLev;     //! The maximum priority of the tdg

  expData_t         _expData;

	    
  struct TeamData : public ScheduleTeamData
  {
    // queues of ready tasks to be executed
    WDPriorityQueue<> *_readyQueues;

#ifdef NANOS_DEBUG_ENABLED
    Lock* _lockfout;
    FILE* _fout;
#endif
		
    TeamData () : ScheduleTeamData() {
      _readyQueues = NEW WDPriorityQueue<>[3];

#ifdef NANOS_DEBUG_ENABLED      
      _lockfout = new Lock();
      _fout = fopen("MUESTRA.m", "w");
      fprintf(_fout, "data0=[];\ndata1=[];\nfreq=[];\n");
      fflush(_fout);
#endif
    }
    virtual ~TeamData () {
      delete[] _readyQueues;

#ifdef NANOS_DEBUG_ENABLED
      fclose(_fout);
      delete _lockfout;
#endif
    }
  };

  
  /* disable copy and assigment */
  explicit BotLev ( const BotLev & );
  const BotLev & operator= ( const BotLev & );

  typedef enum {
    NANOS_SCHED_BLEV_ZERO,
    NANOS_SCHED_BLEV_ATCREATE,
    NANOS_SCHED_BLEV_ATSUBMIT,
    NANOS_SCHED_BLEV_ATIDLE
  } sched_bf_event_value;

 public:
  BotLev();
  virtual ~BotLev();

  virtual size_t getTeamDataSize () const;
  virtual size_t getThreadDataSize () const;
  virtual ScheduleTeamData * createTeamData ();
  virtual ScheduleThreadData * createThreadData ();
	    
  struct WDData : public ScheduleWDData
  {
    int _criticality;

    void setCriticality( int cr )            { _criticality = cr; }
    int getCriticality ( )                   { return _criticality; }
    WDData () : _criticality( 0 )            { }
    virtual ~WDData()                        { }
  };

  virtual void queue ( BaseThread *thread, WD &wd );
  void updateBottomLevels( BotLevDOData *dodata, int botLev );
  void atCreate ( DependableObject &depObj );
  virtual WD * atSubmit ( BaseThread *thread, WD &newWD );
  virtual WD *atIdle( BaseThread *thread, int numSteal );

  virtual void atShutdown( );
  virtual WD *atBeforeExit( BaseThread *thread, WD &current, bool schedule );


  /**********
   ** EXPS **
   *********/
  inline long long getFreqFromIdx_Hz(int idx){ return (10+idx)*100000; }
  /** OJO: asumimos que el mutex está cogido **/
  void calculateNewFreqs();
  void modifyFreqs();
  void wakeUpWorker();

  inline long double estimatePower_idx(unsigned nIdle,
                                       unsigned nRunning, unsigned rFreqIdx){
      return (nIdle*_expData.powerIdle + (nRunning*_expData.power[rFreqIdx]));
  }
    	    
};

class fgcsDynFreqPlugin : public Plugin
{
 public:
  fgcsDynFreqPlugin();
  virtual void config ( Config &config_ );
  virtual void init();
};
} //namespace ext
} //namespace nanos


#endif //_FGCS_DYNFREQ_HPP_

