#include "schedule.hpp"
#include "wddeque.hpp"
#include "plugin.hpp"
#include "system.hpp"
#include "config.hpp"

#include "fgcsMalleable_sched.hpp"
#include "smpplugin_decl.hpp"
#include <cmath>

namespace nanos {
namespace ext {
    float getTimeStamp()
    {
	static double gtod_ref_time_sec = 0.0;

	struct timeval tv;
	gettimeofday(&tv, NULL);
	    
	// If this is the first invocation of through dclock(), then initialize the
	// "reference time" global variable to the seconds field of the tv struct.
	if (gtod_ref_time_sec == 0.0)
	    gtod_ref_time_sec = (double) tv.tv_sec;

	// Normalize the seconds field of the tv struct so that it is relative to the
	// "reference time" that was recorded during the first invocation of dclock().
	const double norm_sec = (double) tv.tv_sec - gtod_ref_time_sec;

	// Compute the number of seconds since the reference time.
	const double t = norm_sec + tv.tv_usec * 1.0e-6;

	return (float) t;
    }
    
FgcsMalleable::FgcsMalleable()  : SchedulePolicy("FGCS malleable"),
				  _NFREQS(10), _powerIdle(1.474423832),
				  _maxBudget(500.0)
{
    /** A) ABOUT NUMBER OF WORKERS, ETC **/
    char *txt=getenv("DISABLE_BINDING");
    if(txt==NULL){
	fprintf(stderr, "ERROR, DISABLE_BINDING is not set! (chapuza para solucionar lo de nanox)\n");
	abort();
    }
	  
    txt=getenv("START_BINDING");
    if(txt==NULL){
	fprintf(stderr, "ERROR, START_BINDING is not defined (chapuza para solucionar lo de nanox)\n");
	abort();
    }
	  
    /** 0) LOG FILE **/
    _expControl.logFile =
	( getenv("NANOS_LOGFILE") == NULL) ?
	fopen("nanos_log.txt","w+") :
	fopen(getenv("NANOS_LOGFILE"),"w+");


    struct timeval tv;
    gettimeofday(&tv, NULL);
    const unsigned long long refTime  = tv.tv_sec*1e6 + tv.tv_usec;
	  
    fprintf(_expControl.logFile,
	    "refTime  = %Lu\n"		   \
	    "currFreq = [];\n"		   \
	    "nActiveThreads = [];\n"	   \
	    "maxBudget = [];\n"		   \
	    "\n"			   \
	    "\n",
	    refTime);
	  

  /** A) INIT CONTROL STRUCTURES & VARIABLES **/
  _expControl.numCores         = sys.getNumThreads();
  _expControl.threadStatus     = new ThreadStatus[_expControl.numCores];
  _expControl.sems             = new sem_t[_expControl.numCores];
  _expControl.coreOccupation   = new int[_expControl.numCores];
  _expControl.threadAssignment = new int[_expControl.numCores];
  _expControl.nActiveThreads   = 0;
  _expControl.futureChange     = false;

  _expControl.currMask.clear();

  /** B) INIT THREAD STATUS. BLOCK THRS. ETC ... **/
  for(int i=0; i<_expControl.numCores; i++){
    _expControl.threadStatus[i] = IDLE;
    _expControl.coreOccupation[i] = -1;
    _expControl.threadAssignment[i] = -1;
    sem_init(&_expControl.sems[i], 0, 0); //init threads blocked.
  }


  /** C) ABOUT FREQUENCIES **/
  _expControl.runningFreq_idx = _NFREQS-1;	
  _powers[0]=2.97; //1000
  _powers[1]=3.09;
  _powers[2]=3.22;
  _powers[3]=3.35;
  _powers[4]=3.43;
  _powers[5]=3.56;
  _powers[6]=3.79;
  _powers[7]=3.97;
  _powers[8]=4.24;
  _powers[9]=4.45; 
  _powers[10]=6.07; //2000 -> not used


  /** D) APP START **/	  
  SS_initLibrary();

  /** D.1) Affinity  **/
  dynamic_bitset bt;
  SS_startExecution(&bt);
  if(bt.count() == 0) SS_askForNewCpuMask(&bt);
  ONLY_DBG( fprintf(stderr, "START mask: %s\n", bt.to_string().c_str()) );

  increaseAffinity(bt);
  //wakeUp_sleepingWorkers(_expControl.currMask.count());
  wakeUp_sleepingWorkers(1);
  

  /** D.2) Frequency **/
  _maxBudget = SS_requestNewBudget (0, getDesiredBudget());
  calculateNewFreq();

  
  _expControl.appStatus=START;

  
  fprintf(_expControl.logFile,
	  "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
  fprintf(_expControl.logFile,
	  "nAssignedCores = [nAssignedCores; [ %f %d ]];\n", getTimeStamp(), _expControl.currMask.count());
  fprintf(_expControl.logFile,
	  "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), _maxBudget);
}

FgcsMalleable::~FgcsMalleable() 
{
  for(int i=0; i<_expControl.numCores; i++)
    sem_destroy(&_expControl.sems[i]);

  delete[] _expControl.sems;
  delete[] _expControl.threadStatus;
  delete[] _expControl.coreOccupation;
  delete[] _expControl.threadAssignment;

  fclose (_expControl.logFile); 
}
	
void FgcsMalleable::queue ( BaseThread *thread, WD &wd )
{
  //puts("QUEUE 1");
	    
  TeamData &tdata = (TeamData &) *thread->getTeam()->getScheduleData();
  //wd.setPriority(0);
  tdata._readyQueue->push_back( &wd );
}
    
void FgcsMalleable::queue ( BaseThread ** threads, WD ** wds, size_t numElems )
{
  //puts("QUEUE 2");
	    
  ThreadTeam* team = threads[0]->getTeam();
  TeamData &tdata = (TeamData &) *team->getScheduleData();
  tdata._readyQueue->push_back( wds, numElems );
}
    
    
void FgcsMalleable::atSuccessor   ( DependableObject &successor, DependableObject &predecessor )
{
  //puts("AT SUCCESSOR");	    
  return;
}
    
WD *FgcsMalleable::atSubmit ( BaseThread *thread, WD &newWD )
{
  //puts( "AT SUBMIT" );	    
  queue( thread, newWD );
  return 0;
}
    
	    
WD * FgcsMalleable::atIdle ( BaseThread *thread, int numSteal )
{
  if(_expControl.appStatus == STOP) return NULL;
  int threadId = thread->getCpuId();


  //A) Check if we need to release the current Core [not Core0])
  if(threadId != 0){
      while(SS_needToFreeCore()>0){
	  freeWorkerThread(thread);
      }
  }

  
  if( _expControl.threadStatus[threadId] == WAKING_UP)
      waitForNewMask(threadId);
    
  //Maybe someone has blocked us
  sem_wait(&_expControl.sems[threadId]);  
  if(_expControl.appStatus == STOP) return NULL;
  
  TeamData &tdata = (TeamData &) *thread->getTeam()->getScheduleData();


  /* Before waking up more threads, we update the budget. Maybe the budget increases
     and we can wake up more threads
  */
  const int aux = _expControl.nActiveThreads + std::max(static_cast<int>(tdata._readyQueue->size()-1), 0);
  const int maxThreadsPossible=std::min( _expControl.currMask.count(), aux);
	  
  float newBudget=SS_requestNewBudget(getCurrentPowerConsumption(), getDesiredBudget(maxThreadsPossible));
  if( newBudget > 0 && newBudget != _maxBudget){
      _maxBudget = newBudget;

      /*@ LOG @*/
      ONLY_DBG( fprintf(stderr, "New Budget %.2f\n", newBudget) );
      fprintf(_expControl.logFile,
	      "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), newBudget);
      /*@-----@*/
	      
      calculateNewFreq();
  }



	  

  
  //Check if we have to wake up more threads.
  /*** Despertamos threads en este punto para evitar overhead:
   **    El propio thread encargado de despertar los hilos es el encargado
   **    de cambiar la frecuencia (y por tanto esperar en los mutex).
   **    Para evitar que una tarea esté bloqueadea en esos mutex, despertamos a 
   **    los Workers antes de coger una tarea nueva!
   ***/
  if(tdata._readyQueue->size() > 1 && _expControl.nActiveThreads < _expControl.currMask.count()){      
    wakeUp_sleepingWorkers(std::min(static_cast<int>(tdata._readyQueue->size()-1), getMaxThreadsWithBudget(_maxBudget) - _expControl.nActiveThreads));

    newBudget=SS_requestNewBudget(getCurrentPowerConsumption(), getDesiredBudget());
    if(newBudget != _maxBudget){
	_maxBudget = newBudget;

	/*@ LOG @*/
	ONLY_DBG( fprintf(stderr, "New Budget %.2f\n", newBudget) );
	fprintf(_expControl.logFile,
		"maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), newBudget);
	/*@-----@*/
    }

    calculateNewFreq();
  }
  

  WD* wd = tdata._readyQueue->pop_front( thread );

  int spins=60000;//10;
  while(wd==NULL && --spins){
    wd = tdata._readyQueue->pop_front( thread );
  }



  if( wd==NULL ){
    if( threadId==0 ){
      sem_post(&_expControl.sems[threadId]);
    } else {
	if(_expControl.appStatus == START){
	    ONLY_DBG( fprintf(stderr, "Sleeping thread %d [nRunning=%d]\n", threadId, _expControl.nActiveThreads) );
	    sleepThread(threadId);
	}
    }
  }
  
  return wd;          
}
    
WD * FgcsMalleable::atPrefetch ( BaseThread *thread, WD &current )
{
  //puts("AT PREFETCH");
  return atIdle(thread,false);
}
    
WD * FgcsMalleable::atBeforeExit ( BaseThread *thread, WD &current, bool schedule )
{
  //puts("AT BEFORE EXIT");
  int coreId=thread->getCpuId();
  sem_post(&_expControl.sems[coreId]);
  return NULL;
}

void FgcsMalleable::atShutdown()
{
//  puts("AT SHUTDOWN");

    //_expControl.threadStatus_lck.acquire();
    _expControl.appStatus=STOP;
    //_expControl.threadStatus_lck.release();
    
    SS_finishLibrary();

    for(int i=0; i<_expControl.numCores; i++){
	if(_expControl.threadStatus[i] == SLEEPING || _expControl.threadStatus[i]==IDLE){
	    sem_post(&_expControl.sems[i]);
	}
    }
}

WD* FgcsMalleable::atWakeUp ( BaseThread *thread, WD &wd )
{
  // for(int i=0; i<_expControl.numCores; i++)
  //   sem_post(&_expControl.sems[i]);

    return NULL;
}




//==============================================================================

void FgcsMalleable::waitForNewMask(int threadId)
{
  if(!_expControl.waking_lock.tryAcquire())
    return;

  _expControl.threadStatus[threadId] = WAKING_UP;


  dynamic_bitset bt;
  while(_expControl.threadStatus[threadId] == WAKING_UP &&
        _expControl.appStatus==START) {

    bt.clear();

    // wait for new mask
    do{
      SS_askForNewCpuMask(&bt);    
    } while(_expControl.appStatus==START && bt.count()==0);

    increaseAffinity(bt);

    fprintf(_expControl.logFile,
	    "nAssignedCores = [nAssignedCores; [ %f %d ]];\n", getTimeStamp(), _expControl.currMask.count());
  }

  _expControl.waking_lock.release();
}


void FgcsMalleable::increaseAffinity(const dynamic_bitset& bt)
{
  if(bt.count()==0) return;
  
  ONLY_DBG( assert( _expControl.currMask.count() < bt.count()) );


  LockBlock l(_expControl.threadStatus_lck); // <<====
  
  const unsigned mask = ( bt.count() == _expControl.numCores)
      ? (IDLE | WAKING_UP)
      : IDLE;
  
  int coreIdx = -1;
  int threadId = 0;
  while((coreIdx = bt.find_next(coreIdx)) >=0){
    //this core is already used by a thread
    if(_expControl.coreOccupation[coreIdx]!=-1) continue;
    
    //find the first thread in mask status (sleeping / wakingUP)
    for(; threadId<_expControl.numCores; threadId++){
      if( (_expControl.threadStatus[threadId] & mask ) !=0 )
        break;
    }
    
    ONLY_DBG( assert(_expControl.threadStatus[threadId] != RUNNING) );
    ONLY_DBG( fprintf(stderr,"Waking up thread %d\n", threadId) );
	    
    //A) set affinity
    ((SMPPlugin*)(sys.getSMPPlugin()))->setAffinity(threadId, coreIdx);
        
    //B) Register association coreId <-> threadId
    _expControl.coreOccupation[coreIdx]     = threadId;
    _expControl.threadAssignment[threadId]  = coreIdx;

    //C) Mark worker as sleeping
    _expControl.threadStatus[threadId] = SLEEPING;

    //D) If sleeping -> set minimum freq
    changeFreq(threadId, MIN_FREQ);
  }


  //no previous threads running. Need to start an additionally thread to
  //wait for the new mask
  if(_expControl.nActiveThreads == 0 && bt.count()!=_expControl.numCores)
    putThread_toWait();
  
  _expControl.currMask=bt;
}



void FgcsMalleable::putThread_toWait() {
  ONLY_DBG( assert( _expControl.nActiveThreads < _expControl.numCores));

  int idx=0;

  //A) There is a thread sleeping -> set it to wait
  if( _expControl.nActiveThreads < _expControl.currMask.count()){
    for(; idx < _expControl.numCores && _expControl.threadStatus[idx] != SLEEPING; ++idx);
  }
  else {   //B) All assigned cores are in use
    for(; idx < _expControl.numCores && _expControl.threadStatus[idx] != IDLE; ++idx);
    ((SMPPlugin*)sys.getSMPPlugin())->setAffinity(idx, _expControl.threadAssignment[0]);
  }
  
  ONLY_DBG( assert(_expControl.threadStatus[idx] != RUNNING) );
  ONLY_DBG( fprintf(stderr, "WAITING thread [%d]\n", idx) );


  _expControl.threadStatus[idx] = WAKING_UP;
  sem_post(&_expControl.sems[idx]);
}




void FgcsMalleable::wakeUp_sleepingWorkers(int numWorkers)
{ 
  {
    LockBlock l(_expControl.threadStatus_lck); // <<====
    for(int threadId=0; threadId<_expControl.numCores && numWorkers>0; ++threadId){
      if( _expControl.threadStatus[threadId] == SLEEPING){
        ONLY_DBG( fprintf(stderr,"Waking up thread (sleep) %d\n", threadId) );

        _expControl.threadStatus[threadId] = RUNNING;
        ++_expControl.nActiveThreads;
        sem_post(&_expControl.sems[threadId]);

        --numWorkers;
      }
    }
  } //lockblock

  /*@ LOG @*/
  fprintf(_expControl.logFile,
	  "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
  /*@-----@*/
    
  ONLY_DBG( fprintf(stderr, "After waking up: %d threads are running\n", _expControl.nActiveThreads) );

}


void FgcsMalleable::sleepThread(int threadId)
{
    if(_expControl.threadStatus[threadId] == WAKING_UP)
	return;
    
    {
	LockBlock l(_expControl.threadStatus_lck);

	ONLY_DBG( assert( _expControl.threadStatus[threadId] == RUNNING ) );
	_expControl.threadStatus[threadId] = SLEEPING;
	--_expControl.nActiveThreads;
	changeFreq(threadId, MIN_FREQ);

	/*@ LOG @*/
	fprintf(_expControl.logFile,
		"nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
	/*@-----@*/
    }

    /**
     * Before doing anything, check if we have to modify the budget 
     **/
    float newBudget=SS_requestNewBudget(getCurrentPowerConsumption(), getDesiredBudget());
    if(newBudget > 0 && newBudget != _maxBudget){
	_maxBudget = newBudget;

	/*@ LOG @*/
	ONLY_DBG( fprintf(stderr, "New Budget %.2f\n", newBudget) );
	fprintf(_expControl.logFile,
		"maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), newBudget);
	/*@-----@*/	
    }
    calculateNewFreq();
}


//ASUMIMOS QUE EL LOCK YA ESTA COGIDO
void FgcsMalleable::freeSleepingWorker()
{

  int coreId=-1;
  int threadId=0;
  for(; threadId < _expControl.numCores && _expControl.threadStatus[threadId]!=SLEEPING; ++threadId);
  ONLY_DBG (  assert (threadId < _expControl.numCores) );
  
  coreId = _expControl.threadAssignment[threadId];
  _expControl.threadStatus[threadId] = IDLE;
  
  _expControl.threadAssignment[threadId] = -1;
  _expControl.coreOccupation[coreId]    = -1;
  _expControl.currMask.reset(coreId);

  ONLY_DBG (fprintf (stderr, "Free sleeping core %d\n", coreId));
      

  dynamic_bitset newMask;
  SS_informFreeCore(coreId, &newMask);
  
}



void FgcsMalleable::freeWorkerThread(BaseThread *thread)
{

   
  _expControl.threadStatus_lck.acquire();
  {    
      if(_expControl.nActiveThreads < _expControl.currMask.count()){
	  //#ERROR: quizas hay que liberar el core que esta WAKING_UP
	  freeSleepingWorker();
      } else { 
	  const int threadId = thread->getCpuId();
  
	  sem_wait(&_expControl.sems[threadId]);
	  if(_expControl.appStatus == STOP) return;
	  _expControl.threadStatus[threadId] = IDLE;

	  int coreId = _expControl.threadAssignment[threadId];
	  ONLY_DBG( assert( coreId == ((SMPThread*)thread)->getRealCpuId()) );
	    
	  _expControl.threadAssignment[threadId] = -1;
	  _expControl.coreOccupation[coreId]    = -1;
	  _expControl.currMask.reset(coreId);
  
	  --_expControl.nActiveThreads;

	  ONLY_DBG ( fprintf (stderr, "Free workerThread id: %d [core:%d]", threadId, coreId) );

	  dynamic_bitset newMask;
	  SS_informFreeCore(coreId, &newMask);
      }
  }
  _expControl.threadStatus_lck.release();


  fprintf(_expControl.logFile,
	  "nAssignedCores = [nAssignedCores; [ %f %d ]];\n", getTimeStamp(), _expControl.currMask.count());
  fprintf(_expControl.logFile,
	  "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);

  /**
   * Before doing anything, check if we have to modify the budget 
   **/
  float newBudget=SS_requestNewBudget(getCurrentPowerConsumption(), getDesiredBudget());
  if(newBudget > 0 && newBudget != _maxBudget){
      _maxBudget = newBudget;

      /*@ LOG @*/
      ONLY_DBG( fprintf(stderr, "New Budget %.2f\n", newBudget) );
      fprintf(_expControl.logFile,
	      "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), newBudget);
      /*@-----@*/	
  }
  calculateNewFreq();

} //freeWorkerThread



/**** FREQ *****/
void FgcsMalleable::calculateNewFreq()
{
    if(_expControl.changeInProgress.tryAcquire()){    //no one is modifying the freq
    }
    else{
	_expControl.futureChange=true;
	return;
    }
	  
    const unsigned int nIdle    = _expControl.currMask.count() - _expControl.nActiveThreads;
    const unsigned int nRunning = _expControl.nActiveThreads;
    static const long double eps2= 0.8; //1.7; //0.001;//Error en las centesimas //2*std::numeric_limits<double>::epsilon();
	    
    const long double maxBudget = _maxBudget;
	    
    _expControl.futureChange=false;
	    
	  
    int fIdx = 0;
    for(int i=_NFREQS-1; i>=0; --i){
	const long double estimation=estimatePower_idx(nIdle, nRunning, i);
	const long double absDiff=std::abs(maxBudget-estimation);

	ONLY_DBG( fprintf(stderr, "maxBudget: %Lf\n", maxBudget) );
	ONLY_DBG( fprintf(stderr, "estimation: %Lf\n", estimation) );
	ONLY_DBG( fprintf(stderr, "diff: %Lf\n", absDiff); );
	ONLY_DBG( fprintf(stderr, "eps2: %Lf\n", eps2));
	ONLY_DBG(	fprintf(stderr, "nActive: %u\n", nRunning));
	      
	if( estimation<= maxBudget || absDiff <=eps2){
	    fIdx = i;
	    _expControl.runningFreq_idx = i;
	    break;
	}
    }
    ONLY_DBG(fprintf(stderr, "\n\n"));
    
    
    ONLY_DBG( fprintf(stderr, "New freq %Lu for %dthreads\n", getFreqFromIdx_Hz(fIdx), _expControl.nActiveThreads) );
    
    /*@ LOG @*/
    fprintf(_expControl.logFile,
    	    "currFreq = [currFreq; [ %f %Lu ]];\n", getTimeStamp(), getFreqFromIdx_Hz(fIdx));
    /*@-----@*/
	
    for(int i=0; i<_expControl.numCores; ++i){
	if(_expControl.threadStatus[i] == RUNNING){
	    ONLY_DBG( fprintf(stderr, "\t Thread %d -> freq %Lu\n", i, getFreqFromIdx_Hz(fIdx)) );
	    changeFreq(i, getFreqFromIdx_Hz(fIdx));
	}
    }
	    
    _expControl.changeInProgress.release();
	  
    if(_expControl.futureChange) //maybe someone change the arch while we were modifing the freq
	calculateNewFreq();	  
}


}  //namespace ext
}//namespace nanos


DECLARE_PLUGIN("fgcsMalleable",nanos::ext::FgcsMalleablePlugin);
