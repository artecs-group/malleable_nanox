#include "schedule.hpp"
#include "wddeque.hpp"
#include "plugin.hpp"
#include "system.hpp"
#include "config.hpp"

#include "fgcsMalleableFreqArch_sched.hpp"
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
    
    FgcsMalleableFreqArch::FgcsMalleableFreqArch()  : SchedulePolicy("FGCS malleable Freq & Arch"),
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
	( getenv("NANOS_LOGFILE") == NULL ) ?
	fopen("nanos_log.txt","w+") :
	fopen(getenv("NANOS_LOGFILE"),"w+");

      
      struct timeval tv;
      gettimeofday(&tv, NULL);
      const unsigned long long refTime  = tv.tv_sec*1e6 + tv.tv_usec;
      
      fprintf(_expControl.logFile,
	      "refTime  = %Lu\n"		   \
	      "currFreq = [];\n"		   \
	      "nActiveThreads = [];\n"		   \
	      "maxBudget = [];\n"		   \
	      "\n"				   \
	      "\n",
	      refTime);
      fflush(_expControl.logFile);

      /** A) INIT CONTROL STRUCTURES & VARIABLES **/
      _expControl.numCores         = sys.getNumThreads() -1 ; /* -1 IS IMPORTANT */
      _expControl.sems             = new sem_t[_expControl.numCores];

      _expControl.threadStatus     = new threadStatus_t[_expControl.numCores + 1];
      _expControl.to_untie         = new bool[_expControl.numCores];
      _expControl.coreOccupation   = new int[_expControl.numCores];
      _expControl.threadAssignment = new int[_expControl.numCores];

      _expControl.nActiveThreads   = 0;
      _expControl.futureChange     = false;
      _expControl.runningFreq_idx  = 0;
      
      _expControl.currMask.clear();

      /** B) INIT THREAD STATUS. BLOCK THRS. ETC ... **/
      for(int i=0; i<_expControl.numCores; i++){
	ONLY_DBG (fprintf (stderr, ">> Setting status [%d] to %s in line %d\n", i, "NOT_TIED", __LINE__) );
	_expControl.threadStatus[i] = NOT_TIED;
	_expControl.to_untie[i] = false;
	
	_expControl.coreOccupation[i] = -1;
	_expControl.threadAssignment[i] = -1;
	sem_init(&_expControl.sems[i], 0, 0); //init threads blocked.
      }


      /** C) ABOUT FREQUENCIES **/
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
      if(bt.count() == 0)  SS_askForNewCpuMask(&bt);
      ONLY_DBG( fprintf(stderr, "START mask: %s\n", bt.to_string().c_str()) );


      //-----------------
      //1) Affinity 
      increaseAffinity(bt);
      initializeWaitingThread();
      wakeUp_sleepingWorkers(1); //only worker 0 (master)
  
      
      //2) Frequency 
      _maxBudget = SS_requestNewBudget (0, getDesiredBudget());
      adjustFrequency();

 
      fprintf(_expControl.logFile,
	      "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
      fprintf(_expControl.logFile,
	      "nAssignedCores = [nAssignedCores; [ %f %d ]];\n", getTimeStamp(), _expControl.currMask.count());
      fprintf(_expControl.logFile,
	      "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), _maxBudget);


      _expControl.appStatus=START;
    }

    FgcsMalleableFreqArch::~FgcsMalleableFreqArch() 
    {
      for(int i=0; i<_expControl.numCores; i++)
	sem_destroy(&_expControl.sems[i]);

      delete[] _expControl.sems;
      delete[] _expControl.threadStatus;
      delete[] _expControl.to_untie;
      delete[] _expControl.coreOccupation;
      delete[] _expControl.threadAssignment;

      fclose (_expControl.logFile); 
    }

    
    void FgcsMalleableFreqArch::queue ( BaseThread *thread, WD &wd )
    {
      TeamData &tdata = (TeamData &) *thread->getTeam()->getScheduleData();
      tdata._readyQueue->push_back( &wd );
    }
    
    void FgcsMalleableFreqArch::queue ( BaseThread ** threads, WD ** wds, size_t numElems )
    {
      ThreadTeam* team = threads[0]->getTeam();
      TeamData &tdata = (TeamData &) *team->getScheduleData();
      tdata._readyQueue->push_back( wds, numElems );
    }
    
    void FgcsMalleableFreqArch::atSuccessor   ( DependableObject &successor, DependableObject &predecessor )
    {
      return;
    }
    
    WD *FgcsMalleableFreqArch::atSubmit ( BaseThread *thread, WD &newWD )
    {
      queue( thread, newWD );
      return 0;
    }
    
    WD * FgcsMalleableFreqArch::atPrefetch ( BaseThread *thread, WD &current )
    {
      return atIdle(thread,false);
    }
    
    WD * FgcsMalleableFreqArch::atBeforeExit ( BaseThread *thread, WD &current, bool schedule )
    {
      int threadId=thread->getCpuId();
      sem_post(&_expControl.sems[threadId]);
      return NULL;
    }

    void FgcsMalleableFreqArch::atShutdown()
    {
      _expControl.appStatus=STOP;

      SS_finishLibrary();
      const int mask = SLEEPING | NOT_TIED | WAKING_UP;
      for(int i=0; i<_expControl.numCores; i++){
	if((_expControl.threadStatus[i] & mask) !=0){
	  sem_post(&_expControl.sems[i]);
	}
      }
    }



    WD * FgcsMalleableFreqArch::atIdle ( BaseThread *thread, int numSteal )
    {
      if(_expControl.appStatus == STOP) return NULL;
      int threadId = thread->getCpuId();

      
      if( _expControl.threadStatus[threadId] == WAKING_UP){
	  waitForNewMask(threadId);
	  return NULL;
      }


      //Maybe someone has blocked us
      sem_wait(&_expControl.sems[threadId]);   /// <<-------------- BLOCK HERE
      if(_expControl.appStatus == STOP) return NULL;      
      

      if(_expControl.to_untie[threadId]){
	  untie_myself (threadId);
	  return NULL; //short-circuit to block earlier
      }
      
	 
//       //A) Check if we need to release a core
//       while(threadId != 0  && SS_needToFreeCore()>0 && untie_workerThread(thread) ){	
// /** TODO: RECALCULAR FRECUENCIAS AQUI por si somos bloqueados, o esperar a algun
//     hilo activo que pase el semaforo y lo solicite?? */
//       }


      TeamData &tdata = (TeamData &) *thread->getTeam()->getScheduleData();
      
      /* Before waking up more threads, we update the budget. Maybe the budget increases
	 and we can wake up more threads */
      const int aux = _expControl.nActiveThreads + std::max(static_cast<int>(tdata._readyQueue->size()-1), 0);
      const int maxThreadsPossible=std::min( _expControl.currMask.count(), aux);
      
      
      float newBudget=SS_requestNewBudget(getCurrentPowerConsumption(), getDesiredBudget(maxThreadsPossible));
      if( newBudget > 0 && newBudget != _maxBudget){
	//@ LOG
	ONLY_DBG( fprintf(stderr, "New Budget %.2f\n", newBudget) );
	fprintf(_expControl.logFile, "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), newBudget);

	_maxBudget = newBudget;
	adjustFrequency();
      }
      
      
      //Check if we have to wake up more threads.
      /*** Despertamos threads en este punto para evitar overhead:
       **    El propio thread encargado de despertar los hilos es el encargado
       **    de cambiar la frecuencia (y por tanto esperar en los mutex).
       **    Para evitar que una tarea esté bloqueadea en esos mutex, despertamos a 
       **    los Workers antes de coger una tarea nueva!
       ***/
      if(tdata._readyQueue->size() > 1 && _expControl.nActiveThreads < _expControl.currMask.count()){
	//	if(_expControl.wakeSleeping_lck.tryAcquire()){
	  wakeUp_sleepingWorkers( std::min( static_cast<int>(tdata._readyQueue->size()-1),
					    getMaxThreadsWithBudget(_maxBudget) - _expControl.nActiveThreads));
	//   _expControl.wakeSleeping_lck.release();
	// }
	
	/* //TODO: Comprobar si esto es necesario o no:
	newBudget=SS_requestNewBudget(getCurrentPowerConsumption(), getDesiredBudget());
	if(newBudget != _maxBudget){
	//@ LOG
	  ONLY_DBG( fprintf(stderr, "New Budget %.2f\n", newBudget) );
	  fprintf(_expControl.logFile,  "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), newBudget);
	  
	  _maxBudget = newBudget;
	}
	*/
	adjustFrequency();
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
	    sleepThread(thread);
	    adjustFrequency();
	  }
	}
      }
  
      return wd;          
    }
    
//==============================================================================

    void FgcsMalleableFreqArch::waitForNewMask(int threadId)
    {
      if(!_expControl.waking_lock.tryAcquire()){
	abort();
      }
      
      ONLY_DBG (assert (_expControl.threadStatus[threadId] == WAKING_UP ));

      
      dynamic_bitset bt;
      while( _expControl.appStatus==START ){
	bool adjustFreq = false;
	//A) Check for new mask
	bt.clear(); 
	
	// wait for new mask
	SS_askForNewCpuMask(&bt);    
	
	if(bt.count()!= 0){
	  ONLY_DBG (puts ("CHANGING AFFINITY") );
	  
	  increaseAffinity(bt);
	  fprintf(_expControl.logFile, "nAssignedCores = [nAssignedCores; [ %f %d ]];\n", getTimeStamp(), _expControl.currMask.count());
	  adjustFreq = true;	  
	} else {
	  int nToFree = SS_needToFreeCore();
	  ONLY_DBG (fprintf(stderr, "RELEASING %d CORES\n", nToFree) );
	  untie_workerThread(nToFree);
	  if(nToFree > 0) adjustFreq = true;
	}
	
	
	// //Si es la misma mascara, quizas sea que tenemos que despertar a alguien
	// LockBlock l (_expControl.threadStatus_lck);
	// while(_expControl.nActiveThreads < _expControl.currMask.count() //there are sleeping threads
	// 	&&  )
	// {
	//   untie_sleepingWorker();  

      
	if(adjustFreq){
	  float newBudget=SS_requestNewBudget(getCurrentPowerConsumption(), getDesiredBudget());
	  if( newBudget > 0 && newBudget != _maxBudget){
	    //@ LOGfile
	    ONLY_DBG( fprintf(stderr, "New Budget %.2f\n", newBudget) );
	    fprintf(_expControl.logFile, "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), newBudget);
	  
	    _maxBudget = newBudget;
	    adjustFrequency();
	  }
	}
      }
      
      _expControl.waking_lock.release();
    }

      


    
    void FgcsMalleableFreqArch::increaseAffinity(const dynamic_bitset& bt)
    {
      if(bt.count()==0) return;
      if(bt.count() <= _expControl.currMask.count()) return;
      
      //ONLY_DBG( assert( _expControl.currMask.count() < bt.count()) );
      ONLY_DBG( assert (bt.count() <= _expControl.numCores) );


      LockBlock l(_expControl.threadStatus_lck); // <<====
  
      int coreIdx = -1;
      int threadId = 0;
      
      while((coreIdx = bt.find_next(coreIdx)) >=0){
	//this core is already being used by a thread
	if(_expControl.coreOccupation[coreIdx] != -1) continue;
    
	//find the first thread in mask status (sleeping / wakingUP)
	for(;
	    threadId<_expControl.numCores && _expControl.threadStatus[threadId] != NOT_TIED;
	    threadId++)
	{}
	

	ONLY_DBG( assert (threadId < _expControl.numCores) );
	ONLY_DBG( assert(_expControl.threadStatus[threadId] == NOT_TIED) );
	ONLY_DBG( fprintf(stderr,"Setting affinity of thread %d to %d\n", threadId, coreIdx) );
	
	//A) set affinity
	((SMPPlugin*)(sys.getSMPPlugin()))->setAffinity(threadId, coreIdx);
        
	//B) Register association coreId <-> threadId
	_expControl.coreOccupation[coreIdx]     = threadId;
	_expControl.threadAssignment[threadId]  = coreIdx;

	//C) Mark worker as sleeping
	ONLY_DBG (fprintf (stderr, ">> Setting status [%d] to %s in line %d\n", threadId, "SLEEPING", __LINE__) );
	_expControl.threadStatus[threadId] = SLEEPING;
	changeFreq (threadId, MIN_FREQ);
      }     
  
      _expControl.currMask=bt;
    }



    void FgcsMalleableFreqArch::initializeWaitingThread() {
      const int threadId = _expControl.numCores;
      

      ONLY_DBG (fprintf (stderr, ">> Setting status [%d] to %s in line %d\n", threadId, "WAKING_UP", __LINE__) );
      _expControl.threadStatus[threadId] = WAKING_UP;
      ((SMPPlugin*)sys.getSMPPlugin())->setAffinity(threadId, _expControl.threadAssignment[0]);

      
      ONLY_DBG( assert(_expControl.threadStatus[threadId] != RUNNING) );
      ONLY_DBG( fprintf(stderr, "WAITING thread [%d]\n", threadId) );
    }



    void FgcsMalleableFreqArch::wakeUp_sleepingWorkers(int numWorkers)
    {
      if(_expControl.currMask.count() == _expControl.nActiveThreads) return;
      
      {
	LockBlock l(_expControl.threadStatus_lck); // <<====
	if(_expControl.appStatus == STOP) return;
	for(int threadId=0; threadId<_expControl.numCores && numWorkers>0; ++threadId){
	  if( _expControl.threadStatus[threadId] == SLEEPING) {
	    ONLY_DBG( fprintf(stderr,"Waking up thread (sleep) %d\n", threadId) );

	    ONLY_DBG (fprintf (stderr, ">> Setting status [%d] to %s in line %d\n", threadId, "RUNNING", __LINE__) );
	    _expControl.threadStatus[threadId] = RUNNING;
	    ++_expControl.nActiveThreads;
	    sem_post(&_expControl.sems[threadId]);

	    --numWorkers;
	  }
	}
	/*@ LOG @*/
	fprintf(_expControl.logFile, "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);    
	ONLY_DBG( fprintf(stderr, "After waking up: %d threads are running\n", _expControl.nActiveThreads) );
      } //lockblock
    }


    void FgcsMalleableFreqArch::untie_workerThread (int nToFree)
    {

      LockBlock l (_expControl.threadStatus_lck);
      if(_expControl.appStatus == STOP) return;
      
      while(nToFree--){
	  if(!untie_sleepingWorker()){
	      order_untie();
	  }
      } 
      
      /*@ LOG @*/
      fprintf(_expControl.logFile, "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
      fprintf(_expControl.logFile, "nAssignedCores = [nAssignedCores; [ %f %d ]];\n", getTimeStamp(), _expControl.currMask.count());
    }

    

    bool FgcsMalleableFreqArch::untie_sleepingWorker(){
      int threadId = 0;
      for(; threadId<_expControl.numCores && _expControl.threadStatus[threadId]!=SLEEPING; ++threadId);

      if(threadId == _expControl.numCores) return false; // There are not sleeping threads
      
      
      int coreId = _expControl.threadAssignment[threadId];
    
      ONLY_DBG ( assert (_expControl.threadStatus[threadId] == SLEEPING) );

      ONLY_DBG (fprintf (stderr, ">> Setting status [%d] to %s in line %d\n", threadId, "NOT_TIED", __LINE__) );
      _expControl.threadStatus[threadId] = NOT_TIED;
      _expControl.coreOccupation[coreId] = -1;
      _expControl.threadAssignment[threadId] = -1;
      _expControl.currMask.reset (coreId);
      
      ONLY_DBG (fprintf (stderr, "Free sleeping core %d\n", coreId));      
      dynamic_bitset newMask;
      SS_informFreeCore(coreId, &newMask);

      return true;
    }
    

    void FgcsMalleableFreqArch::order_untie (){
      ONLY_DBG (assert (_expControl.threadStatus_lck.getState() == NANOS_LOCK_BUSY) );

      int threadId=1;
      for(; threadId<_expControl.numCores && _expControl.threadStatus[threadId] != RUNNING; ++threadId);

      ONLY_DBG (assert (threadId < _expControl.numCores) );
      ONLY_DBG (assert (_expControl.threadStatus[threadId] == RUNNING) );
      _expControl.to_untie[threadId] = true;
    }

    
    void FgcsMalleableFreqArch::untie_myself (int threadId){
      ONLY_DBG ( assert (_expControl.to_untie[threadId]) );
      
      LockBlock l (_expControl.threadStatus_lck);

      int coreId =  _expControl.threadAssignment[threadId];
      threadStatus_t prevStatus = _expControl.threadStatus[threadId];

      ONLY_DBG (fprintf (stderr, ">> Setting status [%d] to %s in line %d\n", threadId, "NOT_TIED", __LINE__) );
      _expControl.threadStatus[threadId] = NOT_TIED;
      _expControl.to_untie[threadId] = false;
      _expControl.coreOccupation[coreId] = -1;
      _expControl.threadAssignment[threadId] = -1;
      _expControl.currMask.reset(coreId);
      
      if(prevStatus == RUNNING) {
	--_expControl.nActiveThreads;
	sem_wait (&_expControl.sems[threadId]);
      }
      
      
      ONLY_DBG (fprintf (stderr, "Free running core %d\n", coreId));      
      dynamic_bitset newMask;
      SS_informFreeCore(coreId, &newMask);
    }




    void FgcsMalleableFreqArch::sleepThread (BaseThread *thread){
      const int threadId = thread->getCpuId();

      if(_expControl.threadStatus[threadId] == SLEEPING) return;    
      LockBlock l(_expControl.threadStatus_lck);

      if(_expControl.to_untie[threadId]){
	untie_myself (threadId);
	return;
      } else {	
	ONLY_DBG( assert( _expControl.threadStatus[threadId] == RUNNING) );
      
      
	ONLY_DBG (fprintf (stderr, ">> Setting status [%d] to %s in line %d\n", threadId, "SLEEPING", __LINE__) );
	_expControl.threadStatus[threadId] = SLEEPING;
	--_expControl.nActiveThreads;

	changeFreq (threadId, MIN_FREQ);
      }

      /*@ LOG @*/
      fprintf(_expControl.logFile, "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
    }
  
  




/**** FREQ *****/
    void FgcsMalleableFreqArch::adjustFrequency()
    {
      if(_expControl.changeInProgress.tryAcquire()){    //no one is modifying the freq
      }
      else{
	_expControl.futureChange=true;
	return;
      }
      
      updateFreqIdx();
      applyFreq();
      

      _expControl.changeInProgress.release();
      if(_expControl.futureChange) adjustFrequency();
    }

    void FgcsMalleableFreqArch::updateFreqIdx()
    {
      const unsigned int nIdle    = _expControl.currMask.count() - _expControl.nActiveThreads;
      const unsigned int nRunning = _expControl.nActiveThreads;
      static const long double eps2= 0.8; //1.7; //0.001;//Error en las centesimas //2*std::numeric_limits<double>::epsilon();
	    
      const long double maxBudget = _maxBudget;
	    

      _expControl.futureChange=false;

	  
      int fIdx = 0;
      for(int i=_NFREQS-1; i>=0; --i){
	const long double estimation=estimatePower_idx(nIdle, nRunning, i);
	const long double absDiff=std::abs(maxBudget-estimation);
	      
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
    }

    void FgcsMalleableFreqArch::applyFreq()
    {
      const long long maxFreq = getFreqFromIdx_Hz(_expControl.runningFreq_idx);      
//      const long long minFreq = MIN_FREQ;

      
      for(int i=0; i<_expControl.numCores; ++i){
	switch (_expControl.threadStatus[i]){
	case RUNNING:
	  ONLY_DBG( fprintf(stderr, "\t Thread %d -> freq %Ld\n", i, maxFreq ) );
	  changeFreq(i, maxFreq);
	  break;
	// case SLEEPING:
	//   ONLY_DBG( fprintf(stderr, "\t Thread %d -> freq %Ld\n", i, minFreq ) );
	//   changeFreq(i, minFreq);
	//   break;
	// case WAKING_UP:
	//   if(_expControl.nActiveThreads == _expControl.currMask.count()){} //oversubscription
	//   else {
	//     ONLY_DBG( fprintf(stderr, "\t Thread %d -> freq %Ld\n", i, minFreq) );
	//     changeFreq(i, minFreq);
	//   }
	//   break;
	default: break;
	}
	
      }
    }
  }  //namespace ext
}//namespace nanos


DECLARE_PLUGIN("fgcsMalleableFreqArch",nanos::ext::FgcsMalleableFreqArchPlugin);
