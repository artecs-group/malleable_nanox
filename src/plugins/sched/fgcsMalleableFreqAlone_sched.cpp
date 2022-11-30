#include "schedule.hpp"
#include "wddeque.hpp"
#include "plugin.hpp"
#include "system.hpp"
#include "config.hpp"

#include "fgcsMalleableFreqAlone_sched.hpp"
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


      
      FgcsMalleableFreqAlone::FgcsMalleableFreqAlone()  :
	  SchedulePolicy("FGCS malleableFreq"),
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
		  "refTime  = %Lu\n"	   \
		  "currFreq = [];\n"	   \
		  "nActiveThreads = [];\n" \
		  "maxBudget = [];\n"	   \
		  "\n"			   \
		  "\n",
		  refTime);	  

	  
	  //_expControl.startBinding    = getenv(((SMPPlugin*)sys.getSMPPlugin())->getBindingStart(); //ERROR: libreria todavía no cargada
	  _expControl.startBinding    = atoi(txt);
	  _expControl.numCores         = sys.getNumThreads();
	  _expControl.threadStatus     = new ThreadStatus[_expControl.numCores];
	  _expControl.sems             = new sem_t[_expControl.numCores];
	  _expControl.nActiveThreads   = 1;//sys.getNumThreads();

	  /** 
	   * Threads RUNNING or IDLE
	   **/
	  int i=0; 	  //master thread
	  ((SMPPlugin*)sys.getSMPPlugin())->setAffinity(i, i + _expControl.startBinding);
	  _expControl.threadStatus[i] = RUNNING;
	  sem_init(&_expControl.sems[i], 0, 1); //masterThread running
	  
	  //other threads
	  for( i=1; i<_expControl.numCores; i++){
	      ((SMPPlugin*)sys.getSMPPlugin())->setAffinity(i, i + _expControl.startBinding);
	      sem_init(&_expControl.sems[i], 0, 0); 
	      _expControl.threadStatus[i] = IDLE;
	      changeFreq(i, MIN_FREQ);
	  }
	  

	  
	  /** B) ABOUT FREQUENCIES **/
	  //PCAP5
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

      

          char* bg=getenv("MAX_BUDGET");
          if(bg==NULL){
            fprintf(stderr, "ERROR: MAX_BUDGET UNDEFINED\n");
            abort();
          }
          
	  _maxBudget=atof(bg);
          //	  calculateNewFreq();
	  
	  _expControl.appStatus=START;

	  {// NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MOD_ARCH1) );
	    ONLY_DBG( fprintf(stderr, "STARTING SCHEDULING POLICY") );
	  }// NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );

          /*@ @*/
          SS_initLibrary();
          /*@ @*/

	  fprintf(_expControl.logFile,
		  "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
	  fprintf(_expControl.logFile,
		  "maxBudget = [maxBudget; [ %f %f ]];\n", getTimeStamp(), _maxBudget);
	/*@-----@*/
      }

      FgcsMalleableFreqAlone::~FgcsMalleableFreqAlone() 
      {
	  for(int i=0; i<_expControl.numCores; i++)
	      sem_destroy(&_expControl.sems[i]);
	  
	  delete[] _expControl.sems;
	  delete[] _expControl.threadStatus;

	  fclose(_expControl.logFile);
      }
      
      void FgcsMalleableFreqAlone::queue ( BaseThread *thread, WD &wd )
      {
	  //puts("QUEUE 1");
	    
	  TeamData &tdata = (TeamData &) *thread->getTeam()->getScheduleData();
	  tdata._readyQueue->push_back( &wd );
      }
    
      void FgcsMalleableFreqAlone::queue ( BaseThread ** threads, WD ** wds, size_t numElems )
      {
	  //puts("QUEUE 2");
	  
	  ThreadTeam* team = threads[0]->getTeam();
	  TeamData &tdata = (TeamData &) *team->getScheduleData();
	  tdata._readyQueue->push_back( wds, numElems );	
      }
    
    
      void FgcsMalleableFreqAlone::atSuccessor   ( DependableObject &successor, DependableObject &predecessor )
      {
	  //puts("AT SUCCESSOR");	    
	  return;
      }
    
      WD *FgcsMalleableFreqAlone::atSubmit ( BaseThread *thread, WD &newWD )
      {
	  //puts( "AT SUBMIT" );	    
	  queue( thread, newWD );
	  return 0;
      }
    
	    
      WD * FgcsMalleableFreqAlone::atIdle ( BaseThread *thread, int numSteal )
      {
	  if(_expControl.appStatus == STOP) return NULL;
	  
	  //puts("AT IDLE");
	  int threadId = thread->getCpuId();
	    	    
	  //Maybe someone has blocked us
	  sem_wait(&_expControl.sems[threadId]);


	  TeamData &tdata = (TeamData &) *thread->getTeam()->getScheduleData();
	    
	  //Check if we have to wake up more threads.
	  /*** Despertamos threads en este punto para evitar overhead:
	   **    El propio thread encargado de despertar los hilos es el encargado
	   **    de cambiar la frecuencia (y por tanto esperar en los mutex).
	   **    Para evitar que una tarea esté bloqueadea en esos mutex, despertamos a 
	   **    los Workers antes de coger una tarea nueva!
	   ***/
	  if(tdata._readyQueue->size() > 1 && _expControl.nActiveThreads < _expControl.numCores)
            wakeUp_threads(std::min(static_cast<int>(tdata._readyQueue->size()-1), getMaxThreadsWithBudget(_maxBudget) - _expControl.nActiveThreads));

	  
	  WD* wd = tdata._readyQueue->pop_front( thread );
	  
	  int spins=60000;//10;
	  while(wd==NULL && --spins){
	      wd = tdata._readyQueue->pop_front( thread );
	  }
	  
	  if( wd==NULL ){
	    if( threadId==0 ){//|| sys.getTaskNum() > _expControl.nActiveThreads ){
		  sem_post(&_expControl.sems[threadId]);
	      } else {
#ifdef NANOS_INSTRUMENTATION_ENABLED
		  static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
		  static nanos_event_key_t key = ID->getEventKey("fgcs_2");
		  nanos_event_value_t value = 1;
		  sys.getInstrumentation()->raisePointEvents(1, &key, &value);
#endif

		  /*** BLOCK ***/
		  ONLY_DBG( fprintf(stderr, "Sleeping thread %d [nRunning=%d]\n", threadId, _expControl.nActiveThreads) );
		  //If we do not post the sem, the thread will be blocked on the next iteration
		  if(_expControl.appStatus == START) sleepThread(threadId);
	      }
	      
	      return wd;		
	  } 	    
	  return wd;
      }
    
      WD * FgcsMalleableFreqAlone::atPrefetch ( BaseThread *thread, WD &current )
      {
	  //puts("AT PREFETCH");
	  return atIdle(thread,false);
      }
    
      WD * FgcsMalleableFreqAlone::atBeforeExit ( BaseThread *thread, WD &current, bool schedule )
      {
	  //puts("AT BEFORE EXIT");
	  //int coreId=thread->runningOn()->getId() - _expControl.startBinding;
	  int coreId=thread->getCpuId();
	  sem_post(&_expControl.sems[coreId]);
	  return NULL;
      }

      void FgcsMalleableFreqAlone::atShutdown()
      {
	  //puts("AT SHUTDOWN");
	  _expControl.appStatus=STOP;
	    
	  for(int i=0; i<_expControl.numCores; i++){
	      sem_post(&_expControl.sems[i]);
	  }

          
          SS_finishLibrary();
      }

      WD* FgcsMalleableFreqAlone::atWakeUp ( BaseThread *thread, WD &wd )
      {
	  //puts("AT WAKE UP");
	  for(int i=0; i<_expControl.numCores; i++)
	      sem_post(&_expControl.sems[i]);
	  //	    sem_post(&_expControl.sems[0]);
	  return NULL;
      }




      //====================================================================
      void FgcsMalleableFreqAlone::calculateNewFreq()
      {
	  if(_expControl.changeInProgress.tryAcquire()){    //no one is modifying the freq
	  }
	  else{
	      _expControl.futureChange=true;
	      return;
	  }
	  
	  const unsigned int nIdle    = _expControl.numCores - _expControl.nActiveThreads;
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
#ifdef NANOS_INSTRUMENTATION_ENABLED
	  static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
	  static nanos_event_key_t key = ID->getEventKey("fgcs_5");
	  nanos_event_value_t value = getFreqFromIdx_Hz(fIdx);
	  sys.getInstrumentation()->raisePointEvents(1, &key, &value);
#endif

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
      

      void FgcsMalleableFreqAlone::changeFreq(int idx, long long freq) const
      {
	  ((SMPPlugin*)sys.getSMPPlugin())->setFrequency(idx, freq);
	  
	  // const long long currFreq = cpufreq_get_freq_kernel(idx);
	  // if     (currFreq > freq) return decreaseFreq(idx, freq);
	  // else if(currFreq < freq) return increaseFreq(idx, freq);
	  // else                     return;
      }
      
      // void FgcsMalleableFreqAlone::increaseFreq(int idx, long long freq) const
      // {
      //   cpufreq_modify_policy_max(idx, freq);
      //   cpufreq_modify_policy_min(idx, freq);
      // }
      
      // void FgcsMalleableFreqAlone::decreaseFreq(int idx, long long freq) const
      // {
      //   cpufreq_modify_policy_min(idx, freq);
      //   cpufreq_modify_policy_max(idx, freq);
      // }
      
	
	
    //====================================================================
      void FgcsMalleableFreqAlone::sleepThread(int threadId)
      {
	  {  LockBlock l(_expControl.threadStatus_lck);

	    ONLY_DBG( assert( _expControl.threadStatus[threadId] == RUNNING ) );
	      _expControl.threadStatus[threadId] = IDLE;
	      --_expControl.nActiveThreads;
	      changeFreq(threadId, MIN_FREQ);
	      
#ifdef NANOS_INSTRUMENTATION_ENABLED
	static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
	static nanos_event_key_t key = ID->getEventKey("fgcs_4");
	nanos_event_value_t value = _expControl.nActiveThreads;
	sys.getInstrumentation()->raisePointEvents(1, &key, &value);
#endif
	/*@ LOG @*/
	fprintf(_expControl.logFile,
		"nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
	/*@-----@*/
	  }
	  calculateNewFreq();
      }

      void FgcsMalleableFreqAlone::wakeUp_threads(int maxThreads)
      {	  
	  { LockBlock l(_expControl.threadStatus_lck);
	    ONLY_DBG( fprintf(stderr, "Waking up %d thrads [nRunning=%d]\n", maxThreads, _expControl.nActiveThreads) );
	      for(int i=0; i<_expControl.numCores && maxThreads; ++i){
		  if(_expControl.threadStatus[i] == IDLE){
		      _expControl.threadStatus[i] = RUNNING;
		      sem_post(&_expControl.sems[i]);
		      --maxThreads;
		      ++_expControl.nActiveThreads;

#ifdef NANOS_INSTRUMENTATION_ENABLED
		      static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
		      static nanos_event_key_t key = ID->getEventKey("fgcs_3");
		      nanos_event_value_t value = 1;
		      
		      sys.getInstrumentation()->raisePointEvents(1, &key, &value);
#endif
		  }
	      }
	      
#ifdef NANOS_INSTRUMENTATION_ENABLED
	      static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
	      static nanos_event_key_t key= ID->getEventKey("fgcs_4");
	      nanos_event_value_t value = _expControl.nActiveThreads;
	      sys.getInstrumentation()->raisePointEvents(1, &key, &value);
#endif
	      /*@ LOG @*/
	      fprintf(_expControl.logFile,
		      "nActiveThreads = [nActiveThreads; [ %f %d ]];\n", getTimeStamp(), _expControl.nActiveThreads);
	      /*@-----@*/
	      
	      ONLY_DBG( fprintf(stderr, "After waking up: %d threads are running\n", _expControl.nActiveThreads) );
	  }//lockBlock - mutex

	  calculateNewFreq();
      }//wakeUp_threads

  } //namespace ext
}//namespace nanox



DECLARE_PLUGIN("fgcsMalleableFreqAlone",nanos::ext::FgcsMalleableFreqAlonePlugin);
