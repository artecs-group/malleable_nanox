#include "fgcs_dynFreq.hpp"

/***
 ** Distintas políticas para repartir la frecuencia y budget entre
 ** los distintos workers:
 **   A.- [6815558] Los little siempre están a la misma frecuencia y se
 **       intenta maximizar la frecuencia de los big. (asumimos que todos
 **       los big van a estar a la misma frecuencia)
 **
       float powerPerBig = ( MaxBudget -  (nL * pow[lFreq]) / nB;
         for(idx=NFREQS-1; idx>0; --idx)
	   if(_expData.power[idx] <= powerPerBig) break;
 ** 
 **  B.- [2a7438a9] Intentamos maximizar la frecuencia de los cores big a 
 **      costa de minimizar la frecuencia de los little. Esta solución 
 **      siempre intentará perjudicar a los little antes que a los big

      for( bIdx=NFREQS-1; bIdx >= 0; --bIdx ){
        for( lIdx=bIdx; lIdx >= 0; --lIdx ){
	  if(estimatePower_idx(nLittle, lIdx, nBig, bIdx) <= maxBudget){
	     ....; break;
     
 **  C.- [1697b517] Suponiendo que los cores big están a la frecuencia bIdx,
 **      considerar las opciones [bIdx-2 ... bIdx+2] que maximiza el consumo
 **      sin superar el máximo budget.
 **
 **  Modificar las frecuencias durante la ejecución no sirve para nada si no se
 **  tiene en cuenta el grafo (ya sea para tomar decisiones sobre él, o para
 **  modificarlo):
 **
 **  D.- [        ] Si no hay tareas listas para un hilo, éste se duerme a si
 **      mismo y se tiene en cuenta solamente su consumo en idle para realizar
 **      el reparto de budget. OJO: como el consumo en idle también depende de
 **      la frecuencia, antes de dormirse se cambia la frecuencia al mínimo.
***/


/*@ @*/
#include "smpplugin_decl.hpp"
#include <cstdlib>
/*@ @*/


#define SUPPORT_UPDATE 1 

namespace nanos {
namespace ext {

/**-------------------**/
/** ____BotLevCfg____ **/
/**-------------------**/
// Static data for botlev sched

// Initialise default values
int BotLevCfg::updateFreq = 0;
int BotLevCfg::numSpins = 300;
int BotLevCfg::hpFrom = 0;
int BotLevCfg::hpTo = 3;
int BotLevCfg::hpSingle = 0;
int BotLevCfg::steal = 0;
int BotLevCfg::maxBL = 1;
int BotLevCfg::strict = 0;
int BotLevCfg::taskNumber = 0;
NANOS_INSTRUMENT( int BotLevCfg::numCritical=0; )   //! The number of critical tasks (for instrumentation)


/**---------------------------**/
/** ____Experiment Config____ **/
/**---------------------------**/
long double ExperimentsGlobalCfg::MaxBudget=99999;

	
/**----------------------**/
/** ____BotLevDOData____ **/
/**----------------------**/
// Helper class for the computation of bottom levels,
//     One instance of this class is stored inside each dependableObject.
BotLevDOData::BotLevDOData(int tnum, int blev) :
    _taskNumber(tnum), _botLevel(blev), _isReady(false), _lock(), _wd(NULL), _predecessors() { }
BotLevDOData::~BotLevDOData() { }

void BotLevDOData::reset () { _wd = NULL; }
int BotLevDOData::getTaskNumber() const { return _taskNumber; }
int BotLevDOData::getBotLevel() const { return _botLevel; }
bool BotLevDOData::setBotLevel( int blev ) {
  if ( blev > _botLevel ) {
    {
      LockBlock lock1( _lock );
      _botLevel = blev;
    }
    return true;
  }
  return false;
}

bool BotLevDOData::getReady() const { return _isReady; }
void BotLevDOData::setReady() { _isReady = true; }
 
void BotLevDOData::setCriticality( short c ) { _isCritical = c; }
short BotLevDOData::getCriticality() { return _isCritical; }

WD* BotLevDOData::getWorkDescriptor() const { return _wd; }
void BotLevDOData::setWorkDescriptor(WD *wd) { _wd = wd; }

std::set<BotLevDOData *>* BotLevDOData::getPredecessors() { return &_predecessors; }
void BotLevDOData::addPredecessor(BotLevDOData *predDod) { _predecessors.insert(predDod); }




/**----------------**/
/** ____BotLev____ **/
/**---------------**/
BotLev::BotLev() : SchedulePolicy ( "fgcsDynFreq" ), _expData(sys.getNumThreads()) {
  sys.setPredecessorLists(true);
  _currMax = _maxBotLev = BotLevCfg::maxBL;


  char* mb;
  if((mb=getenv("MAX_BUDGET"))!=NULL)
    ExperimentsGlobalCfg::MaxBudget=atof(mb);
  
  ONLY_DBG( fprintf(stderr, "MaxBudget: %Lf\n", ExperimentsGlobalCfg::MaxBudget) );
  const unsigned nWorkers = sys.getNumThreads();

  {
    LockBlock l(_expData.mtxArch);
    calculateNewFreqs();
  }

  const long long runningFreq = getFreqFromIdx_Hz(_expData.runningFreq);
  
  for(unsigned i=0; i<nWorkers; i++)
    ((SMPPlugin*)(sys.getSMPPlugin()))->setFrequency(i, runningFreq);
  
}
    
BotLev::~BotLev() {
}
       
size_t BotLev::getTeamDataSize () const { return sizeof(TeamData); }       
size_t BotLev::getThreadDataSize () const { return 0; }

ScheduleTeamData* BotLev::createTeamData () { return NEW TeamData(); }
ScheduleThreadData* BotLev::createThreadData () { return 0; }


/*!
 *  \brief Enqueue a work descriptor in the readyQueue of the passed thread
 *  \param thread pointer to the thread to which readyQueue the task must be appended
 *  \param wd a reference to the work descriptor to be enqueued
 *  \sa ThreadData, WD and BaseThread
 */
void BotLev::queue ( BaseThread *thread, WD &wd ) {
#ifdef NANOS_INSTRUMENTATION_ENABLED
  int criticality; 
#endif
  TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
  // Find the priority
  DependableObject *dos = wd.getDOSubmit();
  BotLevDOData *dodata;
  unsigned int priority = 0; 
  short qId = -1;
  if ( dos ){ //&& numThreads>1 ) {
    dodata = (BotLevDOData *)dos->getSchedulerData();
    dodata->setWorkDescriptor(&wd);
    priority = dodata->getBotLevel();
  }
  else {
    if(wd.getDepth() == 0 ) { 
      wd.setPriority(0);
      data._readyQueues[2].push_back( &wd ); //readyQueue number 2 is for the main and implicit tasks
    }
    else {  //in this case numThreads = 1 inserting in queue number 2
      data._readyQueues[2].push_back( &wd );
    }
#ifdef NANOS_INSTRUMENTATION_ENABLED
    criticality = 3; 
    if(wd.getSchedulerData() == NULL) { 
      WDData * wddata = new WDData();
      wddata->setCriticality(criticality); 
      wd.setSchedulerData((ScheduleWDData*)wddata, true);
    } 
    else { 
      WDData & scData = *dynamic_cast<WDData*>( wd.getSchedulerData() ); 
      scData.setCriticality(criticality); 
    }
#endif
    return;
  }
  wd.setPriority(priority);
  /* Critical tasks' consideration
     1st case: Detection of a new longest path (wdPriority > maxPriority)
     2nd case: Detection of the next critical task in the current lontest path (belongs in _topSuccessors)
     3rd case: The remaining tasks are not critical
  */
  if( ( wd.getPriority() >  _maxBotLev ) || 
      ( !BotLevCfg::strict && (wd.getPriority() ==  _maxBotLev)) ) {
    //The task is critical
    {
      LockBlock l(_botLevLock);
      _maxBotLev = wd.getPriority();
      _currMax = _maxBotLev;
      _topSuccesors = (dos->getSuccessors());
      NANOS_INSTRUMENT( BotLevCfg::numCritical++; )
    }
    dodata->setCriticality(1);
    qId = 1;
    NANOS_INSTRUMENT ( criticality = 1; )
  }
  else if( ((_topSuccesors.find( std::make_pair( wd.getId(), dos ) )) != (_topSuccesors.end()))
           && wd.getPriority() >= _currMax-1 ) {
    //The task is critical
    {
      LockBlock l(_botLevLock);
      _currMax = wd.getPriority();
      _topSuccesors = (dos->getSuccessors());
      NANOS_INSTRUMENT( BotLevCfg::numCritical++; )
    }
    dodata->setCriticality(1);
    qId = 1;
    NANOS_INSTRUMENT ( criticality = 1; )
  }
  else
  {
    //Non-critical task
    dodata->setCriticality(2);
    qId = 0;
    NANOS_INSTRUMENT ( criticality = 2; )
  }
  data._readyQueues[qId].push_back( &wd ); //queues 0 or 1
  dodata->setReady();

#ifdef NANOS_INSTRUMENTATION_ENABLED
  if(wd.getSchedulerData() == NULL) {
    WDData * wddata = new WDData();
    wddata->setCriticality(criticality);
    wd.setSchedulerData((ScheduleWDData*)wddata, true);

  }
  else {
    WDData & scData = *dynamic_cast<WDData*>( wd.getSchedulerData() );
    scData.setCriticality(criticality);
  }
  WDData & wddata = *dynamic_cast<WDData*>( wd.getSchedulerData() );
  if(wddata.getCriticality() == 1) {
    NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
    NANOS_INSTRUMENT ( static nanos_event_key_t critical_wd_id = ID->getEventKey("critical-wd-id"); )
    NANOS_INSTRUMENT ( nanos_event_key_t crit_key[1]; )
    NANOS_INSTRUMENT ( nanos_event_value_t crit_value[1]; )
    NANOS_INSTRUMENT ( crit_key[0] = critical_wd_id; )
    NANOS_INSTRUMENT ( crit_value[0] = (nanos_event_value_t) wd.getId(); )
    NANOS_INSTRUMENT( sys.getInstrumentation()->raisePointEvents(1, crit_key, crit_value); )
  }
      
#endif
  return;
}

void BotLev::updateBottomLevels( BotLevDOData *dodata, int botLev ) {
  std::vector<BotLevDOData *> stack;
	   
  dodata->setBotLevel(botLev); // Set the bottom level, and add to the stack
  stack.push_back(dodata);
	   
  if( dodata->getReady() )  return; //Task is ready so, there are no predecessors to be updated

  while ( !stack.empty() ) { // A depth first traversal
    // Pop an element from the stack and get its level
    BotLevDOData *bd = stack.back();
    stack.pop_back();

    int botLevNext = bd->getBotLevel() + 1;
    bool changed = false;
    // Deal with the predecessors
    std::set<BotLevDOData *> *preds = bd->getPredecessors();
                  
    for ( std::set<BotLevDOData *>::iterator pIter = preds->begin(); pIter != preds->end(); pIter++ ) {
      BotLevDOData *pred = *pIter;
      if(botLevNext > pred->getBotLevel()) {
        changed = pred->setBotLevel( botLevNext );
        stack.push_back( pred );
      }
#if SUPPORT_UPDATE
      //Reorder the readyQueues if the work descriptor with updated priority is ready
      WD * predWD = pred->getWorkDescriptor();
      if ( changed && pred->getReady() && predWD ) {
        TeamData &data = ( TeamData & ) *myThread->getTeam()->getScheduleData();
        short criticality = pred->getCriticality();
        if(criticality == 1)
          data._readyQueues[1].reorderWD(predWD);
        else if(criticality == 2)
          data._readyQueues[0].reorderWD(predWD);
      } 
#endif
    }
  }
}


void BotLev::atCreate ( DependableObject &depObj )
{
  //! Creating DO scheduler data
  BotLevDOData *dodata = new BotLevDOData(++BotLevCfg::taskNumber, 0);
  depObj.setSchedulerData( (DOSchedulerData*) dodata );

  DepObjVector predecessors;
  { 
    LockBlock l(depObj.getLock());
    predecessors = depObj.getPredecessors();
  }
  for ( DepObjVector::iterator it = predecessors.begin(); it != predecessors.end(); it++ ) {
    DependableObject *pred = it->second;
    if (pred) {
      BotLevDOData *predObj = (BotLevDOData *)pred->getSchedulerData();
      if (predObj) {
        dodata->addPredecessor(predObj);
      }
    }
  }

  //! When reaching the threshold we need to update bottom levels,
  //! otherwise we push dodata in a temporal stack
  if ( _blStack.size() + 1 >= (unsigned int) BotLevCfg::updateFreq ) {
    updateBottomLevels(dodata, 0);
    while ( !_blStack.empty() ) {
      BotLevDOData *dodata1 = _blStack.top();
      if (dodata1->getBotLevel() <= 0) 
        updateBottomLevels(dodata1, 0);
      {
        LockBlock l(_stackLock);
        _blStack.pop();
      }
    }
  } else {
    LockBlock l(_stackLock);
    _blStack.push( dodata );
  }
}

/*!
 * Function called when a new task must be created: the new created task
 *          is directly queued (Breadth-First policy)
 * \param thread pointer to the thread to which belongs the new task
 * \param wd a reference to the work descriptor of the new task
 */
WD* BotLev::atSubmit ( BaseThread *thread, WD &newWD )
{
  queue(thread, newWD);
  return 0;
}

/*! 
 *  \brief Function called by the scheduler when a thread becomes idle to schedule it: implements the CILK-scheduler algorithm
 *  \param thread pointer to the thread to be scheduled
 *  \sa BaseThread
 */
WD * BotLev::atIdle ( BaseThread *thread, int numSteal ) {    
  if(_expData.appStatus==APP_STOPPING) return NULL;

  WorkDescriptor * wd = NULL;
  TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
  
#ifdef NANOS_INSTRUMENTATION_ENABLED
  static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
  static nanos_event_key_t critical_id = ID->getEventKey("n-crit");
  static nanos_event_key_t non_critical_id = ID->getEventKey("n-noCrit");
	    
  nanos_event_key_t keys[2];
  nanos_event_value_t values[2];
  keys[0] = critical_id;      values[0] = data._readyQueues[1].size();
  keys[1] = non_critical_id;  values[1] = data._readyQueues[0].size();

  sys.getInstrumentation()->raisePointEvents(2, keys, values);
#endif
    


  const int coreId = thread->runningOn()->getId();
  int spins=BotLevCfg::numSpins;

  //check if we're blocked
  sem_wait(&_expData.sems[coreId]);
  
  while(wd==NULL && --spins){
                   wd=data._readyQueues[1].pop_front(thread);          //critical
      if(wd==NULL) wd=data._readyQueues[0].pop_front(thread);      //non-critical
      if(wd==NULL) wd=data._readyQueues[2].pop_front(thread);            //others
  }
  

  
  if(wd==NULL){ //SLEEP
      if(coreId!=0 && _expData.appStatus==APP_RUNNING){
ONLY_DBG(fprintf(stderr, "SLEEP\n") );
	  LockBlock l(_expData.mtxArch);
	  //Do not release mutex. We're going to sleep
	  // sem_post(&_expControl.sems[threadId]); 
	  
	  _expData.coreStatus[coreId] = CORE_IDLE;
	  --_expData.nRunning;
	  
	  calculateNewFreqs();
      }
  } else if( (_expData.nRunning != _expData.nWorkers) &&
	     (data._readyQueues[1].size() > 0 || data._readyQueues[0].size() > 0)){
ONLY_DBG(fprintf(stderr, "WAKE UP\n"));      
      wakeUpWorker();
  }
  
  return wd;
}
       
WD * BotLev::atBeforeExit(BaseThread *thread, WD &wd, bool schedule) {
    const int coreId=thread->runningOn()->getId();
    sem_post(&_expData.sems[coreId]);
    return NULL;
}

void BotLev::atShutdown() {
    _expData.appStatus = APP_STOPPING;
    
    for(int coreIdx=0; coreIdx<_expData.nWorkers; ++coreIdx)
	sem_post(&_expData.sems[coreIdx]);

}

//=======================================================================
//
// EXPS
//
//=======================================================================
void BotLev::calculateNewFreqs()
{
  ONLY_DBG( assert(_expData.mtxArch.getState() != NANOS_LOCK_FREE) );

  const unsigned nIdle    = _expData.nWorkers - _expData.nRunning;
  const unsigned nRunning = _expData.nRunning;
  const long double maxBudget  = ExperimentsGlobalCfg::MaxBudget;

  const long double eps2=2*std::numeric_limits<long double>::epsilon();
  
  for(int fIdx=NFREQS-1; fIdx >=0; --fIdx){
      const long double estimation=estimatePower_idx(nIdle, nRunning, fIdx);
      const long double absDiff=std::abs(maxBudget-estimation);

      if((estimation<=maxBudget) || absDiff<=eps2){
	  _expData.runningFreq=fIdx;
	  break;
      }
  }

 
  modifyFreqs();

#ifdef NANOS_DEBUG_ENABLED
  fprintf(stderr, "AAA: Setting freq: %d\n", _expData.runningFreq);
#endif
	
}

void BotLev::modifyFreqs(){
    for(int i=0; i<_expData.nWorkers; i++){

	const long long goalFreq =
	    (_expData.coreStatus[i]==CORE_RUNNING)
	    ? getFreqFromIdx_Hz(_expData.runningFreq)
	    : getFreqFromIdx_Hz(_expData.idleFreq);
	const long long currFreq = cpufreq_get_freq_kernel(i);

	
	//fprintf(stderr, "Changing big freq to: %Ld\n", goalFreq);
	
	if(currFreq == goalFreq) return;
	else if(currFreq < goalFreq){
	    cpufreq_modify_policy_max(i, goalFreq);
	    cpufreq_modify_policy_min(i, goalFreq);
	    //      cpufreq_set_frequency(i, goalFreq);
	} else {
	    cpufreq_modify_policy_min(i, goalFreq);
	    cpufreq_modify_policy_max(i, goalFreq);
	    //      cpufreq_set_frequency(i, goalFreq);
	}
    }
}

void BotLev::wakeUpWorker()
{
    LockBlock l(_expData.mtxArch);
    for(int coreIdx=0; coreIdx<_expData.nWorkers; ++coreIdx){
	if(_expData.coreStatus[coreIdx]==CORE_IDLE){

	    _expData.coreStatus[coreIdx]=CORE_RUNNING;
	    ++_expData.nRunning;
	    calculateNewFreqs();

	    sem_post(&_expData.sems[coreIdx]);
	    return;
	}
    }
}



//=======================================================================
//=======================================================================     
/**-------------------------**/
/** ____fgcsSchedPlugin____ **/
/**-------------------------**/
fgcsDynFreqPlugin::fgcsDynFreqPlugin() : Plugin( "Dynamic frequency adaptation based on BotLev",1 ) {}       
void fgcsDynFreqPlugin::config ( Config &config_ ) {
  config_.setOptionsSection( "Bottom level", "Bottom-level scheduling module" );
  config_.registerConfigOption ( "update-freq", new Config::PositiveVar( BotLevCfg::updateFreq ), "Defines how often to update the bottom levels" );
  config_.registerArgOption ( "update-freq", "update-freq" );
  config_.registerEnvOption ( "update-freq", "NX_BL_FREQ" );
      
  config_.registerConfigOption ( "numSpins", new Config::PositiveVar( BotLevCfg::numSpins ), "Defines the number of spins in atIdle (work stealing)" );
  config_.registerArgOption ( "numSpins", "numSpins" );
  config_.registerEnvOption ( "numSpins", "NX_NUM_SPINS" );
      
  config_.registerConfigOption ( "from", new Config::IntegerVar( BotLevCfg::hpFrom ), "Sets the thread id of the first fast core" );
  config_.registerArgOption ( "from", "from" );
  config_.registerEnvOption ( "from", "NX_HP_FROM" );
      
  config_.registerConfigOption ( "to", new Config::PositiveVar( BotLevCfg::hpTo ), "Sets the thread id of the last fast core" );
  config_.registerArgOption ( "to", "hpTo" );
  config_.registerEnvOption ( "to", "NX_HP_TO" );
      
  config_.registerConfigOption ( "single", new Config::IntegerVar( BotLevCfg::hpSingle ), "Sets the thread id of a single fast core" );
  config_.registerArgOption ( "single", "hpSingle" );

  config_.registerConfigOption ( "maxBlev", new Config::IntegerVar( BotLevCfg::maxBL ), "Defines the initial value of maximum bottom level" );
  config_.registerArgOption ( "maxBlev", "maxBlev" );
  config_.registerEnvOption ( "maxBlev", "NX_MAXB" );
      
  config_.registerConfigOption ( "strict", new Config::IntegerVar( BotLevCfg::strict ), "Defines whether we use strict policy. (Strict -- les crit. tasks, Flexible -- more crit. tasks)" );
  config_.registerArgOption ( "strict", "strict" );
  config_.registerEnvOption ( "strict", "NX_STRICTB" );
      
  config_.registerConfigOption ( "steal", new Config::IntegerVar( BotLevCfg::steal ), "Defines if we use bi-directional work stealing fast <--> slow" );
  config_.registerArgOption ( "steal", "steal" );
  config_.registerEnvOption ( "steal", "NX_STEALB" );
     
}

void fgcsDynFreqPlugin::init() {
  sys.setDefaultSchedulePolicy(new BotLev());
}
    
} //ext namespace
} //nanos namespace


DECLARE_PLUGIN("fgcsDynFreq",nanos::ext::fgcsDynFreqPlugin);
