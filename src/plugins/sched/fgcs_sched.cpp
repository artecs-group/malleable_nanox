#include "fgcs_sched.hpp"

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
	int BotLevCfg::experiment_type = 0;



	/**---------------------------**/
	/** ____Experiment Config____ **/
	/**---------------------------**/
	float ExperimentsCfg::MaxBudget   = 11.0;
	float ExperimentsCfg::LittlePower = 1.9754;
	float ExperimentsCfg::BigPower    = 5.043375;
	int   ExperimentsCfg::nCoresUsed  = 4;
	    
	int   ExperimentsCfg::bigFreq     = 2000000;
	int   ExperimentsCfg::littleFreq  = 1000000;

	
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
	BotLev::BotLev() : SchedulePolicy ( "BotLev" ) {
	    sys.setPredecessorLists(true);
	    _currMax = _maxBotLev = BotLevCfg::maxBL;

	    this->exp14_ctr();
	    //this->exp15_ctr();
	}
       
	BotLev::~BotLev() {
	    this->exp14_dtr();
	    //this->exp15_dtr();
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
	void BotLev::queue ( BaseThread *thread, WD &wd )
	{
	    exp14_queue(thread, wd);
	    //exp15_queue(thread, wd);
	}

	void BotLev::updateBottomLevels( BotLevDOData *dodata, int botLev )
	{
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
	WD * BotLev::atIdle ( BaseThread *thread, int numSteal )
	{
	    return this->exp14_atIdle(thread, numSteal);
	    //return this->exp15_atIdle(thread, numSteal);
	}
       
	WD * BotLev::atBeforeExit(BaseThread *thread, WD &wd, bool schedule) {
	    return this->exp14_atBeforeExit(thread, wd, schedule);
	    //return this->exp15_atBeforeExit(thread, wd, schedule);
	}


	void BotLev::atShutdown() {
	    exp14_atShutdown();
	    //exp15_atShutdown();
	}
	//=======================================================================
	//
	// EXPS
	//
	//=======================================================================
       
	//-----------------------------------------------------------------------
	// exp00                                                               //
	//-----------------------------------------------------------------------
	WD* BotLev::exp00_atIdle ( BaseThread *thread, int numSteal )
	{
	    WorkDescriptor * wd;
	    TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
	    unsigned int spins = BotLevCfg::numSpins;

	    //Separation of big and small cores - big cores execute from queue 1 - small cores execute from queue 2
	    if( ((thread->runningOn()->getId() >= BotLevCfg::hpFrom && thread->runningOn()->getId() <= BotLevCfg::hpTo) || 
		 ( BotLevCfg::hpSingle && thread->runningOn()->getId() == BotLevCfg::hpSingle )) ) {
		//Big core
		wd = data._readyQueues[1].pop_front( thread );
		while( wd == NULL && spins )
		{
		    wd = data._readyQueues[1].pop_front( thread );
		    spins--;
		}
		if(!wd ) {
		    // default work stealing: big stealing from small
		    wd = data._readyQueues[0].pop_front( thread );
		}
	    }
	    else {
		//Small core
		wd = data._readyQueues[0].pop_front( thread );
		while( wd == NULL && spins )
		{
		    wd = data._readyQueues[0].pop_front( thread );
		    spins--;
		}
		if(!wd && BotLevCfg::steal) {
		    //optionally: small stealing from big
		    wd = data._readyQueues[1].pop_front( thread );
		}
	       
	    }
	   
	    if(!wd) wd = data._readyQueues[2].pop_front( thread ); 
	    return wd;
	}

	void BotLev::exp00_queue ( BaseThread *thread, WD &wd )
	{
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

		return;
	    }
	   
	    wd.setPriority(priority);
	    /* Critical tasks' consideration
	       1st case: Detection of a new longest path (wdPriority > maxPriority)
	       2nd case: Detection of the next critical task in the current longest path (belongs in _topSuccessors)
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
		}
		dodata->setCriticality(1);
		qId = 1;
	    }
	    else if( ((_topSuccesors.find( std::make_pair( wd.getId(), dos ) )) != (_topSuccesors.end()))
		     && wd.getPriority() >= _currMax-1 ) {
		//The task is critical
		{
		    LockBlock l(_botLevLock);
		    _currMax = wd.getPriority();
		    _topSuccesors = (dos->getSuccessors());
		}
		dodata->setCriticality(1);
		qId = 1;
	    }
	    else
	    {
		//Non-critical task
		dodata->setCriticality(2);
		qId = 0;
	    }
	    data._readyQueues[qId].push_back( &wd ); //queues 0 or 1
	    dodata->setReady();


	    /*@ +++++++++ @*/
	    data._lockfout->acquire();
	    fprintf(data._fout, "data0 = [ data0; %u %lu ];\n", (unsigned)time(NULL), data._readyQueues[0].size());
	    fprintf(data._fout, "data1 = [ data1; %u %lu ];\n", (unsigned)time(NULL), data._readyQueues[1].size());
	    data._lockfout->release();
	    /*@ +++++++++ @*/

	    return;
	}
       

       
	//-----------------------------------------------------------------------
	// exp11                                                               //
	//-----------------------------------------------------------------------
// 	WD* BotLev::exp11_atIdle ( BaseThread *thread, int numSteal ) 
// 	{
// 	    WorkDescriptor * wd;
// 	    TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
// //         unsigned int spins = BotLevCfg::numSpins;

// 	    int coreId = thread->runningOn()->getId();	 
// 	    long long currFreq = cpufreq_get_freq_kernel(coreId);	 
// 	    if( (wd=data._readyQueues[1].pop_front(thread)) != NULL){
// 		/* If there are critical tasks, convert the core into BIG, and run the task
// 		 */
// 		if(currFreq != BIG_FREQ) {
// 		    cpufreq_modify_policy_max(coreId, BIG_FREQ);
// 		    cpufreq_modify_policy_min(coreId, BIG_FREQ);
// 		    cpufreq_set_frequency(coreId, BIG_FREQ);
// 		}
// 	    } else if( (wd=data._readyQueues[0].pop_front(thread)) != NULL){
// 		/* If there are not critical tasks, convert the core into LITTLE, and run the task
// 		 */
// 		long long destFreq =  LITTLE_FREQ_MAX;;
// 		if(currFreq != destFreq) {
// 		    cpufreq_modify_policy_min(coreId, destFreq);
// 		    cpufreq_modify_policy_max(coreId, destFreq);
// 		    cpufreq_set_frequency(coreId, destFreq);
// 		}
// 	    } else {
// 		wd = data._readyQueues[2].pop_front( thread ); 
// 	    }
	 
// 	    return wd;
// 	}

// 	void BotLev::exp11_queue ( BaseThread *thread, WD &wd )
// 	{
// 	    exp00_queue(thread,wd);
// 	}
       
	//-----------------------------------------------------------------------
	// exp12                                                               //
	//-----------------------------------------------------------------------
// 	WD* BotLev::exp12_atIdle ( BaseThread *thread, int numSteal )
// 	{
// 	    WorkDescriptor * wd;
// 	    TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
// //         unsigned int spins = BotLevCfg::numSpins;

// 	    int coreId = thread->runningOn()->getId();	 
// 	    long long currFreq = cpufreq_get_freq_kernel(coreId);	 
// 	    if( (wd=data._readyQueues[1].pop_front(thread)) != NULL){
// 		/* If there are critical tasks, convert the core into BIG, and run the task
// 		 */
// 		if(currFreq != BIG_FREQ) {
// 		    cpufreq_modify_policy_max(coreId, BIG_FREQ);
// 		    cpufreq_modify_policy_min(coreId, BIG_FREQ);
// 		    cpufreq_set_frequency(coreId, BIG_FREQ);
// 		}
// 	    } else if( (wd=data._readyQueues[0].pop_front(thread)) != NULL){
// 		long long destFreq = LITTLE_FREQ_MIN;
// 		if(currFreq != destFreq) {
// 		    cpufreq_modify_policy_min(coreId, destFreq);
// 		    cpufreq_modify_policy_max(coreId, destFreq);
// 		    cpufreq_set_frequency(coreId, destFreq);
// 		}
// 	    } else {
// 		wd = data._readyQueues[2].pop_front( thread ); 
// 	    }
// 	    return wd;
// 	}

// 	void BotLev::exp12_queue ( BaseThread *thread, WD &wd )
// 	{
// 	    exp00_queue(thread,wd);
// 	}
       
//        //-----------------------------------------------------------------------
//        // exp13                                                               //
//        //-----------------------------------------------------------------------       
//        void BotLev::exp13_ctr()
//        {
// 	   for(int i=0; i<NCORES_PER_APP; i++){
// 	       sem_init(&_expControl._sems[i], 0, 1);
// 	       _expControl._coreType[i] = CORELITTLE;	       
// 	   }

	   
// 	   _expControl._archIdx = 0;
// 	   _expControl._archMaxIdx = 2;
// 	   //                          BIG                           LITTLE
// 	   _expControl._archConf[0][0] = 0; _expControl._archConf[0][1] = 4;
// 	   _expControl._archConf[1][0] = 1; _expControl._archConf[1][1] = 3;
// 	   _expControl._archConf[2][0] = 2; _expControl._archConf[2][1] = 1;
// //	   _expControl._archConf[3][0] = 4; _expControl._archConf[3][1] = 1;
//        }
       
//        void BotLev::exp13_dtr()
//        {
// 	   for(int i=0; i<NCORES_PER_APP; i++)
// 	       sem_destroy(&_expControl._sems[i]);
//        }
       
//        WD*  BotLev::exp13_atIdle ( BaseThread *thread, int numSteal )
//        {
// 	   WorkDescriptor * wd;
// 	   TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
// 	   //unsigned int spins = BotLevCfg::numSpins;
// 	   unsigned int spins = 0;
// 	   int coreId = thread->runningOn()->getId();

//            //A) get our own semaphore
// 	   sem_wait(&_expControl._sems[coreId]);
	   

// 	   //B) Execute
// 	   //Separation of big and small cores - big cores execute from queue 1 - small cores execute from queue 2
// 	   if( _expControl._coreType[thread->runningOn()->getId()] == COREBIG ) {
// 	       //Big core
// 	       wd = data._readyQueues[1].pop_front( thread );
// 	       while( wd == NULL && spins )
// 	       {
// 		   wd = data._readyQueues[1].pop_front( thread );
// 		   spins--;
// 	       }
// 	       if(!wd ) {
// 		   // default work stealing: big stealing from small
// 		   wd = data._readyQueues[0].pop_front( thread );
// 	       }
// 	   }
// 	   else {
// 	       //Small core
// 	       wd = data._readyQueues[0].pop_front( thread );
// 	       while( wd == NULL && spins )
// 	       {
// 		   wd = data._readyQueues[0].pop_front( thread );
// 		   spins--;
// 	       }
// 	       if(!wd && BotLevCfg::steal) {
// 		   //optionally: small stealing from big
// 		   wd = data._readyQueues[1].pop_front( thread );
// 	       }
	       
// 	   }
	   
// 	   if(!wd) wd = data._readyQueues[2].pop_front( thread ); 
// 	   if(wd == NULL)
// 	       sem_post(&_expControl._sems[coreId]);

// 	   return wd;
//        }
       

//        WD*  BotLev::exp13_atBeforeExit( BaseThread *thread, WD &wd, bool schedule)
//        {
// 	   // Free the resource
// 	   sem_post(&_expControl._sems[thread->runningOn()->getId()]);

// 	   // Check if we should modify the architecture
// 	   TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
// 	   int nCrit    = data._readyQueues[1].size();
// 	   int nNotCrit = data._readyQueues[0].size();
	   
// 	   exp13_modArch( nCrit, nNotCrit);

	   
// 	   return NULL; 
//        }

//        void BotLev::exp13_queue ( BaseThread *thread, WD &wd )
//        {
// 	   //A) queue the work
// 	   exp00_queue(thread,wd);

// 	   //B) modify the architecture if needed
// 	   TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
// 	   int nCrit    = data._readyQueues[1].size();
// 	   int nNotCrit = data._readyQueues[0].size();
	   
// 	   exp13_modArch( nCrit, nNotCrit);
//        }

       
//        void BotLev::exp13_modArch(int nCrit, int nNotCrit)
//        {
// #ifdef NANOS_INSTRUMENTATION_ENABLED
// 	 static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
// 	 static nanos_event_key_t critical_id = ID->getEventKey("n-crit");
// 	 static nanos_event_key_t non_critical_id = ID->getEventKey("n-noCrit");
// 	 nanos_event_key_t keys[2];
// 	 nanos_event_value_t values[2];
// 	 keys[0] = critical_id; values[0] = nCrit;
// 	 keys[1] = non_critical_id; values[1] = nNotCrit;
// 	 sys.getInstrumentation()->raisePointEvents(2, keys, values);
// #endif
// 	   //fprintf(stderr, "Status [%d, %d]\n", nCrit, nNotCrit);
// 	   //fprintf(stderr, "archIdx [%d]\n", _expControl._archIdx);

	   
// 	   if(nCrit > _expControl._archConf[_expControl._archIdx][0] && //more critic tasks than big Cores
// 	      _expControl._archIdx < _expControl._archMaxIdx         && //we can have more BIG cores
// 	      _expControl._mtxArch.tryAcquire() )                         //we're in charge of the change
// 	   {
// 	       if(nCrit <= _expControl._archConf[_expControl._archIdx][0]){
// 		   _expControl._mtxArch.release();
// 		   return;
// 	       }

// 	       NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MODARCH) );

// 	       /**                               **/
// 	       /* _INCREMENT_ NUMBER OF BIG CORES */
// 	       /**                               **/
	       
// 	       int currentBig    = _expControl._archConf[_expControl._archIdx][0];
// 	       int currentLittle = _expControl._archConf[_expControl._archIdx][1];
	       
// 	       int destBig    = _expControl._archConf[_expControl._archIdx + 1 ][0];
// 	       int destLittle = _expControl._archConf[_expControl._archIdx + 1 ][1];

// 	       int toIncrease = destBig - currentBig;
// 	       // int toDecrease = currentLittle - destLittle;
// 	       // int toModify   = toIncrease + toDecrease;

// 	       //A.- Lock all the cores affected:
// 	       {
// 		   //A1.- Little to Big
// 		   for(int i=currentBig; i<currentBig + toIncrease; i++){
// 		       sem_wait(&_expControl._sems[i]);
// 		       //fprintf(stderr, "A1- little to big [%d]\n", i);		   
// 		   }
		   
// 		   //A2.- Little to Idle
// 		   for(int i=destBig+destLittle; i<currentBig + currentLittle; i++){
// 		       //fprintf(stderr, "A2- little to idle [%d] - pre", i);		   
// 		       sem_wait(&_expControl._sems[i]);
// 		       //fprintf(stderr, " post\n");		   

// 		   }
// 	       }
	       
// 	       //B.- Change frequency Little to big
// 	       {
		   
// 		   long long currFreq;
// 		   for(int i=currentBig; i<currentBig + toIncrease; i++){
// 		       currFreq = cpufreq_get_freq_kernel(i);

// 		       //fprintf(stderr, "B.- Changing freq to BIG (from %Ld) [%d]\n", currFreq, i);
// 		       if(currFreq != BIG_FREQ) {
// 			   cpufreq_modify_policy_max(i, BIG_FREQ);
// 			   cpufreq_modify_policy_min(i, BIG_FREQ);
// 			   cpufreq_set_frequency(i, BIG_FREQ);
// 		       }
// 		       _expControl._coreType[i] = COREBIG;
// 		   }
		   
// 	       }
	       	       
//  	       //C.- Release the new Big cores
// 	       {
// 		   for(int i=currentBig; i<currentBig + toIncrease; i++){
// 		       sem_post(&_expControl._sems[i]);
// 		       //fprintf(stderr, "C- Release BIG [%d]\n", i);		   
// 		   }
// 		   //fprintf(stderr, "\n");
// 	       }

// 	       //D.- Update and release the Arch Struct
// 	       ++ _expControl._archIdx;
// 	       _expControl._mtxArch.release();

//    	       NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
// 	   }
// 	   else if(
// 	       nCrit < _expControl._archConf[_expControl._archIdx][0]      && //less critic tasks than big Cores
// 	       _expControl._archIdx > 0                                    && //we can have more LITTLE cores
// 	       nCrit <= _expControl._archConf[_expControl._archIdx -1 ][0] && //less critic tasks than big Cores
// 	       _expControl._mtxArch.tryAcquire() )                            //we're in charge of the change
// 	   {
// 	       if(nCrit >= _expControl._archConf[_expControl._archIdx][0]){
// 		   _expControl._mtxArch.release();
// 		   return ;
// 	       }

// 	       NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MODARCH) );	       
// 	       /**                               **/
// 	       /* _DECREMENT_ NUMBER OF BIG CORES */
// 	       /**                               **/
// 	       int currentBig    = _expControl._archConf[_expControl._archIdx][0];
// 	       int currentLittle = _expControl._archConf[_expControl._archIdx][1];
	       
// 	       int destBig    = _expControl._archConf[_expControl._archIdx - 1 ][0];
// 	       int destLittle = _expControl._archConf[_expControl._archIdx - 1 ][1];

// 	       //int toIncrease = destBig - currentBig;
// 	       int toDecrease = currentBig - destBig;
// 	       // int toModify   = toIncrease + toDecrease;

// 	       //A.- Lock all the cores affected:
// 	       {
// 		   //A1.- Big to Little
// 		   for(int i=currentBig-1; i>currentBig - toDecrease - 1; i--){
// 		       sem_wait(&_expControl._sems[i]);
// 		       //fprintf(stderr, "A1- big to Little [%d]\n", i);		   
// 		   }
// 	       }
       
// 	       //B.- Change frequency big to little
// 	       {
		   
// 		   long long currFreq;
// 		   for(int i=currentBig-1; i>currentBig - toDecrease -1 ; i--){
// 		       currFreq = cpufreq_get_freq_kernel(i);

// 		       //fprintf(stderr, "B.- Changing freq to LITTLE (from %Ld) [%d]\n", currFreq, i);
// 		       if(currFreq != LITTLE_FREQ_MIN) {
// 			   cpufreq_modify_policy_max(i, LITTLE_FREQ_MIN);
// 			   cpufreq_modify_policy_min(i, LITTLE_FREQ_MIN);
// 			   cpufreq_set_frequency(i, LITTLE_FREQ_MIN);
// 		       }
// 		       _expControl._coreType[i] = CORELITTLE;
// 		   }
		   
// 	       }
	       
//  	       //C.- Release the new cores
// 	       {
// 		   //C1.- Big to little
// 		   for(int i=currentBig-1; i>currentBig - toDecrease -1 ; i--){
// 		       sem_post(&_expControl._sems[i]);
// 		       //fprintf(stderr, "C1- Release LITTLE [%d]\n", i);		   
// 		   }

// 		   //C2.- Idle to little
// 		   for(int i=currentBig + currentLittle; i<destBig+destLittle; i++){
// 		       sem_post(&_expControl._sems[i]);
// 		       //fprintf(stderr, "C2- idle to little [%d]\n", i);		   
// 		   }
// 		   //fprintf(stderr, "\n");
// 	       }

// 	       //D.- Update and release the Arch Struct
// 	       _expControl._archIdx --;
// 	       _expControl._mtxArch.release();
//    	       NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );	       
// 	   }
//        }

	//-----------------------------------------------------------------------
	// exp14                                                               //
	//-----------------------------------------------------------------------       
	void BotLev::exp14_ctr()
	{
	    ExperimentsCfg::nCoresUsed = sys.getNumThreads();
	    
	    _expControl.coreType = new E_coreType[ExperimentsCfg::nCoresUsed];
	    _expControl.sems     = new sem_t[ExperimentsCfg::nCoresUsed];
	    
	    for(int i=0; i<ExperimentsCfg::nCoresUsed; i++){
		sem_init(&_expControl.sems[i], 0, 1);
		_expControl.coreType[i] = LITTLE;
	    }
	   
	    _expControl.currentBudget = ExperimentsCfg::nCoresUsed * ExperimentsCfg::LittlePower;
	    NANOS_INSTRUMENT (_expControl.nBig=0);
	    NANOS_INSTRUMENT (_expControl.nFutureBig=0);
	    NANOS_INSTRUMENT (_expControl.nLITTLE=ExperimentsCfg::nCoresUsed);
	    NANOS_INSTRUMENT (_expControl.nSleeping=0);
	    

	    char *buf;
	    if((buf=getenv("MAX_BUDGET"))!=NULL)   ExperimentsCfg::MaxBudget   = atof(buf);
	    if((buf=getenv("LITTLE_POWER"))!=NULL) ExperimentsCfg::LittlePower = atof(buf);
	    if((buf=getenv("BIG_POWER"))!=NULL)    ExperimentsCfg::BigPower    = atof(buf);
	    if((buf=getenv("NCORES_USED"))!=NULL)  ExperimentsCfg::nCoresUsed  = atoi(buf);
	    if((buf=getenv("BIG_FREQ"))!=NULL)     ExperimentsCfg::bigFreq     = atoi(buf);
	    if((buf=getenv("LITTLE_FREQ"))!=NULL)  ExperimentsCfg::littleFreq  = atoi(buf);
	    

	    //_expControl.ImSleeping();
	    _expControl.IwantToPromote = -1;

	}
       
	void BotLev::exp14_dtr()
	{
	    for(int i=0; i<ExperimentsCfg::nCoresUsed; i++)
		sem_destroy(&_expControl.sems[i]);

	    delete[] _expControl.coreType;
	    delete[] _expControl.sems;
	}
       
	WD*  BotLev::exp14_atIdle ( BaseThread *thread, int numSteal )
	{
	    WorkDescriptor * wd;
	    TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
	    //unsigned int spins = BotLevCfg::numSpins;
	    int coreId = thread->runningOn()->getId();

	    int nCrit   = data._readyQueues[1].size(),
		nNoCrit = data._readyQueues[0].size();
	    
	    //A) Convert core if needed.
	    //	    NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MODARCH) );
	    this->exp14_modArch(coreId, nCrit, nNoCrit);
	    
	    //	    NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );

	    //B) get our own semaphore.
	    //+Maybe we are blocked by the previous call (modArch)
	    sem_wait(&_expControl.sems[coreId]);

	    //C) Execute
	    if( _expControl.coreType[coreId] == BIG ||
		_expControl.coreType[coreId] == FUTURE_BIG ) {        //Big core

		wd = data._readyQueues[1].pop_front( thread );
		// big steals from little queue
		if( wd == NULL )   wd = data._readyQueues[0].pop_front( thread );
	    }
	    else {    		                                    //Small core
		wd = data._readyQueues[0].pop_front( thread );
		if(!wd && BotLevCfg::steal)  wd = data._readyQueues[1].pop_front( thread );
	    }
	   
	    if( wd == NULL ) wd = data._readyQueues[2].pop_front( thread ); 
	    if( wd == NULL )
		sem_post(&_expControl.sems[coreId]);

	    return wd;
	}
       

	WD*  BotLev::exp14_atBeforeExit( BaseThread *thread, WD &wd, bool schedule)
	{
	    // Free the resource
	    int coreId = thread->runningOn()->getId();
	    sem_post(&_expControl.sems[coreId]);
	    
	    return NULL; 
	}
	
	void BotLev::exp14_queue ( BaseThread *thread, WD &wd )
	{
	    exp00_queue(thread,wd);
	}

	void BotLev::exp14_modArch(int coreId, int nCrit, int nNotCrit)
	{
#ifdef NANOS_INSTRUMENTATION_ENABLED
	    static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary();
	    static nanos_event_key_t critical_id = ID->getEventKey("n-crit");
	    static nanos_event_key_t non_critical_id = ID->getEventKey("n-noCrit");
	    static nanos_event_key_t n_big_id = ID->getEventKey("n-BIG");
	    static nanos_event_key_t n_futureBig_id = ID->getEventKey("n-FUTURE_BIG");
	    static nanos_event_key_t n_little_id = ID->getEventKey("n-LITTLE");
	    static nanos_event_key_t n_sleeping_id = ID->getEventKey("n-SLEEPING");
	    static nanos_event_key_t n_total_id = ID->getEventKey("n-TOTAL");
	    
	    nanos_event_key_t keys[7];
	    nanos_event_value_t values[7];
	    keys[0] = critical_id;      values[0] = nCrit;
	    keys[1] = non_critical_id;  values[1] = nNotCrit;

	    keys[2] = n_big_id;         values[2] = _expControl.nBig;
	    keys[3] = n_futureBig_id;   values[3] = _expControl.nFutureBig;
	    keys[4] = n_little_id;      values[4] = _expControl.nLITTLE;
	    keys[5] = n_sleeping_id;    values[5] = _expControl.nSleeping;
	    keys[6] = n_total_id;       values[6] = _expControl.nBig + _expControl.nFutureBig + _expControl.nLITTLE + _expControl.nSleeping;
	    
	    sys.getInstrumentation()->raisePointEvents(7, keys, values);
#endif

	    const float MAX_BUDGET   = ExperimentsCfg::MaxBudget;
	    const float LITTLE_POWER = ExperimentsCfg::LittlePower;
	    const float BIG_POWER    = ExperimentsCfg::BigPower;
	    const int NCORES_PER_APP = ExperimentsCfg::nCoresUsed;

	    
	    
    
	    //A.- Check if there's a core waiting to promote (sleep ourselves).
	    // A.1.- Maybe, we can promote it without going to sleep.
	    //   A.1.a Do we need still to promote ?
	    // A.2.- If not:
	    //   A.2.a.- Sleep ourselves
	    //   A.2.b.- Try to promote again
	    //B.- Convert to big:
	    // B.1.- We have enough free budget: try to promote ourselves
	    // B.2.- We need to wait another thread to finish first
	    //C.- Auto-convert to LITTLE.
	    // C.1.-  Wake up other threads?

	    
            //                                                                  //A.- There's another thread waiting to promote
	    if( _expControl.IwantToPromote != -1 &&
		_expControl.coreType[coreId] == LITTLE ) {
		
		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MOD_ARCH1) );
		LockBlock l(_expControl.mtxArch);

		if(_expControl.IwantToPromote != -1 &&                          // A.1.a Do we need still to promote ?
		   (nCrit == 0 && (nCrit + nNotCrit) >= NCORES_PER_APP )) {
		    
		    _expControl.coreType[_expControl.IwantToPromote] = LITTLE;
		    
		    NANOS_INSTRUMENT (_expControl.nLITTLE++);
		    NANOS_INSTRUMENT (_expControl.nFutureBig--);
		    NANOS_INSTRUMENT(					\
			fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
				nCrit, nNotCrit,			\
				coreId, _expControl.IwantToPromote, "De-promoting", \
 				_expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping));
		    _expControl.IwantToPromote = -1;
		}
		
		if(_expControl.IwantToPromote != -1){
		    
		    //LockBlock l2(_expControl.mtxBudget);
		    //                                                          // A.1.- Maybe, we can promote it without going to sleep.
		    if(_expControl.currentBudget + BIG_POWER - LITTLE_POWER <= MAX_BUDGET){
			_expControl.currentBudget += BIG_POWER - LITTLE_POWER;
			this->exp14_convertToBig(_expControl.IwantToPromote);
			_expControl.coreType[_expControl.IwantToPromote] = BIG;
			
			NANOS_INSTRUMENT (_expControl.nBig++);
			NANOS_INSTRUMENT (_expControl.nFutureBig--);
			NANOS_INSTRUMENT(				\
			    fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
				    nCrit, nNotCrit,			\
				    coreId, _expControl.IwantToPromote, "Promoting", \
				    _expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping));
			_expControl.IwantToPromote = -1;
		    } else if(coreId !=0) {                                                    // A.2.- If not:
			sem_wait(&_expControl.sems[coreId]);                    //  A.2.a.- Sleep ourselves (if we're not the thread 0!!!
			_expControl.ImSleeping.push_back(coreId);
			_expControl.currentBudget -= LITTLE_POWER;
			_expControl.coreType[coreId] = SLEEP;
			
			NANOS_INSTRUMENT (_expControl.nLITTLE--);
			NANOS_INSTRUMENT (_expControl.nSleeping++);
			
			NANOS_INSTRUMENT(				\
			    fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
				    nCrit, nNotCrit,			\
				    coreId, coreId, "SLEEP",		\
				    _expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping));
			//                                                      //   A.2.b.- Try to promote again
			if(_expControl.currentBudget + BIG_POWER - LITTLE_POWER <= MAX_BUDGET) {
			    _expControl.currentBudget += BIG_POWER - LITTLE_POWER;
			    this->exp14_convertToBig(_expControl.IwantToPromote);
			    _expControl.coreType[_expControl.IwantToPromote] = BIG;
			    
			    NANOS_INSTRUMENT (_expControl.nBig++);
			    NANOS_INSTRUMENT (_expControl.nFutureBig--);
			    NANOS_INSTRUMENT(				\
				fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
					nCrit, nNotCrit,		\
					coreId, _expControl.IwantToPromote, "Promoting", \
					_expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping));
			    _expControl.IwantToPromote = -1;
			}
		    }
		    return;
		}
		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
	    }
	    
	    
	    //                                                                  
	    if( (nCrit + nNotCrit > 0) &&
		(nCrit > 0 || (nCrit + nNotCrit) < NCORES_PER_APP) &&           //B.- There're critical tasks on the queue
		_expControl.IwantToPromote == -1 ) {                            //There's not any th waiting to promote:
		if( _expControl.coreType[coreId] == BIG ||
		    _expControl.coreType[coreId] == FUTURE_BIG )
		    return;
		
		LockBlock l2(_expControl.mtxArch);
		if(_expControl.IwantToPromote != -1)   return;		    

		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MOD_ARCH2) );
		//LockBlock l (_expControl.mtxBudget);                        
		{                                                               //B.1.- There's enough free budget
		    if( _expControl.currentBudget + BIG_POWER - LITTLE_POWER <= MAX_BUDGET){ 
			_expControl.currentBudget += (BIG_POWER - LITTLE_POWER);
			
			this->exp14_convertToBig(coreId);
			_expControl.coreType[coreId] = BIG;
			
			NANOS_INSTRUMENT (_expControl.nBig++);
			NANOS_INSTRUMENT (_expControl.nLITTLE--);
			
			NANOS_INSTRUMENT(				\
 			    fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
 				    nCrit, nNotCrit,			\
 				    coreId, coreId, "BIG",		\
 				    _expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping); )
		     } else {                                                    //B.2.- we need to wait other threads to finish
			_expControl.coreType[coreId] = FUTURE_BIG;		    
			_expControl.IwantToPromote = coreId;
			
			NANOS_INSTRUMENT (_expControl.nFutureBig++);
			NANOS_INSTRUMENT (_expControl.nLITTLE--);
			
			NANOS_INSTRUMENT (				\
			    fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
				    nCrit, nNotCrit,			\
				    coreId, coreId, "FUTURE_BIG",	\
				    _expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping);)
		     }
		}
		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
		return;
	    }


	    
	    if( _expControl.coreType[coreId] != LITTLE &&
		( nCrit == 0 || (nCrit + nNotCrit) >= NCORES_PER_APP )) {       // C.- Auto-convert to LITTLE.
		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MOD_ARCH3) );
		LockBlock l2 (_expControl.mtxArch);
		
		if(_expControl.coreType[coreId] == BIG){
		    //LockBlock l (_expControl.mtxBudget);
		    _expControl.currentBudget += (LITTLE_POWER - BIG_POWER);
		    this->exp14_convertToLITTLE(coreId);
		    _expControl.coreType[coreId] = LITTLE;				    
		    
		    NANOS_INSTRUMENT (_expControl.nBig-- );
		    NANOS_INSTRUMENT (_expControl.nLITTLE++);
		    NANOS_INSTRUMENT(					\
			fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
				nCrit, nNotCrit,			\
				coreId, coreId, "LITTLE",		\
				_expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping));		    
			
		} else if(_expControl.coreType[coreId] == FUTURE_BIG){
		    _expControl.IwantToPromote = -1;
		    _expControl.coreType[coreId] = LITTLE;
		    
		    NANOS_INSTRUMENT ( _expControl.nFutureBig-- );
		    NANOS_INSTRUMENT (_expControl.nLITTLE++);
		    NANOS_INSTRUMENT(					\
			fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
				nCrit, nNotCrit,			\
				coreId, coreId, "LITTLE",		\
				_expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping));
		}		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
	    }	    
	    { 	                                                                // C.1.-wake up other threads
		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MOD_ARCH4) );
		LockBlock l2(_expControl.mtxArch);
		//LockBlock l(_expControl.mtxBudget);
		while(_expControl.IwantToPromote == -1 && //none wants to promote
		      _expControl.ImSleeping.size() > 0 &&
		      _expControl.currentBudget + LITTLE_POWER <= MAX_BUDGET)
		{
		    _expControl.currentBudget += LITTLE_POWER;
		    std::vector<int>::iterator it = _expControl.ImSleeping.begin();
		    sem_post(&_expControl.sems[*it]);
		    _expControl.coreType[*it] = LITTLE;
		    NANOS_INSTRUMENT( int cc = *it );
		    
		    _expControl.ImSleeping.erase(it);

		    NANOS_INSTRUMENT (_expControl.nLITTLE++);
		    NANOS_INSTRUMENT (_expControl.nSleeping--);

		    NANOS_INSTRUMENT (					\
			fprintf(stderr, "[ %d - %d ],[(%d) %d -> %s ], little: %d, big: %d, fbig: %d, sleeping: %d//n", \
				nCrit, nNotCrit,			\
				coreId, cc, "LITTLE",		 \
				_expControl.nLITTLE, _expControl.nBig, _expControl.nFutureBig, _expControl.nSleeping);)		}
		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
		return;
	    }
	}

	void BotLev::exp14_convertToBig(int coreId)
	{
	    long long currFreq = cpufreq_get_freq_kernel(coreId);	 
	    if(currFreq != ExperimentsCfg::bigFreq) {
		cpufreq_modify_policy_max(coreId, ExperimentsCfg::bigFreq);
		cpufreq_modify_policy_min(coreId, ExperimentsCfg::bigFreq);
		cpufreq_set_frequency(coreId, ExperimentsCfg::bigFreq);
	    }	    
	}
	
	void BotLev::exp14_convertToLITTLE(int coreId)
	{
	    long long currFreq = cpufreq_get_freq_kernel(coreId);	 
	    if(currFreq != ExperimentsCfg::littleFreq) {
		cpufreq_modify_policy_max(coreId, ExperimentsCfg::littleFreq);
		cpufreq_modify_policy_min(coreId, ExperimentsCfg::littleFreq);
		cpufreq_set_frequency(coreId, ExperimentsCfg::littleFreq);
	    }	
	}
	
	void BotLev::exp14_atShutdown()
	{
	    for(int i=0; i<ExperimentsCfg::nCoresUsed; ++i)
		sem_post(&_expControl.sems[i]);
	}



//	// -----------------------------------------------------------------------
//	// exp15                                                               //
//	// -----------------------------------------------------------------------       
// 	void BotLev::exp15_ctr()
// 	{
// 	    ExperimentsCfg::nCoresUsed = sys.getNumThreads();
// 	    _expControl.sems     = new sem_t[ExperimentsCfg::nCoresUsed];
// 	    _expControl.coreType = new E_coreType[ExperimentsCfg::nCoresUsed];
	    
// 	    for(int i=0; i<ExperimentsCfg::nCoresUsed; i++){
// 		sem_init(&_expControl.sems[i], 0, 0); //all sems are init to 0!!!
// 		_expControl.coreType[i] = SLEEP;
// 	    }


// //	    NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MOD_ARCH1) );
// 	    {
// 		SS_initLibrary();
// 		SS_startExecution(&_expControl.mask);

// 		exp15_setMask(_expControl.mask);
// 	    }
// //	    NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
// 	}
       
// 	void BotLev::exp15_dtr()
// 	{
// 	    for(int i=0; i<ExperimentsCfg::nCoresUsed; i++)
// 		sem_destroy(&_expControl.sems[i]);

// 	    delete[] _expControl.sems;

// 	    SS_finishLibrary();
// 	}
       
// 	WD*  BotLev::exp15_atIdle ( BaseThread *thread, int numSteal )
// 	{
// 	    WorkDescriptor * wd;
// 	    TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
// 	    //unsigned int spins = BotLevCfg::numSpins;
// 	    int coreId = thread->runningOn()->getId();



// 	    int nCrit    = data._readyQueues[0].size();
// 	    int nNotCrit = data._readyQueues[1].size();
// 	    int nOthers  = data._readyQueues[2].size();

// 	    if( nCrit == 0 && nNotCrit == 0 && nOthers == 0)
// 		return NULL;
	    
	    
// 	    if(_expControl.ss_com.tryAcquire()){
// 		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_MOD_ARCH2) );

// 	    	SS_updateNtasks(nCrit, nNotCrit);

// 	    	SS_askForNewCpuMask(&_expControl.mask);
// 	    	exp15_setMask(_expControl.mask);
		
// 	    	_expControl.ss_com.release();
// 		NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
// 	    }

	    
// 	    //B) get our own semaphore.
// 	    //+Maybe we are blocked by another call
// 	    sem_wait(&_expControl.sems[coreId]);

// 	    //BIG
// 	    wd = data._readyQueues[1].pop_front( thread );

// 	    //If not, LITTLE
// 	    if( wd == NULL )   wd = data._readyQueues[0].pop_front( thread );
	    
// 	    //If not, others
// 	    if( wd == NULL ) wd = data._readyQueues[2].pop_front( thread ); 
// 	    if( wd == NULL )
// 		sem_post(&_expControl.sems[coreId]);

// 	    return wd;
// 	}
       

// 	WD*  BotLev::exp15_atBeforeExit( BaseThread *thread, WD &wd, bool schedule)
// 	{
// 	    // Free the resource
// 	    int coreId = thread->runningOn()->getId();
// 	    sem_post(&_expControl.sems[coreId]);
	    
// 	    return NULL; 
// 	}
	
// 	void BotLev::exp15_queue ( BaseThread *thread, WD &wd )
// 	{
// 	    exp00_queue(thread,wd);
// 	}
	
// 	void BotLev::exp15_atShutdown()
// 	{
// 	    for(int i=0; i<ExperimentsCfg::nCoresUsed; ++i)
// 		sem_post(&_expControl.sems[i]);
// 	}
	
	    
// 	void BotLev::exp15_setMask(dynamic_bitset &mask)
// 	{
// 	    for(int i=0; i<ExperimentsCfg::nCoresUsed; ++i){
		
// 		if(mask.test(i)){
// 		    if(_expControl.coreType[i] == SLEEP){
// 			sem_post(&_expControl.sems[i]);
// 			_expControl.coreType[i] = RUNNING;
// 		    }
// 		} else {
// 		    if(_expControl.coreType[i] != SLEEP){
// 			sem_wait(&_expControl.sems[i]);
// 			_expControl.coreType[i] = SLEEP;

// 			SS_confirmCoreFree(i);
// 		    }
// 		}
// 	    }
// 	}
				   

        //=======================================================================
	//=======================================================================


	
       
	/**-------------------------**/
	/** ____fgcsSchedPlugin____ **/
	/**-------------------------**/
	fgcsSchedPlugin::fgcsSchedPlugin() : Plugin( "Distributed Breadth-First scheduling Plugin",1 ) {}       
	void fgcsSchedPlugin::config ( Config &config_ ) {
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


	    config_.registerConfigOption ( "changeFreq", new Config::IntegerVar( BotLevCfg::experiment_type ), "Defines which experiment we shouls launch" );
	    config_.registerArgOption ( "changeFreq", "changeFreq" );
	    config_.registerEnvOption ( "changeFreq", "NX_CHANGEFREQ" );



	    config_.registerConfigOption ( "changeFreq", new Config::IntegerVar( BotLevCfg::experiment_type ), "Defines which experiment we shouls launch" );
	    config_.registerArgOption ( "changeFreq", "changeFreq" );
	    config_.registerEnvOption ( "changeFreq", "NX_CHANGEFREQ" );



	    // config_.registerConfigOption ( "maxBudget", new Config::PositiveVar( ExperimentsCfg::MaxBudget ), "Max power budget for the experiments" );
	    // config_.registerArgOption ( "maxBudget", "maxBudget" );
	    // config_.registerEnvOption ( "maxBudget", "maxBudget" );

	}

	void fgcsSchedPlugin::init() {
	    sys.setDefaultSchedulePolicy(new BotLev());
	}
       
    } //ext namespace
} //nanos namespace


DECLARE_PLUGIN("fgcs",nanos::ext::fgcsSchedPlugin);



















// WD * BotLev::atBeforeExit(BaseThread *thread, WD &wd, bool schedule) {
//    NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
//    NANOS_INSTRUMENT ( static nanos_event_key_t criticalityKey = ID->getEventKey("wd-criticality"); )
//    NANOS_INSTRUMENT ( nanos_event_value_t wd_criticality; )
//    NANOS_INSTRUMENT ( WDData & wddata = *dynamic_cast<WDData*>( wd.getSchedulerData() ); )
//    NANOS_INSTRUMENT ( wd_criticality = wddata.getCriticality(); )
//    NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseBurstEvent ( criticalityKey, wd_criticality ); )
//    return NULL;
// }


