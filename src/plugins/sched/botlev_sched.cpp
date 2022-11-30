#include <unistd.h>
#include "schedule.hpp"
#include "wddeque.hpp"
#include "plugin.hpp"
#include "system.hpp"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <unistd.h>
#include <signal.h>

#include <cpufreq.h>

#define MAX(a,b) ((a > b) ? (a) : (b))

#define SUPPORT_UPDATE 1

#ifdef ODROID
#define CORE_BIG  4
#define CORE_LITTLE 0
#define N_CORE_BIG 4
#define N_CORE_LITTLE 4

#elif JUNO
#define CORE_BIG  4
#define CORE_LITTLE 0
#define N_CORE_BIG 2
#define N_CORE_LITTLE 4

#elif KELVIN
#define CORE_BIG  2
#define CORE_LITTLE 0
#define N_CORE_BIG 2
#define N_CORE_LITTLE 2

#elif VOLTA1
#define CORE_BIG 4
#define CORE_LITTLE 0
#define N_CORE_BIG 4
#define N_CORE_LITTLE 4

#else
ERROR: FALTA POR DEFINIR LA MAQUINA DESTINO (-D...);
#endif


namespace nanos {
    namespace ext {
	struct BotLevCfg {
	public:
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

	// Initialise default values
	int BotLevCfg::updateFreq = 0;
	int BotLevCfg::numSpins = 300;
	int BotLevCfg::hpFrom = 0;
	int BotLevCfg::hpTo = 0;
	int BotLevCfg::hpSingle = 0;
	int BotLevCfg::steal = 0;
	int BotLevCfg::maxBL = 1;
	int BotLevCfg::strict = 0;
	int BotLevCfg::taskNumber = 0;
	NANOS_INSTRUMENT( int BotLevCfg::numCritical; )


	/*@ -- @*/
	/** Struct para guardar el estado actual de los cores, freq, ... **/
	struct FreqCfg {
	public:
	    static int   changeFreq;    // Env option. activate/deactivate

	    static int   maxIdxLittle;
	    static int   maxIdxBig;
#ifdef ODROID
	    static int   littleFreq[6];
	    static int   bigFreq[6];
#elif JUNO
	    static int   littleFreq[5];
	    static int   bigFreq[5];
#elif KELVIN
	    static int   littleFreq[6];
	    static int   bigFreq[6];
#elif VOLTA1
	    static int   littleFreq[3];
	    static int   bigFreq[3];
#endif
	    static int   currLittleIdx; //Current little freq.
	    static int   currBigIdx;    //Current bif freq.

	    static Lock  *_lock;          //! Structure lock

	    static int   maxQueueSize;
	    ~FreqCfg();

	    //Metodos auxiliares para cambiar la frecuencia
	    static unsigned long changeBigFreq(int idxFreq);
	    static unsigned long changeLittleFreq(int idxFreq);
	};

	FreqCfg::~FreqCfg(){
	    delete _lock;
	}
    
	int FreqCfg::changeFreq    = 1; // default = active (Raw)
	int FreqCfg::maxQueueSize  = 0;
	int FreqCfg::currLittleIdx = 0;
	int FreqCfg::currBigIdx    = 0;
	Lock *FreqCfg::_lock       = new Lock();

#ifdef ODROID
	int FreqCfg::maxIdxLittle  = 5;
	int FreqCfg::maxIdxBig     = 5;
	int FreqCfg::littleFreq[6] = {1300000, 1200000, 1100000, 1000000, 900000, 800000}; // A7
	int FreqCfg::bigFreq[6]    = {1300000, 1200000, 1100000, 1000000, 900000, 800000}; // A15
#elif JUNO
	int FreqCfg::maxIdxLittle  = 4;
	int FreqCfg::maxIdxBig     = 4;
	int FreqCfg::littleFreq[5] = {850000,  775000, 700000, 575000, 450000}; // A53
	int FreqCfg::bigFreq[5]    = {1100000, 950000, 800000, 625000, 450000}; // A57
#elif KELVIN
	int FreqCfg::maxIdxLittle  = 5;
	int FreqCfg::maxIdxBig     = 5;
	int FreqCfg::littleFreq[6] = {3500000, 3000000, 2500000, 2000000, 1500000, 1000000};
	int FreqCfg::bigFreq[6]    = {3500000, 3000000, 2500000, 2000000, 1500000, 1000000};
#elif VOLTA1
	int FreqCfg::maxIdxLittle  = 2;
	int FreqCfg::maxIdxBig     = 2;                                  
	int FreqCfg::littleFreq[3] = {1200000, 1100000, 1000000}; //0-3
      int FreqCfg::bigFreq[3]    = //{2001000, 2000000, 1800000}; //4-7
	                             {2000000, 1900000, 1800000}; //4-7
#else
    ERROR: NO SE ENCUENTRA DEFINIDA LA PLACA; //Error en tiempo de compilacion
#endif

	unsigned long FreqCfg::changeBigFreq(int idxFreq){
	    long long currFreq = cpufreq_get_freq_kernel(CORE_BIG);

	    if( currFreq == FreqCfg::bigFreq[idxFreq] )
		return currFreq;
	
	

	    FreqCfg::_lock->acquire();
	    unsigned long ret;
      
	    for(int i=0; i<N_CORE_BIG; i++){
		currFreq = cpufreq_get_freq_kernel(CORE_BIG+i);

		if(currFreq > FreqCfg::currBigIdx) { //disminuimos
		    cpufreq_modify_policy_min(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_modify_policy_max(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);

		}
		else { //aumentamos
		    cpufreq_modify_policy_max(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_modify_policy_min(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);

		}
	    }

	    FreqCfg::currBigIdx = idxFreq;
	    ret = cpufreq_get_freq_kernel(CORE_BIG);      

	    FreqCfg::_lock->release();
      
	    return ret;
	}

	unsigned long FreqCfg::changeLittleFreq(int idxFreq){
	    long long currFreq = cpufreq_get_freq_kernel(CORE_LITTLE);

	    if( currFreq == FreqCfg::littleFreq[idxFreq])
		return currFreq;

	
	    FreqCfg::_lock->acquire();
	    unsigned long ret;


	    for(int i=0; i<N_CORE_LITTLE; i++){
		currFreq = cpufreq_get_freq_kernel(CORE_LITTLE + i);
		if( currFreq > FreqCfg::currLittleIdx){ //disminuimos
		    cpufreq_modify_policy_min(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_modify_policy_max(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
	    
		}
		else { //aumentamos
		    cpufreq_modify_policy_max(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_modify_policy_min(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);

		}
	    }

	    FreqCfg::currLittleIdx = idxFreq;
	    ret = cpufreq_get_freq_kernel(CORE_LITTLE);
      
	    FreqCfg::_lock->release();
 
	    return ret;      
	}
	/*@ -- @*/


	//<--------------------------------------------------------------------------->
	/* ESTRUCTURA PARA APAGAR/ENCENDER CORES */
	struct CoresCfg {
	public:
	    static int  level; //Nivel de las colas a las que apagar/encender los cores
	    static bool apagadoLittle;
	    static bool apagadoBig;
	    static int  coresBig [N_CORE_BIG];
	    static int  coresLittle [N_CORE_LITTLE];
	    static bool coresApagados [(N_CORE_BIG + N_CORE_LITTLE)];
      
	    static bool puedoEjecutar(int coreId);
	    static void enciendeBig();
	    static void enciendeLittle();

	    static void apagaLittle();
	    static void apagaBig();
	    CoresCfg();
	    ~CoresCfg();
	private:
	    static void apagaCore(int coreId);
	    static void enciendeCore(int coreId);
      
	    static Lock *_lock [(N_CORE_BIG + N_CORE_LITTLE)];
      
	};
      
	int CoresCfg::level           = 30;
	bool CoresCfg::apagadoLittle  = false;
	bool CoresCfg::apagadoBig     = false;
#ifdef ODROID
	int CoresCfg::coresBig [N_CORE_BIG]       = {4, 5, 6, 7};
	int CoresCfg::coresLittle [N_CORE_LITTLE] = {0, 1, 2, 3};
	bool CoresCfg::coresApagados [(N_CORE_BIG + N_CORE_LITTLE)] = {false, false, false, false,
								       false, false, false, false};
	Lock* CoresCfg::_lock [(N_CORE_BIG + N_CORE_LITTLE)] = {new Lock(), new Lock(), new Lock(), new Lock(),
								new Lock(), new Lock(), new Lock(), new Lock()};
#elif JUNO
	int CoresCfg::coresBig [N_CORE_BIG]       = {4, 5};
	int CoresCfg::coresLittle [N_CORE_LITTLE] = {0, 1, 2, 3};
	bool CoresCfg::coresApagados [(N_CORE_BIG + N_CORE_LITTLE)] = {false, false, false, false,
								       false, false};
	Lock* CoresCfg::_lock [(N_CORE_BIG + N_CORE_LITTLE)] = {new Lock(), new Lock(), new Lock(), new Lock(),
								new Lock(), new Lock()};
#elif KELVIN
	int CoresCfg::coresBig [N_CORE_BIG]       = {2, 3};
	int CoresCfg::coresLittle [N_CORE_LITTLE] = {0, 1};
	bool CoresCfg::coresApagados [(N_CORE_BIG + N_CORE_LITTLE)] = {false, false,
								       false, false};

	Lock* CoresCfg::_lock [(N_CORE_BIG + N_CORE_LITTLE)] = {new Lock(), new Lock()
								new Lock(), new Lock()};
#elif VOLTA1
	int CoresCfg::coresBig [N_CORE_BIG]       = {4,5,6,7};
	int CoresCfg::coresLittle [N_CORE_LITTLE] = {0,1,2,3};
	bool CoresCfg::coresApagados [(N_CORE_BIG + N_CORE_LITTLE)] = {false, false, false, false,
								       false, false, false, false};
	Lock* CoresCfg::_lock [(N_CORE_BIG + N_CORE_LITTLE)] = {new Lock(), new Lock(), new Lock(), new Lock(),
								new Lock(), new Lock(), new Lock(), new Lock()};
#endif    
    
	bool CoresCfg::puedoEjecutar(int coreId){
	    bool ret = true;

	    //Comprobar cores little
	    if(apagadoLittle){
		for(int i=0; i<N_CORE_LITTLE; i++){
		    if(coresLittle[i] == coreId){
			if(!coresApagados[coreId]){
//		      FreqCfg::changeLittleFreq(FreqCfg::maxIdxLittle);
			    apagaCore(coreId);
			}
			ret = false;
			break;
		    }
		}
	    }
    
	    //Comprobar que es un core big
	    if(apagadoBig){
		for(int i=0; i<N_CORE_BIG; i++){
		    if(coresBig[i] == coreId){
			if(!coresApagados[coreId]){
//		      FreqCfg::changeBigFreq(FreqCfg::maxIdxBig);
			    apagaCore(coreId);
			}
			ret = false;
			break;
		    }
		}
	    }
      

#ifdef NANOS_INSTRUMENTATION_ENABLED
	    nanos_event_key_t clave[1];
	    nanos_event_value_t valor[1];
	    clave[0] = 5555;
	    if( !ret) valor[0] = 0;
	    else valor[0] = 10;
	    NANOS_INSTRUMENT(sys.getInstrumentation()->raisePointEvents(1, clave, valor); )
#endif
	  
		return ret;
	}

	void CoresCfg::apagaCore(int coreId){
	    //printf("Voy a apagar el core %d\n", coreId);
	    coresApagados[coreId] = true;
      
	
	    if(FreqCfg::changeFreq==7 || FreqCfg::changeFreq==8){
		//printf("Voy a apagar el core %d\n", coreId);
		return;
	    }

      
	    char path[50];
	    sprintf(path, "/sys/devices/system/cpu/cpu%d/online", coreId);

	    FILE* f;
	    int d, a;

	    (CoresCfg::_lock[coreId])->acquire();
	    {
		f = fopen(path, "r+");
		a = fscanf(f, "%d", &d);
	
		if(d!=0 && a!=0){
		    fseek(f, 0, SEEK_SET);
		    fprintf(f, "0"); //Apagado
		}
	
		fclose(f);
	    }
	    (CoresCfg::_lock[coreId])->release();

	    printf("He acabado de apagar el core %d\n", coreId);
	}

	void CoresCfg::enciendeCore(int coreId){
	    //      printf("Voy a encender el core %d\n", coreId);
	    coresApagados[coreId] = false;

      
	    if(FreqCfg::changeFreq==7 || FreqCfg::changeFreq==8){
		return;
	    }


	    char path[50];
	    sprintf(path, "/sys/devices/system/cpu/cpu%d/online", coreId);

	    FILE* f;
	    int d, a;

	    (CoresCfg::_lock[coreId])->acquire();
	    {
		f = fopen(path, "r+");
		a = fscanf(f, "%d", &d);
	
		if(d!=1 && a!=0){
		    fseek(f, 0, SEEK_SET);
		    fprintf(f, "1"); //Encendido
		}
	
		fclose(f);
		coresApagados[coreId] = false;
	    }
	    (CoresCfg::_lock[coreId])->release();
	    printf("He acabado de encender el core %d\n", coreId);
	}

    
	void CoresCfg::enciendeBig(){
	    if(!apagadoBig) return;

	    FreqCfg::changeBigFreq(0); //Freq al máximo
      
	    for (int i=0; i<N_CORE_BIG; i++)
		enciendeCore(coresBig[i]);
	    apagadoBig = false;
	}

    
	void CoresCfg::enciendeLittle(){
	    if( !apagadoLittle ) return;

	    FreqCfg::changeLittleFreq(0); //Freq al máximo
      
	    for(int i=0; i<N_CORE_LITTLE; i++)
		enciendeCore(coresLittle[i]);
	    apagadoLittle = false;
	}

    
	void CoresCfg::apagaLittle(){
	    apagadoLittle=true;
	}

    
	void CoresCfg::apagaBig(){
	    apagadoBig = true;
	}

	//<--------------------------------------------------------------------------->
    

	class BotLevDOData : public DOSchedulerData {
	    /* Helper class for the computation of bottom levels 
	       One instance of this class is stored inside each dependableObject.  */
	public:
	    typedef std::set<BotLevDOData *> predecessors_t;

	private:
	    int               _taskNumber;    //! \todo Keep task number just for debugging, remove?
	    int               _botLevel;      //! Bottom Level Value
	    bool              _isReady;       //! Is the task ready
	    Lock              _lock;          //! Structure lock
	    WD               *_wd;            //! Related WorkDescriptor
	    predecessors_t    _predecessors;  //! List of BotLevDOData predecessors 
	    short             _isCritical;    //! Is the task critical -- needed for reordering the ready queues

	public:
	    BotLevDOData(int tnum, int blev) : _taskNumber(tnum), _botLevel(blev), _isReady(false), _lock(), _wd(NULL), _predecessors() { }
	    ~BotLevDOData() { }
	    void reset () { _wd = NULL; }

	    int getTaskNumber() const { return _taskNumber; }

	    int getBotLevel() const { return _botLevel; }
	    bool setBotLevel( int blev )
		{
		    if ( blev > _botLevel ) {
			{
			    LockBlock lock1( _lock );
			    _botLevel = blev;
			}
			return true;
		    }
		    return false;
		}

	    bool getReady() const { return _isReady; }
	    void setReady() { _isReady = true; }

	    void setCriticality( short c ) { _isCritical = c; }
	    short getCriticality()         { return _isCritical; }

	    WD* getWorkDescriptor() const { return _wd; }
	    void setWorkDescriptor(WD *wd) { _wd = wd; }

	    std::set<BotLevDOData *> *getPredecessors() { return &_predecessors; }
	    void addPredecessor(BotLevDOData *predDod) { _predecessors.insert(predDod); }
	};

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

	    struct TeamData : public ScheduleTeamData
	    {
		/*! queues of ready tasks to be executed */
		WDPriorityQueue<> *_readyQueues;

		/*@ -- @*/
		Lock* _lockfout;
		FILE* _fout;

		TeamData () : ScheduleTeamData() {
		    _readyQueues = NEW WDPriorityQueue<>[3];

		    _lockfout = new Lock();
		    _fout = fopen("MUESTRA.m", "w");
		    fprintf(_fout, "data0=[];\ndata1=[];\nfreq=[];\n");
		    fflush(_fout);
		    /*@ -- @*/
		}

		virtual ~TeamData () {
		    delete[] _readyQueues;
		    /*@ -- @*/

		    NANOS_INSTRUMENT( fprintf(_fout,"Num Critical Detected ===> %d\n",  BotLevCfg::numCritical); )

			fclose(_fout);
		    delete _lockfout;
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
	    // constructor
	    BotLev() : SchedulePolicy ( "BotLev" ) {
		sys.setPredecessorLists(true);
		_currMax = _maxBotLev = BotLevCfg::maxBL;
		NANOS_INSTRUMENT( BotLevCfg::numCritical = 0; )
		    }

	    // destructor
	    virtual ~BotLev() {
	    }

	    virtual size_t getTeamDataSize () const { return sizeof(TeamData); }
	    virtual size_t getThreadDataSize () const { return 0; }

	    virtual ScheduleTeamData * createTeamData () {
		return NEW TeamData();
	    }

	    virtual ScheduleThreadData * createThreadData () {
		return 0;
	    }

	    struct WDData : public ScheduleWDData {
		int _criticality;

		void setCriticality( int cr ) { _criticality = cr; }
		int getCriticality ( ) { return _criticality; }
		WDData () : _criticality( 0 ) {}
		virtual ~WDData() {}
	    };


	    //---------------------------------------------------------------------------
	    // QUEUE
	    //---------------------------------------------------------------------------

	    /*!
	     *  \brief Enqueue a work descriptor in the readyQueue of the passed thread
	     *  \param thread pointer to the thread to which readyQueue the task must be appended
	     *  \param wd a reference to the work descriptor to be enqueued
	     *  \sa ThreadData, WD and BaseThread
	     */
	    virtual void queue ( BaseThread *thread, WD &wd ) {
#ifdef NANOS_INSTRUMENTATION_ENABLED
		int criticality;
#endif
		TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();


		/*@ -- @*/
		//	int a0, a1;
		//a0 = data._readyQueues[0].size();
		//a1 = data._readyQueues[1].size();
		/*@ -- @*/


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

#endif // NANOS_INSTRUMENTATION_ENABLED

		/*@ --------------------------------- @*/
		int b0, b1;
		b0 = data._readyQueues[0].size();
		b1 = data._readyQueues[1].size();
	
		int f =  changeFrequency(b0, b1) / 1000; //small - big
	
		//-- Salida por paraver --//
#ifdef NANOS_INSTRUMENTATION_ENABLED
		NANOS_INSTRUMENT(nanos_event_key_t clave[3];)
		    NANOS_INSTRUMENT(nanos_event_value_t valor[3];)
		    NANOS_INSTRUMENT(clave[0] = 6666; clave[1] = 6667; clave[2] = 6668;)
		    NANOS_INSTRUMENT(valor[0] = b0;   valor[1] = b1;   valor[2] = f;)
	
		    NANOS_INSTRUMENT(sys.getInstrumentation()->raisePointEvents(3, clave, valor); )
#else
		    //-- Salida por fichero de muestras --//
		    data._lockfout->acquire();
		fprintf(data._fout, "data0 = [ data0; %u %d ];\n", (unsigned)time(NULL), b0);
		fprintf(data._fout, "data1 = [ data1; %u %d ];\n", (unsigned)time(NULL), b1);
		fprintf(data._fout, "freq  = [ freq;  %u %d ];\n", (unsigned)time(NULL), f);
		// fflush(data._fout);
		data._lockfout->release();
#endif
		/*@ --------------------------------- @*/

		return;
	    }


	    void updateBottomLevels( BotLevDOData *dodata, int botLev ) {
		std::vector<BotLevDOData *> stack;
		//! Set the bottom level, and add to the stack
		dodata->setBotLevel(botLev);
		stack.push_back(dodata);
		//Task is ready so, there are no predecessors to be updated
		if( dodata->getReady() )  return;
		// A depth first traversal
		while ( !stack.empty() ) {

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

      

	    void atCreate ( DependableObject &depObj ) {
		NANOS_INSTRUMENT( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
		    NANOS_INSTRUMENT( static nanos_event_key_t wd_atCreate  = ID->getEventKey("wd-atCreate"); )
		    NANOS_INSTRUMENT( WD * relObj = (WD*)depObj.getRelatedObject(); )
		    NANOS_INSTRUMENT( unsigned wd_id = ((WD *)relObj)->getId( ); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( wd_atCreate, wd_id ); )
              
		    NANOS_INSTRUMENT( static nanos_event_key_t blev_overheads  = ID->getEventKey("blev-overheads"); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( blev_overheads, 1 ); ) 
		    NANOS_INSTRUMENT( static nanos_event_key_t blev_overheads_br  = ID->getEventKey("blev-overheads-breakdown"); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( blev_overheads_br, NANOS_SCHED_BLEV_ATCREATE ); )

		    //! Creating DO scheduler data
		    BotLevDOData *dodata = new BotLevDOData(++BotLevCfg::taskNumber, 0);
		depObj.setSchedulerData( (DOSchedulerData*) dodata );

		//std::set<DependableObject *> predecessors;
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
		NANOS_INSTRUMENT( sys.getInstrumentation()->raiseCloseBurstEvent ( wd_atCreate, wd_id ); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseCloseBurstEvent ( blev_overheads, 1 ); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( blev_overheads_br, NANOS_SCHED_BLEV_ATCREATE ); )
		    }

      
	    /*!
	     *  \brief Function called when a new task must be created: the new created task
	     *          is directly queued (Breadth-First policy)
	     *  \param thread pointer to the thread to which belongs the new task
	     *  \param wd a reference to the work descriptor of the new task
	     *  \sa WD and BaseThread
	     */
	    virtual WD * atSubmit ( BaseThread *thread, WD &newWD ) {
		NANOS_INSTRUMENT( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
		    NANOS_INSTRUMENT( static nanos_event_key_t blev_overheads  = ID->getEventKey("blev-overheads"); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( blev_overheads, 1 ); )
		    NANOS_INSTRUMENT( static nanos_event_key_t blev_overheads_br  = ID->getEventKey("blev-overheads-breakdown"); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( blev_overheads_br, NANOS_SCHED_BLEV_ATSUBMIT ); )

		    /* Intentamos cambiar el main thread al 4 */
		    if(thread->runningOn()->getId() == 4 && ! thread->isMainThread())
			thread->setMainThread(true);
		/* Intentamos que el 0 no lo sea */
		if(thread->runningOn()->getId() == 0 && thread->isMainThread())
		    thread->setMainThread(false);

	
		queue(thread,newWD);

		NANOS_INSTRUMENT( sys.getInstrumentation()->raiseCloseBurstEvent ( blev_overheads, 1 ); )
		    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseCloseBurstEvent ( blev_overheads_br, NANOS_SCHED_BLEV_ATSUBMIT ); )

		    return 0;
	    }

	    virtual WD *atIdle( BaseThread *thread, int numSteal );

	    //---------------------------------------------------------------------------
	    unsigned long changeFrequency(int littleQueueSize, int bigQueueSize);
	    unsigned long changeFrequencyRaw(int littleQueueSize, int bigQueueSize);
	    unsigned long changeFrequencySoft(int littleQueueSize, int bigQueueSize);
	    unsigned long changeFrequencySoft2Levels(int littleQueueSize, int bigQueueSize);
	    unsigned long changeFrequencySoftBig(int littleQueueSize, int bigQueueSize);
      
	    unsigned long changeLittleFreq(int idxFreq);
	    unsigned long changeBigFreq(int idxFreq);


	    unsigned long apagaCoresLittle(int littleQueueSize, int bigQueueSize);
	    unsigned long apagaCoresBig(int littleQueueSize, int bigQueueSize);

	    unsigned long noAsignarCoresLittle(int littleQueueSize, int bigQueueSize);
	    unsigned long noAsignarCoresBig(int littleQueueSize, int bigQueueSize);
	    //---------------------------------------------------------------------------


	}; //-- Fin clase Botlev --//

	//---------------------------------------------------------------------------
	// CHANGE FREQUENCY
	//---------------------------------------------------------------------------

	//-------------------------
	unsigned long BotLev::changeFrequency(int littleQueueSize, int bigQueueSize){

	    //wrapper para los distintos metodos implementados
	    switch(FreqCfg::changeFreq){
	    case 1:
		return changeFrequencyRaw(littleQueueSize, bigQueueSize);
	    case 2:
		return changeFrequencySoft(littleQueueSize, bigQueueSize);
	    case 3:
		return changeFrequencySoft2Levels(littleQueueSize, bigQueueSize);
	    case 4:
		return changeFrequencySoftBig(littleQueueSize, bigQueueSize);
	    case 5:
	    case 7: //OJO: Chapuza. Mirar apagaCore y enciendeCore
		return apagaCoresLittle(littleQueueSize, bigQueueSize);
	    case 6:
	    case 8: //OJO: Chapuza. Mirar apagaCore y enciendeCore
		return apagaCoresBig(littleQueueSize, bigQueueSize);

	    default:
		return 0;
	    }
	}

    
	unsigned long BotLev::changeFrequencyRaw(int littleQueueSize, int bigQueueSize){
	    /* Cambia la frecuencia segÃºn la relaciÃ³n del tamaÃ±o de las colas big/little:
	     *   Doble:  freq[2]
	     *   Triple: freq[3]
	     *    ...
	     */
      	
	    int nxtIdx;
      
      
	    if( littleQueueSize==0 )
		nxtIdx = FreqCfg::maxIdxLittle;
	    else if( bigQueueSize==0 )
		nxtIdx = 0;
	    else
		nxtIdx = (int) (bigQueueSize / littleQueueSize);

	    if( nxtIdx > FreqCfg::maxIdxLittle )
		nxtIdx = FreqCfg::maxIdxLittle;
      

	    return changeLittleFreq(nxtIdx);
	}

	//---------------------------------------------------------------------------
	unsigned long BotLev::changeFrequencySoft(int littleQueueSize, int bigQueueSize){      
	    //Hemos encontrado un pico. Guardmos y ponemos la frecuencia al mÃ¡ximo
	    if(littleQueueSize >= FreqCfg::maxQueueSize){
		FreqCfg::maxQueueSize = littleQueueSize; //Lock??
	
		return changeLittleFreq(0);
	    }

      
	    //Las colas no se encuentran al max.
	    float tamStep = (FreqCfg::maxQueueSize*1.0 / (FreqCfg::maxIdxLittle+1)*1.0);
	    int step = (int) (littleQueueSize / tamStep);
	    step = FreqCfg::maxIdxLittle - step;
      

	    return changeLittleFreq(step);
	}

	//---------------------------------------------------------------------------
	unsigned long BotLev::changeFrequencySoft2Levels(int littleQueueSize, int bigQueueSize){      
	    //Hemos encontrado un pico. Guardmos y ponemos la frecuencia al mÃ¡ximo
	    if(littleQueueSize >= FreqCfg::maxQueueSize){
		FreqCfg::maxQueueSize = littleQueueSize; //Lock??
	
		return changeLittleFreq(0);
	    }

      
	    //Las colas no se encuentran al max.
	    float tamStep = (FreqCfg::maxQueueSize*1.0 / 2.0);

	    if(littleQueueSize >= tamStep) //al maximo
		return changeLittleFreq(0);
	    else
		return changeLittleFreq(FreqCfg::maxIdxLittle);
	  
	}

	//---------------------------------------------------------------------------
	unsigned long BotLev::changeFrequencySoftBig(int littleQueueSize, int bigQueueSize){
	    //Hemos encontrado un pico. Guardmos y ponemos la frecuencia al mÃ¡ximo
	    if(littleQueueSize >= FreqCfg::maxQueueSize){
		FreqCfg::maxQueueSize = littleQueueSize; //Lock??
	
		return changeBigFreq(0);
	    }

      
	    //Las colas no se encuentran al max.
	    float tamStep = (FreqCfg::maxQueueSize*1.0 / (FreqCfg::maxIdxBig+1)*1.0);
	    int step = (int) (littleQueueSize / tamStep);
	    step = FreqCfg::maxIdxBig - step;
      

	    return changeBigFreq(step);
	}


	//---------------------------------------------------------------------------
	unsigned long BotLev::changeLittleFreq(int idxFreq){
	    long long currFreq = cpufreq_get_freq_kernel(CORE_LITTLE);

	    if( FreqCfg::littleFreq[idxFreq] == currFreq)
		return currFreq;

	    FreqCfg::_lock->acquire();
	    unsigned long ret;

	    for(int i=0; i<N_CORE_LITTLE; i++){      
		currFreq = cpufreq_get_freq_kernel(CORE_LITTLE+i);

		if( currFreq > FreqCfg::currLittleIdx){ //disminuimos
		    cpufreq_modify_policy_min(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_modify_policy_max(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
	  
		}
		else if(currFreq < FreqCfg::currLittleIdx){ //aumentamos
		    cpufreq_modify_policy_max(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_modify_policy_min(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_LITTLE+i, FreqCfg::littleFreq[idxFreq]);
		}
	    }

	    FreqCfg::currLittleIdx = idxFreq;
	    ret = cpufreq_get_freq_kernel(CORE_LITTLE);
      
	    FreqCfg::_lock->release();

	    return ret;      
	}

	//---------------------------------------------------------------------------
    
	unsigned long BotLev::changeBigFreq(int idxFreq){
	    long long currFreq = cpufreq_get_freq_kernel(CORE_BIG);

	    if( currFreq == FreqCfg::bigFreq[idxFreq])
		return currFreq;
	  
	    FreqCfg::_lock->acquire();
	    unsigned long ret;


      
	    for(int i=0; i<N_CORE_BIG; i++){
		currFreq = cpufreq_get_freq_kernel(CORE_BIG);

		if(currFreq > FreqCfg::currBigIdx){ //disminuimos
		    cpufreq_modify_policy_min(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_modify_policy_max(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);

		}
		else if(currFreq < FreqCfg::currBigIdx) { //aumentamos
		    cpufreq_modify_policy_max(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_modify_policy_min(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);
		    cpufreq_set_frequency(CORE_BIG+i, FreqCfg::bigFreq[idxFreq]);

		}
	    }


	    FreqCfg::currBigIdx = idxFreq;
	    ret = cpufreq_get_freq_kernel(CORE_BIG);
      
	    FreqCfg::_lock->release();

	    return ret;
	}

	//---------------------------------------------------------------------------
	unsigned long BotLev::apagaCoresLittle(int littleQueueSize, int bigQueueSize){
	    //      return 2000;
      
      
	    //Hemos encontrado un pico. Guardmos y ponemos la frecuencia al maximo
	    if(FreqCfg::maxQueueSize <= 10 || littleQueueSize >= FreqCfg::maxQueueSize ){
	
		FreqCfg::maxQueueSize = MAX(littleQueueSize, FreqCfg::maxQueueSize); //Lock??
		CoresCfg::enciendeLittle();
	
		return 2000;
	    }

      
	    //Las colas no se encuentran al max. 
	    float tamStep = FreqCfg::maxQueueSize * ((CoresCfg::level*1.0f)/ 100.0f); //JUGAR CON ESTE PARAMETRO
      
	    if(littleQueueSize <= tamStep - (0.01*FreqCfg::maxQueueSize)){//por encima del limite
		CoresCfg::apagaLittle();
		return 1000;

	    } else if(littleQueueSize >= tamStep){// + (0.01*FreqCfg::maxQueueSize)) {
		CoresCfg::enciendeLittle();
		return 2000;
	    } else {
		return (CoresCfg::apagadoLittle) ? 1000 : 2000;
	    }
	
	}

	//---------------------------------------------------------------------------
	unsigned long BotLev::apagaCoresBig(int littleQueueSize, int bigQueueSize){
      
	    //Hemos encontrado un pico. Guardmos y ponemos la frecuencia al maximo
	    if(FreqCfg::maxQueueSize <= 10 || littleQueueSize >= FreqCfg::maxQueueSize ){

		FreqCfg::maxQueueSize = MAX(littleQueueSize, FreqCfg::maxQueueSize); //Lock??
		CoresCfg::enciendeBig();
	
		return 2000;
	    }

      
	    //Las colas no se encuentran al max. 
	    float tamStep = FreqCfg::maxQueueSize * ((CoresCfg::level*1.0f)/ 100.0f); //JUGAR CON ESTE PARAMETRO
      
	    if(littleQueueSize <= tamStep - (0.01*FreqCfg::maxQueueSize)){//por encima del limite
		CoresCfg::apagaBig();
		return 1000;

	    } else if(littleQueueSize >= 0.9* tamStep){// + (0.01*FreqCfg::maxQueueSize)) {
		CoresCfg::enciendeBig();
		return 2000;
	    } else {
		return (CoresCfg::apagadoBig) ? 1000 : 2000;
	    }
	
	}
	//---------------------------------------------------------------------------

    
	//---------------------------------------------------------------------------
	// AT IDLE
	//---------------------------------------------------------------------------

	WD * BotLev::atIdle ( BaseThread *thread, int numSteal ) {
	    NANOS_INSTRUMENT( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
		NANOS_INSTRUMENT( static nanos_event_key_t blev_overheads  = ID->getEventKey("blev-overheads"); )
		NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( blev_overheads, 1 ); )
		NANOS_INSTRUMENT( static nanos_event_key_t blev_overheads_br  = ID->getEventKey("blev-overheads-breakdown"); )
		NANOS_INSTRUMENT( sys.getInstrumentation()->raiseOpenBurstEvent ( blev_overheads_br, NANOS_SCHED_BLEV_ATIDLE ); )


		WorkDescriptor * wd;
	    TeamData &data = ( TeamData & ) *thread->getTeam()->getScheduleData();
	    unsigned int spins = BotLevCfg::numSpins;

	    /*@ -------- @*/
	    //#ifdef ODROID
	    if(thread->runningOn()->getId()!=0 && !CoresCfg::puedoEjecutar(thread->runningOn()->getId()) ){
		//#elif JUNO
		//if(!CoresCfg::puedoEjecutar(thread->runningOn()->getId()) ){
		//#else
		//      if( false ){
		//#endif

		return NULL;
	    }
	    //Chapuza para hacer que 0 no ejecute nada tampoco
	    if(thread->runningOn()->getId()==0 && !CoresCfg::puedoEjecutar(1)){
		return NULL;
	    }

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
	    NANOS_INSTRUMENT( sys.getInstrumentation()->raiseCloseBurstEvent ( blev_overheads, 1 ); )
		NANOS_INSTRUMENT( sys.getInstrumentation()->raiseCloseBurstEvent ( blev_overheads_br, NANOS_SCHED_BLEV_ATIDLE ); )

	 
		if(!wd) wd = data._readyQueues[2].pop_front( thread ); 

#ifdef NANOS_INSTRUMENTATION_ENABLED
	    if(wd) {
		NANOS_INSTRUMENT ( static nanos_event_key_t criticalityKey = ID->getEventKey("wd-criticality"); )
		    NANOS_INSTRUMENT ( nanos_event_value_t wd_criticality; )
		    NANOS_INSTRUMENT ( WDData & wddata = *dynamic_cast<WDData*>( wd->getSchedulerData() ); )
		    NANOS_INSTRUMENT ( wd_criticality = wddata.getCriticality(); )
		    NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenBurstEvent ( criticalityKey, wd_criticality ); )	      
		    }

#endif

	    /*@ --------------------------------- @*/
	    if(wd) {
		int b0, b1;
		b0 = data._readyQueues[0].size();
		b1 = data._readyQueues[1].size();
	
		int f =  changeFrequency(b0, b1) / 1000; //small - big


	
		//-- Salida por paraver --//
#ifdef NANOS_INSTRUMENTATION_ENABLED
		NANOS_INSTRUMENT(nanos_event_key_t clave[3];)
		    NANOS_INSTRUMENT(nanos_event_value_t valor[3];)
		    NANOS_INSTRUMENT(clave[0] = 6666; clave[1] = 6667; clave[2] = 6668;)
		    NANOS_INSTRUMENT(valor[0] = b0;   valor[1] = b1;   valor[2] = f;)
	
		    NANOS_INSTRUMENT(sys.getInstrumentation()->raisePointEvents(3, clave, valor); )
#else
		    //-- Salida por fichero de muestras --//
		    data._lockfout->acquire();
		fprintf(data._fout, "data0 = [ data0; %u %d ];\n", (unsigned)time(NULL), b0);
		fprintf(data._fout, "data1 = [ data1; %u %d ];\n", (unsigned)time(NULL), b1);
		fprintf(data._fout, "freq  = [ freq;  %u %d ];\n", (unsigned)time(NULL), f);
		// fflush(data._fout);
		data._lockfout->release();
#endif
	    }
	    /*@ --------------------------------- @*/

	    return wd;
	}



	//---------------------------------------------------------------------------
	// CONFIG
	//---------------------------------------------------------------------------
	class BotLevSchedPlugin : public Plugin
	{
	public:
	    BotLevSchedPlugin() : Plugin( "Distributed Breadth-First scheduling Plugin",1 ) {}

	    virtual void config ( Config &config_ )
		{
		    config_.setOptionsSection( "Bottom level", "Bottom-level scheduling module" );
		    config_.registerConfigOption ( "update-freq", new Config::PositiveVar( BotLevCfg::updateFreq ), "Defines how often to update the bottom levels" );
		    config_.registerArgOption ( "update-freq", "update-freq" );
		    config_.registerEnvOption ( "update-freq", "NX_BL_FREQ" );

		    config_.registerConfigOption ( "numSpins", new Config::PositiveVar( BotLevCfg::numSpins ), "Defines the number of spins in atIdle (ork stealing)" );
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

		    config_.registerConfigOption ( "changeFreq", new Config::IntegerVar( FreqCfg::changeFreq ), "Defines if the frequency should change relative to the queue size" );
		    config_.registerArgOption ( "changeFreq", "changeFreq" );
		    config_.registerEnvOption ( "changeFreq", "NX_CHANGEFREQ" );

		    config_.registerConfigOption( "levelQueue", new Config::IntegerVar( CoresCfg::level),
						  "Porcentaje (sobre 100, entero) del tamanyo de las colas para encender y apagar cores");
		    config_.registerArgOption ( "levelQueue", "levelQueue" );
		    config_.registerEnvOption ( "levelQueue", "NX_LEVELQUEUE" );
	
		}

	    virtual void init() {
		sys.setDefaultSchedulePolicy(new BotLev());
	    }
	};

    }
}

DECLARE_PLUGIN("botlev",nanos::ext::BotLevSchedPlugin);
