#ifndef _FGCS_SCHED_HPP_
#define _FGCS_SCHED_HPP_

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

#ifndef VOLTA1
#error "OJO, NO ES VOLTA1, ERROR"
#else
 #define BIG_FREQ 2000000
 #define LITTLE_FREQ_MAX 1200000
 #define LITTLE_FREQ_MIN 1000000
#endif

// #define NCORES_PER_APP 4
#define COREBIG 1
#define CORELITTLE 0


// #define MAX_BUDGET 11.0
// #define LITTLE_POWER 1.9754
// #define BIG_POWER 5.043375
/*@-- @*/


namespace nanos
{
    namespace ext 
    {
	typedef enum {
	    LITTLE,
	    BIG,
	    FUTURE_BIG,
	    SLEEP,
	    RUNNING,
	    UNDEF
	} E_coreType;
	
	/**-------------------**/
	/** ____BotLevCfg____ **/
	/**-------------------**/
	/* Configuration options for the BotLev scheduling policy */
	struct BotLevCfg
	{
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
	    static int   experiment_type;
	};


	/**---------------------------**/
	/** ____Experiment Config____ **/
	/**---------------------------**/
	struct ExperimentsCfg
	{
	public:
	    static float MaxBudget;
	    static float LittlePower;
	    static float BigPower;
	    static int   nCoresUsed;
	    
	    static int   bigFreq;
	    static int   littleFreq;
	};

	/**----------------------**/
	/** ____BotLevDOData____ **/
	/**----------------------**/
        /* Helper class for the computation of bottom levels 
	   One instance of this class is stored inside each dependableObject.
	*/
	class BotLevDOData : public DOSchedulerData
	{
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


	/**--------------------**/
	/** ____ExpControl____ **/
	/**--------------------**/
	/* Data gathered for controlling the experiments */
	// typedef struct {
	// public:
	//     //Architecture 
	//     Lock _mtxArch;
	//     int _archIdx; // = 0;
	//     int _archMaxIdx; //=2
	    
	//     int _archConf[3][2];
		
	//     // ...

	//     //Thread syncronization
	//     sem_t _sems[NCORES_PER_APP];
	//     int _coreType[NCORES_PER_APP];
	// } exp13_Control;


	
	typedef struct {
	public:
	    //Power budget
	    Lock  mtxBudget;
	    float currentBudget;

	    //dynamic architecture
	    Lock mtxArch;
	    int IwantToPromote;
	    std::vector<int> ImSleeping;
	    NANOS_INSTRUMENT (int nBig);
	    NANOS_INSTRUMENT (int nFutureBig);
	    NANOS_INSTRUMENT (int nLITTLE);
	    NANOS_INSTRUMENT (int nSleeping);
	    
	    
	    E_coreType *coreType;//[NCORES_PER_APP];	    

	    //Thread syncronization
	    sem_t *sems;//[NCORES_PER_APP];



	    
	} exp14_Control;


	// typedef struct {
	// public:
	//     //comunication with ss
	//     Lock ss_com;

	//     E_coreType *coreType; //[NCORES_PER_APP];	    

	//     //Thread syncronization
	//     sem_t *sems; //[NCORES_PER_APP];
	//     dynamic_bitset mask;
	    
	// } exp15_Control;
	
	       
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

	    exp14_Control        _expControl;
	    //exp15_Control        _expControl;
	    
            struct TeamData : public ScheduleTeamData
            {
		// queues of ready tasks to be executed
		WDPriorityQueue<> *_readyQueues;

		Lock* _lockfout;
		FILE* _fout;
		
		TeamData () : ScheduleTeamData() {
		    _readyQueues = NEW WDPriorityQueue<>[3];

		    _lockfout = new Lock();
		    _fout = fopen("MUESTRA.m", "w");
		    fprintf(_fout, "data0=[];\ndata1=[];\nfreq=[];\n");
		    fflush(_fout);
		}
		virtual ~TeamData () {
		    delete[] _readyQueues;

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
            BotLev();
            virtual ~BotLev();

            virtual size_t getTeamDataSize () const;
            virtual size_t getThreadDataSize () const;
            virtual ScheduleTeamData * createTeamData ();
	    virtual ScheduleThreadData * createThreadData ();
	    
            struct WDData : public ScheduleWDData
            {
		int _criticality;

		void setCriticality( int cr )              { _criticality = cr; }
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
	    //exp original
	    inline WD*  exp00_atIdle ( BaseThread *thread, int numSteal );
	    inline void exp00_queue ( BaseThread *thread, WD &wd);
	    // //exp11
	    // inline WD*  exp11_atIdle ( BaseThread *thread, int numSteal );
	    // inline void exp11_queue ( BaseThread *thread, WD &wd );
	    // //exp12
	    // inline WD*  exp12_atIdle ( BaseThread *thread, int numSteal );
	    // inline void exp12_queue ( BaseThread *thread, WD &wd );
	    //exp13
	    // inline void exp13_ctr();
	    // inline void exp13_dtr();
	    // inline WD*  exp13_atIdle ( BaseThread *thread, int numSteal );
	    // inline WD*  exp13_atBeforeExit( BaseThread *thread, WD &wd, bool schedule);
	    // inline void exp13_queue ( BaseThread *thread, WD &wd );
	    //        void exp13_modArch ( int nCrit, int nNotCrit );
	    //exp14
	    inline void exp14_ctr();
	    inline void exp14_dtr();
	    inline WD*  exp14_atIdle ( BaseThread *thread, int numSteal );
	    inline WD*  exp14_atBeforeExit( BaseThread *thread, WD &wd, bool schedule);
	    inline void exp14_queue ( BaseThread *thread, WD &wd );
	           void exp14_modArch ( int coreId, int nCrit, int nNotCrit );
	    inline void exp14_convertToBig(int coreId);
	    inline void exp14_convertToLITTLE(int coreId);
	    inline void exp14_atShutdown();

	    // //exp15
	    // inline void exp15_ctr();
	    // inline void exp15_dtr();
	    // inline WD*  exp15_atIdle ( BaseThread *thread, int numSteal );
	    // inline WD*  exp15_atBeforeExit( BaseThread *thread, WD &wd, bool schedule);
	    // inline void exp15_queue ( BaseThread *thread, WD &wd );
	    // inline void exp15_atShutdown();

	    // inline void exp15_setMask(dynamic_bitset &mask);
	    
	};

	class fgcsSchedPlugin : public Plugin
	{
	public:
            fgcsSchedPlugin();
            virtual void config ( Config &config_ );
            virtual void init();
	};
    }
}


#endif //_FGCS_SCHED_HPP_

