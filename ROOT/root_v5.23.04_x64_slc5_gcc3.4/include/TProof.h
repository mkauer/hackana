// @(#)root/proof:$Id: TProof.h 27933 2009-03-25 15:19:55Z ganis $
// Author: Fons Rademakers   13/02/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProof
#define ROOT_TProof


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProof                                                               //
//                                                                      //
// This class controls a Parallel ROOT Facility, PROOF, cluster.        //
// It fires the worker servers, it keeps track of how many workers are  //
// running, it keeps track of the workers running status, it broadcasts //
// messages to all workers, it collects results, etc.                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TProof
#include "TProof.h"
#endif
#ifndef ROOT_TProofMgr
#include "TProofMgr.h"
#endif
#ifndef ROOT_TProofDebug
#include "TProofDebug.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif
#ifndef ROOT_MessageTypes
#include "MessageTypes.h"
#endif
#ifndef ROOT_TMD5
#include "TMD5.h"
#endif
#ifndef ROOT_TSysEvtHandler
#include "TSysEvtHandler.h"
#endif
#ifndef ROOT_TThread
#include "TThread.h"
#endif
#ifndef ROOT_TUrl
#include "TUrl.h"
#endif

#include <map>

#ifdef R__GLOBALSTL
namespace std { using ::map; }
#endif

#define CANNOTUSE(x) Info(x,"Not manager: cannot use this method")

class TChain;
class TCondor;
class TCondorSlave;
class TDrawFeedback;
class TDSet;
class TEventList;
class THashList;
class TList;
class TCollection;
class TMessage;
class TMonitor;
class TPluginHandler;
class TProof;
class TProofInputHandler;
class TProofInterruptHandler;
class TProofLockPath;
class TVirtualProofPlayer;
class TProofPlayer;
class TProofPlayerRemote;
class TProofProgressDialog;
class TProofServ;
class TQueryResult;
class TSignalHandler;
class TSlave;
class TSemaphore;
class TSocket;
class TTree;
class TVirtualMutex;
class TFileCollection;
class TMap;
class TProofDataSetManager;

// protocol changes:
// 1 -> 2: new arguments for Process() command, option added
// 2 -> 3: package manager enabling protocol changed
// 3 -> 4: introduction of multi-level-master support
// 4 -> 5: added friends support
// 5 -> 6: drop TFTP, support for asynchronous queries
// 6 -> 7: support for multisessions, archieve, retrieve, ...
// 7 -> 8: return number of entries in GetNextPacket
// 8 -> 9: support for stateless connection via xproofd
// 9 -> 10: new features requested, tested at CAF
// 10 -> 11: new merging strategy
// 11 -> 12: new progress message
// 12 -> 13: exchange version/architecture/compiler info
// 13 -> 14: new proofserv environment setting
// 14 -> 15: add support for entry lists; new version of TFileInfo
// 15 -> 16: add support for generic non-data based processing
// 16 -> 17: new dataset handling system; support for TFileCollection processing 
// 17 -> 18: support for reconnection on daemon restarts
// 18 -> 19: TProofProgressStatus used in kPROOF_PROGRESS, kPROOF_STOPPROCESS
//           and kPROOF_GETNEXTPACKET messages in Master - worker communication
// 19 -> 20: Fix the asynchronous mode (required changes in some messages)
// 20 -> 21: Add support for session queuing

// PROOF magic constants
const Int_t       kPROOF_Protocol        = 21;            // protocol version number
const Int_t       kPROOF_Port            = 1093;          // IANA registered PROOF port
const char* const kPROOF_ConfFile        = "proof.conf";  // default config file
const char* const kPROOF_ConfDir         = "/usr/local/root";  // default config dir
const char* const kPROOF_WorkDir         = ".proof";      // default working directory
const char* const kPROOF_CacheDir        = "cache";       // file cache dir, under WorkDir
const char* const kPROOF_PackDir         = "packages";    // package dir, under WorkDir
const char* const kPROOF_QueryDir        = "queries";     // query dir, under WorkDir
const char* const kPROOF_DataSetDir      = "datasets";    // dataset dir, under WorkDir
const char* const kPROOF_CacheLockFile   = "proof-cache-lock-";   // cache lock file
const char* const kPROOF_PackageLockFile = "proof-package-lock-"; // package lock file
const char* const kPROOF_QueryLockFile   = "proof-query-lock-";   // query lock file
const char* const kPROOF_TerminateWorker = "+++ terminating +++"; // signal worker termination in MarkBad
const char* const kPROOF_InputDataFile   = "inputdata.root";      // Default input data file name

#ifndef R__WIN32
const char* const kCP     = "/bin/cp -fp";
const char* const kRM     = "/bin/rm -rf";
const char* const kLS     = "/bin/ls -l";
const char* const kUNTAR  = "%s -c %s/%s | (cd %s; tar xf -)";
const char* const kUNTAR2 = "%s -c %s | (cd %s; tar xf -)";
const char* const kUNTAR3 = "%s -c %s | (tar xf -)";
const char* const kGUNZIP = "gunzip";
#else
const char* const kCP     = "copy";
const char* const kRM     = "delete";
const char* const kLS     = "dir";
const char* const kUNTAR  = "...";
const char* const kUNTAR2 = "...";
const char* const kUNTAR3 = "...";
const char* const kGUNZIP = "gunzip";
#endif

R__EXTERN TVirtualMutex *gProofMutex;

typedef void (*PrintProgress_t)(Long64_t tot, Long64_t proc, Float_t proctime);

// PROOF Interrupt signal handler
class TProofInterruptHandler : public TSignalHandler {
private:

   TProofInterruptHandler(const TProofInterruptHandler&); // Not implemented
   TProofInterruptHandler& operator=(const TProofInterruptHandler&); // Not implemented

   TProof *fProof;
public:
   TProofInterruptHandler(TProof *p)
      : TSignalHandler(kSigInterrupt, kFALSE), fProof(p) { }
   Bool_t Notify();
};

// Input handler for messages from TProofServ
class TProofInputHandler : public TFileHandler {
private:

   TProofInputHandler(const TProofInputHandler&); // Not implemented
   TProofInputHandler& operator=(const TProofInputHandler&); // Not implemented

   TSocket *fSocket;
   TProof  *fProof;
public:
   TProofInputHandler(TProof *p, TSocket *s);
   Bool_t Notify();
   Bool_t ReadNotify() { return Notify(); }
};

// Slaves info class
class TSlaveInfo : public TObject {
public:
   enum ESlaveStatus { kActive, kNotActive, kBad };

   TString      fOrdinal;      //slave ordinal
   TString      fHostName;     //hostname this slave is running on
   TString      fMsd;          //mass storage domain slave is in
   Int_t        fPerfIndex;    //relative performance of this slave
   ESlaveStatus fStatus;       //slave status

   TSlaveInfo(const char *ordinal = "", const char *host = "", Int_t perfidx = 0,
              const char *msd = "") :
              fOrdinal(ordinal), fHostName(host), fMsd(msd),
              fPerfIndex(perfidx), fStatus(kNotActive) { }

   const char *GetMsd() const { return fMsd; }
   const char *GetName() const { return fHostName; }
   const char *GetOrdinal() const { return fOrdinal; }
   void        SetStatus(ESlaveStatus stat) { fStatus = stat; }

   Int_t  Compare(const TObject *obj) const;
   Bool_t IsSortable() const { return kTRUE; }
   void   Print(Option_t *option="") const;

   ClassDef(TSlaveInfo,2) //basic info on slave
};

class TProof : public TNamed, public TQObject {

friend class TPacketizer;
friend class TPacketizerDev;
friend class TPacketizerAdaptive;
friend class TProofLite;
friend class TProofDataSetManager;
friend class TProofServ;
friend class TProofInputHandler;
friend class TProofInterruptHandler;
friend class TProofPlayer;
friend class TProofPlayerLite;
friend class TProofPlayerRemote;
friend class TProofProgressDialog;
friend class TSlave;
friend class TSlaveLite;
friend class TVirtualPacketizer;
friend class TXSlave;
friend class TXSocket;        // to access kPing
friend class TXSocketHandler; // to access fCurrentMonitor and CollectInputFrom
friend class TXProofMgr;      // to access EUrgent
friend class TXProofServ;     // to access EUrgent

public:
   // PROOF status bits
   enum EStatusBits {
      kUsingSessionGui     = BIT(14),
      kNewInputData        = BIT(15),
      kIsClient            = BIT(16),
      kIsMaster            = BIT(17)
   };
   enum EQueryMode {
      kSync                = 0,
      kAsync               = 1
   };
   enum EUploadOpt {
      kAppend              = 0x1,
      kOverwriteDataSet    = 0x2,
      kNoOverwriteDataSet  = 0x4,
      kOverwriteAllFiles   = 0x8,
      kOverwriteNoFiles    = 0x10,
      kAskUser             = 0x0
   };
   enum ERegisterOpt {
      kFailIfExists        = 0,
      kOverwriteIfExists   = 1,
      kMergeIfExists       = 2
   };
   enum EUploadDataSetAnswer {
      kError               = -1,
      kDataSetExists       = -2,
      kFail		   = -3
   };
   enum EUploadPackageOpt {
      kUntar               = 0x0,  //Untar over existing dir [default]
      kRemoveOld           = 0x1   //Remove existing dir with same name
   };
   enum ERunStatus {
      kRunning             = 0,    // Normal status
      kStopped             = 1,    // After the stop button has been pressed
      kAborted             = 2     // After the abort button has been pressed
   };

private:
   enum EUrgent {
      kLocalInterrupt      = -1,
      kPing                = 0,
      kHardInterrupt       = 1,
      kSoftInterrupt,
      kShutdownInterrupt
   };
   enum EProofCacheCommands {
      kShowCache           = 1,
      kClearCache          = 2,
      kShowPackages        = 3,
      kClearPackages       = 4,
      kClearPackage        = 5,
      kBuildPackage        = 6,
      kLoadPackage         = 7,
      kShowEnabledPackages = 8,
      kShowSubCache        = 9,
      kClearSubCache       = 10,
      kShowSubPackages     = 11,
      kDisableSubPackages  = 12,
      kDisableSubPackage   = 13,
      kBuildSubPackage     = 14,
      kUnloadPackage       = 15,
      kDisablePackage      = 16,
      kUnloadPackages      = 17,
      kDisablePackages     = 18,
      kListPackages        = 19,
      kListEnabledPackages = 20,
      kLoadMacro           = 21
   };
   enum EProofDataSetCommands {
      kUploadDataSet       = 1,  //Upload a dataset
      kCheckDataSetName    = 2,  //Check wheter dataset of this name exists
      kGetDataSets         = 3,  //List datasets saved on  the master node
      kRegisterDataSet     = 4,  //Save a TList object as a dataset
      kGetDataSet          = 5,  //Get a TFileCollection of TFileInfo objects
      kVerifyDataSet       = 6,  //Try open all files from a dataset and report results
      kRemoveDataSet       = 7,  //Remove a dataset but leave files belonging to it
      kMergeDataSet        = 8,  //Add new files to an existing dataset
      kShowDataSets        = 9,  //Shows datasets, returns formatted output
      kGetQuota            = 10, //Get quota info per group
      kShowQuota           = 11  //Show quotas
   };
   enum ESendFileOpt {
      kAscii               = 0x0,
      kBinary              = 0x1,
      kForce               = 0x2,
      kForward             = 0x4,
      kCpBin               = 0x8,
      kCp                  = 0x10
   };
   enum EProofWrkListAction {
      kActivateWorker      = 1,
      kDeactivateWorker    = 2
   };
   enum EBuildPackageOpt {
      kDontBuildOnClient   = -2,
      kBuildOnSlavesNoWait = -1,
      kBuildAll            = 0,
      kCollectBuildResults = 1
   };
   enum EProofShowQuotaOpt {
      kPerGroup = 0x1,
      kPerUser = 0x2
   };

   Bool_t          fValid;           //is this a valid proof object
   TString         fMaster;          //master server ("" if a master); used in the browser
   TString         fWorkDir;         //current work directory on remote servers
   Int_t           fLogLevel;        //server debug logging level
   Int_t           fStatus;          //remote return status (part of kPROOF_LOGDONE)
   Int_t           fCheckFileStatus; //remote return status after kPROOF_CHECKFILE
   TList          *fRecvMessages;    //Messages received during collect not yet processed
   TList          *fSlaveInfo;       //!list returned by kPROOF_GETSLAVEINFO
   Bool_t          fSendGroupView;   //if true send new group view
   TList          *fActiveSlaves;    //list of active slaves (subset of all slaves)
   TList          *fInactiveSlaves;  //list of inactive slaves (good but not used for processing)
   TList          *fUniqueSlaves;    //list of all active slaves with unique file systems
   TList          *fAllUniqueSlaves;  //list of all active slaves with unique file systems, including all submasters
   TList          *fNonUniqueMasters; //list of all active masters with a nonunique file system
   TMonitor       *fActiveMonitor;   //monitor activity on all active slave sockets
   TMonitor       *fUniqueMonitor;   //monitor activity on all unique slave sockets
   TMonitor       *fAllUniqueMonitor; //monitor activity on all unique slave sockets, including all submasters
   TMonitor       *fCurrentMonitor;  //currently active monitor
   Long64_t        fBytesRead;       //bytes read by all slaves during the session
   Float_t         fRealTime;        //realtime spent by all slaves during the session
   Float_t         fCpuTime;         //CPU time spent by all slaves during the session
   TSignalHandler *fIntHandler;      //interrupt signal handler (ctrl-c)
   TPluginHandler *fProgressDialog;  //progress dialog plugin
   Bool_t          fProgressDialogStarted; //indicates if the progress dialog is up
   TVirtualProofPlayer *fPlayer;     //current player
   TList          *fFeedback;        //list of names to be returned as feedback
   TList          *fChains;          //chains with this proof set
   struct MD5Mod_t {
      TMD5   fMD5;                   //file's md5
      Long_t fModtime;               //file's modification time
   };
   typedef std::map<TString, MD5Mod_t> FileMap_t;
   FileMap_t       fFileMap;         //map keeping track of a file's md5 and mod time
   TDSet          *fDSet;            //current TDSet being validated

   Int_t           fNotIdle;         //Number of non-idle sub-nodes
   Bool_t          fSync;            //true if type of currently processed query is sync
   ERunStatus      fRunStatus;       //run status
   Bool_t          fIsWaiting;       //true if queries have been enqueued

   Bool_t          fRedirLog;        //redirect received log info
   TString         fLogFileName;     //name of the temp file for redirected logs
   FILE           *fLogFileW;        //temp file to redirect logs
   FILE           *fLogFileR;        //temp file to read redirected logs
   Bool_t          fLogToWindowOnly; //send log to window only

   TList          *fWaitingSlaves;   //stores a TPair of the slaves's TSocket and TMessage
   TList          *fQueries;         //list of TProofQuery objects
   Int_t           fOtherQueries;    //number of queries in list from previous sessions
   Int_t           fDrawQueries;     //number of draw queries during this sessions
   Int_t           fMaxDrawQueries;  //max number of draw queries kept
   Int_t           fSeqNum;          //Remote sequential # of the last query submitted

   Int_t           fSessionID;       //remote ID of the session

   Bool_t          fEndMaster;       //true for a master in direct contact only with workers

   TString         fPackageDir;      //package directory (used on client)
   THashList      *fGlobalPackageDirList;//list of directories containing global packages libs
   TProofLockPath *fPackageLock;     //package lock
   TList          *fEnabledPackagesOnClient; //list of packages enabled on client

   TList          *fInputData;       //Input data objects sent over via file
   TString         fInputDataFile;   //File with input data objects

   PrintProgress_t fPrintProgress;   //Function function to display progress info in batch mode

   TVirtualMutex  *fCloseMutex;      // Avoid crashes in MarkBad or alike while closing

   TList          *fLoadedMacros;    // List of loaded macros (just file names)
   static TList   *fgProofEnvList;   // List of TNameds defining environment
                                     // variables to pass to proofserv
protected:
   enum ESlaves { kAll, kActive, kUnique, kAllUnique };

   Bool_t          fMasterServ;     //true if we are a master server
   TUrl            fUrl;            //Url of the master
   TString         fConfFile;       //file containing config information
   TString         fConfDir;        //directory containing cluster config information
   TString         fImage;          //master's image name
   Int_t           fProtocol;       //remote PROOF server protocol version number
   TList          *fSlaves;         //list of all slave servers as in config file
   TList          *fBadSlaves;      //dead slaves (subset of all slaves)
   TMonitor       *fAllMonitor;     //monitor activity on all valid slave sockets
   Bool_t          fDataReady;      //true if data is ready to be analyzed
   Long64_t        fBytesReady;     //number of bytes staged
   Long64_t        fTotalBytes;     //number of bytes to be analyzed
   TList          *fAvailablePackages; //list of available packages
   TList          *fEnabledPackages;   //list of enabled packages
   TList          *fRunningDSets;   // Temporary datasets used for async running

   Int_t           fCollectTimeout; // Timeout for (some) collect actions

   TString         fDataPoolUrl;    // default data pool entry point URL
   TProofMgr::EServType fServType;  // type of server: proofd, XrdProofd
   TProofMgr      *fManager;        // manager to which this session belongs (if any)
   EQueryMode      fQueryMode;      // default query mode
   Bool_t          fDynamicStartup; // are the workers started dynamically?

   static TSemaphore *fgSemaphore;   //semaphore to control no of parallel startup threads

private:
   TProof(const TProof &);           // not implemented
   void operator=(const TProof &);   // idem

   void     CleanGDirectory(TList *ol);

   Int_t    Exec(const char *cmd, ESlaves list, Bool_t plusMaster);
   Int_t    SendCommand(const char *cmd, ESlaves list = kActive);
   Int_t    SendCurrentState(ESlaves list = kActive);
   Bool_t   CheckFile(const char *file, TSlave *sl, Long_t modtime, Int_t cpopt = (kCp | kCpBin));
   Int_t    SendFile(const char *file, Int_t opt = (kBinary | kForward | kCp | kCpBin),
                     const char *rfile = 0, TSlave *sl = 0);
   Int_t    SendObject(const TObject *obj, ESlaves list = kActive);
   Int_t    SendGroupView();
   Int_t    SendInitialState();
   Int_t    SendPrint(Option_t *option="");
   Int_t    Ping(ESlaves list);
   void     Interrupt(EUrgent type, ESlaves list = kActive);
   void     AskStatistics();
   void     AskParallel();
   Int_t    GoParallel(Int_t nodes, Bool_t accept = kFALSE, Bool_t random = kFALSE);
   Int_t    SetParallelSilent(Int_t nodes, Bool_t random = kFALSE);
   void     RecvLogFile(TSocket *s, Int_t size);
   void     NotifyLogMsg(const char *msg, const char *sfx = "\n");
   Int_t    BuildPackage(const char *package, EBuildPackageOpt opt = kBuildAll);
   Int_t    BuildPackageOnClient(const TString &package);
   Int_t    LoadPackage(const char *package, Bool_t notOnClient = kFALSE);
   Int_t    LoadPackageOnClient(const TString &package);
   Int_t    UnloadPackage(const char *package);
   Int_t    UnloadPackageOnClient(const char *package);
   Int_t    UnloadPackages();
   Int_t    UploadPackageOnClient(const TString &package, EUploadPackageOpt opt, TMD5 *md5);
   Int_t    DisablePackage(const char *package);
   Int_t    DisablePackageOnClient(const char *package);
   Int_t    DisablePackages();

   void     Activate(TList *slaves = 0);
   Int_t    Broadcast(const TMessage &mess, TList *slaves);
   Int_t    Broadcast(const TMessage &mess, ESlaves list = kActive);
   Int_t    Broadcast(const char *mess, Int_t kind, TList *slaves);
   Int_t    Broadcast(const char *mess, Int_t kind = kMESS_STRING, ESlaves list = kActive);
   Int_t    Broadcast(Int_t kind, TList *slaves) { return Broadcast(0, kind, slaves); }
   Int_t    Broadcast(Int_t kind, ESlaves list = kActive) { return Broadcast(0, kind, list); }
   Int_t    BroadcastFile(const char *file, Int_t opt, const char *rfile, TList *wrks);
   Int_t    BroadcastFile(const char *file, Int_t opt, const char *rfile = 0, ESlaves list = kAllUnique);
   Int_t    BroadcastGroupPriority(const char *grp, Int_t priority, ESlaves list = kAllUnique);
   Int_t    BroadcastGroupPriority(const char *grp, Int_t priority, TList *workers);
   Int_t    BroadcastObject(const TObject *obj, Int_t kind, TList *slaves);
   Int_t    BroadcastObject(const TObject *obj, Int_t kind = kMESS_OBJECT, ESlaves list = kActive);
   Int_t    BroadcastRaw(const void *buffer, Int_t length, TList *slaves);
   Int_t    BroadcastRaw(const void *buffer, Int_t length, ESlaves list = kActive);
   Int_t    Collect(const TSlave *sl, Long_t timeout = -1, Int_t endtype = -1);
   Int_t    Collect(TMonitor *mon, Long_t timeout = -1, Int_t endtype = -1);
   Int_t    CollectInputFrom(TSocket *s, Int_t endtype = -1);
   Int_t    HandleInputMessage(TSlave *wrk, TMessage *m);
   void     SetMonitor(TMonitor *mon = 0, Bool_t on = kTRUE);

   void     ReleaseMonitor(TMonitor *mon);

   void     FindUniqueSlaves();
   TSlave  *FindSlave(TSocket *s) const;
   TList   *GetListOfSlaves() const { return fSlaves; }
   TList   *GetListOfInactiveSlaves() const { return fInactiveSlaves; }
   TList   *GetListOfUniqueSlaves() const { return fUniqueSlaves; }
   TList   *GetListOfBadSlaves() const { return fBadSlaves; }
   Int_t    GetNumberOfSlaves() const;
   Int_t    GetNumberOfActiveSlaves() const;
   Int_t    GetNumberOfInactiveSlaves() const;
   Int_t    GetNumberOfUniqueSlaves() const;
   Int_t    GetNumberOfBadSlaves() const;

   Bool_t   IsEndMaster() const { return fEndMaster; }
   void     ModifyWorkerLists(const char *ord, Bool_t add);

   Bool_t   IsSync() const { return fSync; }
   void     InterruptCurrentMonitor();

   void     SetRunStatus(ERunStatus rst) { fRunStatus = rst; }

   void     MarkBad(TSlave *wrk, const char *reason = 0);
   void     MarkBad(TSocket *s, const char *reason = 0);
   void     TerminateWorker(TSlave *wrk);
   void     TerminateWorker(const char *ord);

   void     ActivateAsyncInput();
   void     DeActivateAsyncInput();

   Int_t    GetQueryReference(Int_t qry, TString &ref);

   void     PrintProgress(Long64_t total, Long64_t processed, Float_t procTime = -1.);

protected:
   TProof(); // For derived classes to use
   Int_t           Init(const char *masterurl, const char *conffile,
                        const char *confdir, Int_t loglevel,
                        const char *alias = 0);
   virtual Bool_t  StartSlaves(Bool_t attach = kFALSE);
   Int_t AddWorkers(TList *wrks);
   Int_t RemoveWorkers(TList *wrks);

   void                         SetPlayer(TVirtualProofPlayer *player);
   TVirtualProofPlayer         *GetPlayer() const { return fPlayer; }
   virtual TVirtualProofPlayer *MakePlayer(const char *player = 0, TSocket *s = 0);

   void    UpdateDialog();

   void    HandleLibIncPath(const char *what, Bool_t add, const char *dirs);

   TList  *GetListOfActiveSlaves() const { return fActiveSlaves; }
   TSlave *CreateSlave(const char *url, const char *ord,
                       Int_t perf, const char *image, const char *workdir);
   TSlave *CreateSubmaster(const char *url, const char *ord,
                           const char *image, const char *msd);

   virtual void SaveWorkerInfo();

   Int_t    Collect(ESlaves list = kActive, Long_t timeout = -1, Int_t endtype = -1);
   Int_t    Collect(TList *slaves, Long_t timeout = -1, Int_t endtype = -1);

   void         SetDSet(TDSet *dset) { fDSet = dset; }
   virtual void ValidateDSet(TDSet *dset);

   TPluginHandler *GetProgressDialog() const { return fProgressDialog; }

   Int_t AssertPath(const char *path, Bool_t writable);

   void PrepareInputDataFile(TString &dataFile);
   virtual void SendInputDataFile();

   static void *SlaveStartupThread(void *arg);

   static Int_t AssertDataSet(TDSet *dset, TList *input,
                              TProofDataSetManager *mgr, TString &emsg);
   // Input data handling
   static Int_t GetInputData(TList *input, const char *cachedir, TString &emsg);
   static Int_t SaveInputData(TQueryResult *qr, const char *cachedir, TString &emsg);
   static Int_t SendInputData(TQueryResult *qr, TProof *p, TString &emsg);

   // Parse CINT commands
   static Bool_t GetFileInCmd(const char *cmd, TString &fn);

public:
   TProof(const char *masterurl, const char *conffile = kPROOF_ConfFile,
          const char *confdir = kPROOF_ConfDir, Int_t loglevel = 0,
          const char *alias = 0, TProofMgr *mgr = 0);
   virtual ~TProof();

   void        cd(Int_t id = -1);

   Int_t       Ping();
   void        Touch();
   Int_t       Exec(const char *cmd, Bool_t plusMaster = kFALSE);

   virtual Long64_t Process(TDSet *dset, const char *selector,
                            Option_t *option = "", Long64_t nentries = -1,
                            Long64_t firstentry = 0);
   virtual Long64_t Process(TFileCollection *fc, const char *selector,
                            Option_t *option = "", Long64_t nentries = -1,
                            Long64_t firstentry = 0);
   virtual Long64_t Process(const char *dsetname, const char *selector,
                            Option_t *option = "", Long64_t nentries = -1,
                            Long64_t firstentry = 0, TObject *enl = 0);
   virtual Long64_t Process(const char *selector, Long64_t nentries,
                            Option_t *option = "");

   virtual Long64_t DrawSelect(TDSet *dset, const char *varexp,
                          const char *selection = "",
                          Option_t *option = "", Long64_t nentries = -1,
                          Long64_t firstentry = 0);
   Long64_t    DrawSelect(const char *dsetname, const char *varexp,
                          const char *selection = "",
                          Option_t *option = "", Long64_t nentries = -1,
                          Long64_t firstentry = 0, TObject *enl = 0);
   Int_t       Archive(Int_t query, const char *url);
   Int_t       Archive(const char *queryref, const char *url = 0);
   Int_t       CleanupSession(const char *sessiontag);
   Long64_t    Finalize(Int_t query = -1, Bool_t force = kFALSE);
   Long64_t    Finalize(const char *queryref, Bool_t force = kFALSE);
   Int_t       Remove(Int_t query, Bool_t all = kFALSE);
   Int_t       Remove(const char *queryref, Bool_t all = kFALSE);
   Int_t       Retrieve(Int_t query, const char *path = 0);
   Int_t       Retrieve(const char *queryref, const char *path = 0);

   void        StopProcess(Bool_t abort, Int_t timeout = -1);
   void        Browse(TBrowser *b);

   Int_t       SetParallel(Int_t nodes = 9999, Bool_t random = kFALSE);
   void        SetLogLevel(Int_t level, UInt_t mask = TProofDebug::kAll);

   void        Close(Option_t *option="");
   virtual void Print(Option_t *option="") const;

   //-- cache and package management
   virtual void ShowCache(Bool_t all = kFALSE);
   virtual void ClearCache(const char *file = 0);
   TList      *GetListOfPackages();
   TList      *GetListOfEnabledPackages();
   void        ShowPackages(Bool_t all = kFALSE);
   void        ShowEnabledPackages(Bool_t all = kFALSE);
   Int_t       ClearPackages();
   Int_t       ClearPackage(const char *package);
   Int_t       EnablePackage(const char *package, Bool_t notOnClient = kFALSE);
   Int_t       UploadPackage(const char *par, EUploadPackageOpt opt = kUntar);
   Int_t       Load(const char *macro, Bool_t notOnClient = kFALSE, Bool_t uniqueOnly = kTRUE,
                    TList *wrks = 0);

   Int_t       AddDynamicPath(const char *libpath, Bool_t onClient = kFALSE, TList *wrks = 0);
   Int_t       AddIncludePath(const char *incpath, Bool_t onClient = kFALSE, TList *wrks = 0);
   Int_t       RemoveDynamicPath(const char *libpath, Bool_t onClient = kFALSE);
   Int_t       RemoveIncludePath(const char *incpath, Bool_t onClient = kFALSE);

   //-- dataset management
   Int_t       UploadDataSet(const char *dataset,
                             TList *files,
                             const char *dest = 0,
                             Int_t opt = kAskUser,
                             TList *skippedFiles = 0);
   Int_t       UploadDataSet(const char *dataset,
                             const char *files,
                             const char *dest = 0,
                             Int_t opt = kAskUser,
                             TList *skippedFiles = 0);
   Int_t       UploadDataSetFromFile(const char *dataset,
                                     const char *file,
                                     const char *dest = 0,
                                     Int_t opt = kAskUser,
                                     TList *skippedFiles = 0);
   virtual Bool_t  RegisterDataSet(const char *name,
                               TFileCollection *dataset, const char* optStr = "");
   virtual TMap *GetDataSets(const char *uri = "", const char* optStr = "");
   virtual void  ShowDataSets(const char *uri = "", const char* optStr = "");

   TMap       *GetDataSetQuota(const char* optStr = "");
   void        ShowDataSetQuota(Option_t* opt = 0);

   void        ShowDataSet(const char *dataset = "", const char* opt = "M");
   virtual Int_t RemoveDataSet(const char *dataset, const char* optStr = "");
   virtual Int_t VerifyDataSet(const char *dataset, const char* optStr = "");
   virtual TFileCollection *GetDataSet(const char *dataset, const char* optStr = "");
   TList       *FindDataSets(const char *searchString, const char* optStr = "");

   const char *GetMaster() const { return fMaster; }
   const char *GetConfDir() const { return fConfDir; }
   const char *GetConfFile() const { return fConfFile; }
   const char *GetUser() const { return fUrl.GetUser(); }
   const char *GetWorkDir() const { return fWorkDir; }
   const char *GetSessionTag() const { return GetName(); }
   const char *GetImage() const { return fImage; }
   const char *GetUrl() { return fUrl.GetUrl(); }
   Int_t       GetPort() const { return fUrl.GetPort(); }
   Int_t       GetRemoteProtocol() const { return fProtocol; }
   Int_t       GetClientProtocol() const { return kPROOF_Protocol; }
   Int_t       GetStatus() const { return fStatus; }
   Int_t       GetLogLevel() const { return fLogLevel; }
   Int_t       GetParallel() const;
   Int_t       GetSessionID() const { return fSessionID; }
   TList      *GetListOfSlaveInfos();
   Bool_t      UseDynamicStartup() const { return fDynamicStartup; }

   EQueryMode  GetQueryMode(Option_t *mode = 0) const;
   void        SetQueryMode(EQueryMode mode);

   void        SetRealTimeLog(Bool_t on = kTRUE);

   Long64_t    GetBytesRead() const { return fBytesRead; }
   Float_t     GetRealTime() const { return fRealTime; }
   Float_t     GetCpuTime() const { return fCpuTime; }

   Bool_t      IsLite() const { return (fServType == TProofMgr::kProofLite); }
   Bool_t      IsProofd() const { return (fServType == TProofMgr::kProofd); }
   Bool_t      IsFolder() const { return kTRUE; }
   Bool_t      IsMaster() const { return fMasterServ; }
   Bool_t      IsValid() const { return fValid; }
   Bool_t      IsParallel() const { return GetParallel() > 0 ? kTRUE : kFALSE; }
   Bool_t      IsIdle() const { return (fNotIdle <= 0) ? kTRUE : kFALSE; }
   Bool_t      IsWaiting() const { return fIsWaiting; }

   ERunStatus  GetRunStatus() const { return fRunStatus; }
   TList      *GetLoadedMacros() const { return fLoadedMacros; }

   //-- input list parameter handling
   void        SetParameter(const char *par, const char *value);
   void        SetParameter(const char *par, Int_t value);
   void        SetParameter(const char *par, Long_t value);
   void        SetParameter(const char *par, Long64_t value);
   void        SetParameter(const char *par, Double_t value);
   TObject    *GetParameter(const char *par) const;
   void        DeleteParameters(const char *wildcard);
   void        ShowParameters(const char *wildcard = "PROOF_*") const;

   void        AddInput(TObject *obj);
   void        ClearInput();
   TList      *GetInputList();
   TObject    *GetOutput(const char *name);
   TList      *GetOutputList();

   void        AddInputData(TObject *obj, Bool_t push = kFALSE);
   void        SetInputDataFile(const char *datafile);
   void        ClearInputData(TObject *obj = 0);
   void        ClearInputData(const char *name);

   void        AddFeedback(const char *name);
   void        RemoveFeedback(const char *name);
   void        ClearFeedback();
   void        ShowFeedback() const;
   TList      *GetFeedbackList() const;

   virtual TList *GetListOfQueries(Option_t *opt = "");
   Int_t       GetNumberOfQueries();
   Int_t       GetNumberOfDrawQueries() { return fDrawQueries; }
   TList      *GetQueryResults();
   TQueryResult *GetQueryResult(const char *ref = 0);
   void        GetMaxQueries();
   void        SetMaxDrawQueries(Int_t max);
   void        ShowQueries(Option_t *opt = "");

   Bool_t      IsDataReady(Long64_t &totalbytes, Long64_t &bytesready);

   void        SetActive(Bool_t /*active*/ = kTRUE) { }

   void        LogMessage(const char *msg, Bool_t all); //*SIGNAL*
   void        Progress(Long64_t total, Long64_t processed); //*SIGNAL*
   void        Progress(Long64_t total, Long64_t processed, Long64_t bytesread,
                        Float_t initTime, Float_t procTime,
                        Float_t evtrti, Float_t mbrti); // *SIGNAL*
   void        Feedback(TList *objs); //*SIGNAL*
   void        QueryResultReady(const char *ref); //*SIGNAL*
   void        CloseProgressDialog(); //*SIGNAL*
   void        ResetProgressDialog(const char *sel, Int_t sz,
                                   Long64_t fst, Long64_t ent); //*SIGNAL*
   void        StartupMessage(const char *msg, Bool_t status, Int_t done,
                              Int_t total); //*SIGNAL*
   void        DataSetStatus(const char *msg, Bool_t status,
                             Int_t done, Int_t total); //*SIGNAL*

   void        SendDataSetStatus(const char *msg, UInt_t n, UInt_t tot, Bool_t st);

   void        GetLog(Int_t start = -1, Int_t end = -1);
   void        PutLog(TQueryResult *qr);
   void        ShowLog(Int_t qry = -1);
   void        ShowLog(const char *queryref);
   Bool_t      SendingLogToWindow() const { return fLogToWindowOnly; }
   void        SendLogToWindow(Bool_t mode) { fLogToWindowOnly = mode; }

   void        ResetProgressDialogStatus() { fProgressDialogStarted = kFALSE; }

   virtual TTree *GetTreeHeader(TDSet *tdset);
   TList      *GetOutputNames();

   void        AddChain(TChain *chain);
   void        RemoveChain(TChain *chain);

   TDrawFeedback *CreateDrawFeedback();
   void           SetDrawFeedbackOption(TDrawFeedback *f, Option_t *opt);
   void           DeleteDrawFeedback(TDrawFeedback *f);

   void        Detach(Option_t *opt = "");

   virtual void SetAlias(const char *alias="");

   TProofMgr  *GetManager() { return fManager; }
   void        SetManager(TProofMgr *mgr);

   void        ActivateWorker(const char *ord);
   void        DeactivateWorker(const char *ord);

   const char *GetDataPoolUrl() const { return fDataPoolUrl; }
   void        SetDataPoolUrl(const char *url) { fDataPoolUrl = url; }

   void        SetPrintProgress(PrintProgress_t pp) { fPrintProgress = pp; }

   // Opening and managing PROOF connections
   static TProof       *Open(const char *url = 0, const char *conffile = 0,
                             const char *confdir = 0, Int_t loglevel = 0);
   static TProofMgr    *Mgr(const char *url);
   static void          Reset(const char *url, Bool_t hard = kFALSE);

   static void          AddEnvVar(const char *name, const char *value);
   static void          DelEnvVar(const char *name);
   static const TList  *GetEnvVars();
   static void          ResetEnvVars();

   // Input/output list utilities
   static Int_t         GetParameter(TCollection *c, const char *par, TString &value);
   static Int_t         GetParameter(TCollection *c, const char *par, Int_t &value);
   static Int_t         GetParameter(TCollection *c, const char *par, Long_t &value);
   static Int_t         GetParameter(TCollection *c, const char *par, Long64_t &value);
   static Int_t         GetParameter(TCollection *c, const char *par, Double_t &value);

   ClassDef(TProof,0)  //PROOF control class
};

// Global object with default PROOF session
R__EXTERN TProof *gProof;

#endif
