// @(#)root/physics:$Id: TOracleServer.h 23102 2008-04-09 22:18:41Z pcanal $
// Author: Yan Liu and Shaowen Wang   23/11/04

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TOracleServer
#define ROOT_TOracleServer

#ifndef ROOT_TSQLServer
#include "TSQLServer.h"
#endif

#if !defined(__CINT__)
#ifndef R__WIN32
#include <sys/time.h>
#endif
#include <occi.h>
using namespace std;
using namespace oracle::occi;
#ifdef CONST
#undef CONST
#endif
#else
class Environment;
class Connection;
#endif


class TOracleServer : public TSQLServer {

private:
   Environment  *fEnv;    // environment of Oracle access
   Connection   *fConn;   // connection to Oracle server
   TString       fInfo;  // info string with Oracle version information

public:
   TOracleServer(const char *db, const char *uid, const char *pw);
   ~TOracleServer();

   void        Close(Option_t *opt="");
   TSQLResult *Query(const char *sql);
   Bool_t      Exec(const char* sql);
   TSQLStatement *Statement(const char *sql, Int_t niter = 100);
   Bool_t      IsConnected() const { return (fConn!=0) && (fEnv!=0); }
   Bool_t      HasStatement() const { return kTRUE; }
   Int_t       SelectDataBase(const char *dbname);
   TSQLResult *GetDataBases(const char *wild = 0);
   TSQLResult *GetTables(const char *dbname, const char *wild = 0);
   TList      *GetTablesList(const char* wild = 0);
   TSQLTableInfo *GetTableInfo(const char* tablename);
   TSQLResult *GetColumns(const char *dbname, const char *table, const char *wild = 0);
   Int_t       GetMaxIdentifierLength() { return 30; }
   Int_t       CreateDataBase(const char *dbname);
   Int_t       DropDataBase(const char *dbname);
   Int_t       Reload();
   Int_t       Shutdown();
   const char *ServerInfo();

   Bool_t      StartTransaction();
   Bool_t      Commit();
   Bool_t      Rollback();

   ClassDef(TOracleServer,0)  // Connection to Oracle server
};

#endif
