#ifndef CWB_PLUGIN_H
#define CWB_PLUGIN_H

// plugin macros, used to share parameters between Plugin and configPlugin

#define CWB_PLUGIN_IMPORT(TYPE,VAR) {                                           \
  TGlobal* _global = (TGlobal*)gROOT->GetGlobal("__p"#VAR,true);                \
  if(_global!=NULL) {                                                           \
    void* gPOINTER=NULL;                                                        \
    memcpy((void*)&gPOINTER,(void*)_global->GetAddress(),sizeof(void*));        \
    VAR = (TYPE)gPOINTER;                                                       \
  } else {                                                                      \
    cout << "CWB_PLUGIN_IMPORT : global variable not found " << #VAR <<  endl;  \
    exit(1);                                                                    \
  }                                                                             \
}

#define CWB_PLUGIN_EXPORT(VAR) {                                                \
  char __cmdline[128];                                                          \
  TGlobal* _global = (TGlobal*)gROOT->GetGlobal("__p"#VAR,true);                \
  if(_global==NULL) sprintf(__cmdline,"void* __p"#VAR" = (void*)%p;",&VAR);     \
  else              sprintf(__cmdline,"__p"#VAR" = (void*)%p;",&VAR);           \
  gROOT->ProcessLine(__cmdline);                                                \
}

#define CWB_PLUGIN_CHECK(TYPE,VAR,CHECK) {                                      \
  TGlobal* _global = (TGlobal*)gROOT->GetGlobal("__p"#VAR,true);                \
  if(_global!=NULL) {                                                           \
    void* gPOINTER=NULL;                                                        \
    memcpy((void*)&gPOINTER,(void*)_global->GetAddress(),sizeof(void*));        \
    VAR = (TYPE)gPOINTER;                                                       \
    CHECK = true;                                                               \
  } else {                                                                      \
    CHECK = false;                                                              \
  }                                                                             \
}

#endif
