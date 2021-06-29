// -------------------------------------------------------------------------
// Definitions of the search bins used for the GW150914 analysis
// -------------------------------------------------------------------------

// definition of the selection selection cuts 
TCut dqveto("dqveto","veto_cat2_H1 && veto_cat2_L1 && !veto_cat3_H1 && !veto_cat3_L1 && !veto_hveto_H1 && !veto_hveto_L1");
TCut fcut1("fcut1","frequency[0]>32 && frequency[0]<=992");
TCut fcut2("fcut2","frequency[0]>48 && frequency[0]<=992");
TCut netcc("netcc","netcc[0]>0.7");
TCut qveto("qveto","Qveto[0]>0.3 && Qveto[1]>0.3");
TCut elveto("elveto","!(bandwidth[0]<5 || (Lveto[1]<5 && Lveto[2]>0.8))");
TCut cbc("cbc","log10(penalty)<0.5 && chirp[1]>1");

// definition of the inclusive bins
TCut xunmodeled    = dqveto+fcut1+netcc;
TCut xconstrained  = dqveto+fcut2+netcc+qveto+elveto;
TCut xchirp        = dqveto+fcut2+netcc+qveto+elveto+cbc;

//cout << "xunmodeled    : " << xunmodeled.GetTitle() << endl;
//cout << "xconstrained  : " << xconstrained.GetTitle() << endl;
//cout << "xchirp        : " << xchirp.GetTitle() << endl;

// definition of the exclusive bins
TCut exunmodeled   = xunmodeled+!xconstrained;
TCut exconstrained = xconstrained+!xchirp;

//cout << "exunmodeled   : " << exunmodeled.GetTitle() << endl;
//cout << "exconstrained : " << exconstrained.GetTitle() << endl;

