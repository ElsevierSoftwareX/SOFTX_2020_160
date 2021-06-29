// -------------------------------------------------------------------------
// Definitions of the search bins used for the GW150914 analysis
// -------------------------------------------------------------------------

// definition of the selection selection cuts 
TCut dqveto("dqveto","veto_cat2_H1 && veto_cat2_L1 && !veto_cat3_H1 && !veto_cat3_L1 && !veto_hveto_H1 && !veto_hveto_L1");
TCut fcut1("fcut1","frequency[0]>32 && frequency[0]<=992");
TCut fcut2("fcut2","frequency[0]>48 && frequency[0]<=992");
TCut netcc("netcc","netcc[0]>0.7");
TCut qveto("qveto","Qveto[0]>0.3 && Qveto[1]>0.3 && Qveto[2]>0.3 && Qveto[3]>0.3");
TCut lveto("lveto","!(bandwidth[0]<5 || (Lveto[1]<5 && Lveto[2]>0.8))");
TCut cbc("cbc","log10(penalty)<0.5 && chirp[1]>1");

// definition of the inclusive bins
TCut unmodeled    = dqveto+fcut1+netcc;
TCut constrained  = dqveto+fcut2+netcc+qveto+lveto;
TCut chirp        = dqveto+fcut2+netcc+qveto+lveto+cbc;

//cout << "unmodeled    : " << unmodeled.GetTitle() << endl;
//cout << "constrained  : " << constrained.GetTitle() << endl;
//cout << "chirp        : " << chirp.GetTitle() << endl;

// definition of the exclusive bins
TCut eunmodeled   = unmodeled+!constrained;
TCut econstrained = constrained+!chirp;

//cout << "eunmodeled   : " << eunmodeled.GetTitle() << endl;
//cout << "econstrained : " << econstrained.GetTitle() << endl;

