// -------------------------------------------------------------------------
// Definitions of the search bin used for the IMBH analysis
// -------------------------------------------------------------------------

// definition of the selection selection cuts 
TCut dqveto("dqveto","veto_cat2_H1 && veto_cat2_L1 && !veto_cat3_H1 && !veto_cat3_L1 && !veto_hveto_H1 && !veto_hveto_L1");
TCut fcut1("fcut1","frequency[0]>24 && frequency[0]<=256");
TCut netcc_cut("netcc_cut","netcc[0]>0.7");
TCut qveto1("qveto1","!(Qveto[0]<0.3 && Qveto[1]<0.3 && frequency[0] > 60.0)");
TCut qveto2("qveto2","!(Qveto[2]<0.3 && Qveto[3]<0.3 && frequency[0] > 60.0)");
TCut lveto("lveto","!((bandwidth[0]<5 || (Lveto[1]<5 && Lveto[2]>0.8)) && abs(frequency[0]-60.0)<2)");
TCut cbc_cut("cbc_cut","chirp[1]>1");

// definition of the inclusive bins
//TCut imbh    = fcut1+netcc+qveto1+qveto2+lveto+cbc+dqveto;
TCut imbh    = TCut("imbh",(fcut1+netcc_cut+qveto1+qveto2+lveto+cbc_cut+dqveto).GetTitle());
TCut imbh_noqveto    = TCut("imbh_noqveto",(fcut1+netcc_cut+lveto+cbc_cut+dqveto).GetTitle());

//TCut imbh_noqveto    = fcut1+netcc+lveto+cbc+dqveto;
//TCut constrained  = dqveto+fcut2+netcc+qveto+lveto;
//TCut chirp        = dqveto+fcut2+netcc+qveto+lveto+cbc;

//cout << "unmodeled    : " << unmodeled.GetTitle() << endl;
//cout << "constrained  : " << constrained.GetTitle() << endl;
//cout << "chirp        : " << chirp.GetTitle() << endl;

// definition of the exclusive bins
//TCut eunmodeled   = unmodeled+!constrained;
//TCut econstrained = constrained+!chirp;

//cout << "imbh   : " << imbh.GetTitle() << endl;
//cout << "econstrained : " << econstrained.GetTitle() << endl;

