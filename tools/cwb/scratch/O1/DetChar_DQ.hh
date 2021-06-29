// -------------------------------------------------------------------------
// Definitions of the dq sets for DetChar
// -------------------------------------------------------------------------

// definition of the selection selection cuts 
TCut hwinj("hwinj","!veto_cat3_H1 && !veto_cat3_L1");
TCut cat2("cat2","veto_cat2_H1 && veto_cat2_L1");
TCut cat3("cat3","!veto_hveto_H1 && !veto_hveto_L1");
TCut fcut1("fcut1","frequency[0]>32 && frequency[0]<=992");
TCut fcut2("fcut2","frequency[0]>48 && frequency[0]<=992");
TCut netcc("netcc","netcc[0]>0.7");
TCut qveto("qveto","Qveto[0]>0.3 && Qveto[1]>0.3");
TCut lveto("lveto","!(min(bandwidth[0],Lveto[1])<5 && Lveto[2]>0.8)");
TCut cbc("cbc","log10(penalty)<0.5 && chirp[1]>1");

// definition of the dq sets

TCut unmodeled_after_cat1   = hwinj+fcut1+netcc;
TCut unmodeled_after_cat2   = unmodeled_after_cat1+cat2;
TCut unmodeled_after_cat3   = unmodeled_after_cat2+cat3;

TCut constrained_after_cat1 = hwinj+fcut2+netcc+qveto+lveto;
TCut constrained_after_cat2 = constrained_after_cat1+cat2;
TCut constrained_after_cat3 = constrained_after_cat2+cat3;

TCut chirp_after_cat1       = hwinj+fcut2+netcc+qveto+lveto+cbc;
TCut chirp_after_cat2       = chirp_after_cat1+cat2;
TCut chirp_after_cat3       = chirp_after_cat2+cat3;

