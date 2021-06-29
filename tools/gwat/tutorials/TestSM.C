{

  skymap sm1;
  cout << "L1 " << sm1.size() << endl;

  skymap sm2(int(7));
  cout << "L2 " << sm2.size() << endl;

  sm1=sm2;
  cout << "L1 " << sm1.size() << endl;
  exit(0);
}
