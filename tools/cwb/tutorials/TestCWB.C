{
  //
  // Instantiation of the mdc object using a user_parameter.C config macro file
  // Author : Gabriele Vedovato

  cwb iCWB("config/user_parameters.C");
  int irunID=1;
  iCWB.run(irunID);
}
