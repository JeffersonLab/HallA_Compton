{
  const int nMacro = 3;
   TString macros[nMacro] =
    {
      "VJones.C+",
      "MJones.C+",
      "general.C+"
    };

  cout << "Loading macro's:" << endl;
  for(UInt_t i=0; i<nMacro; i++) {
    cout << "\t " << macros[i] << endl;
    gROOT->LoadMacro(macros[i]);
  }
}
