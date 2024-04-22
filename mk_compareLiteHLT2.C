{

  //gROOT->ProcessLine(".exception");

  gROOT->ProcessLine(".L tools.C+g");

  // Link JEC libraries that are modified to work stand-alone
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L compareLiteHLT2.C+g");
  //gROOT->ProcessLine(".L drawCompareLiteHLT2.C+g");

  compareLiteHLT2("2023D");
  //drawCompareLiteHLT2("2023D");
}
