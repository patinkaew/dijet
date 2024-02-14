{
  // For JEC
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L testJERSF.C+g");
  testJERSF("pdf/jerSF/Summer20Run2_ZB_v26c_JRV3_MC_SF_AK4PFchs.txt","Run2_ZB_v26c");
  testJERSF("pdf/jerSF/Summer20UL2018_ZB_v26c_JRV3_MC_SF_AK4PFchs.txt","UL2018_ZB_v26c");
  testJERSF("pdf/jerSF/Summer20UL2017_ZB_v26c_JRV3_MC_SF_AK4PFchs.txt","UL2017_ZB_v26c");
  testJERSF("pdf/jerSF/Summer20UL2016GH_ZB_v26c_JRV3_MC_SF_AK4PFchs.txt","UL2016GH_ZB_v26c");
  testJERSF("pdf/jerSF/Summer20UL2016APV_ZB_v26c_JRV3_MC_SF_AK4PFchs.txt","UL2016APV_ZB_v26c");
  //testJERSF("pdf/jerSF/Summer20UL18_JRV3_MC_SF_AK4PFchs.txt");
}
