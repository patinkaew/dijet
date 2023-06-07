// Purpose: make plots for MPFX introduction talk
// Need to plug in MC truth

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetResolutionObject_cc)
R__LOAD_LIBRARY(JetMETCorrections/Modules/src/JetResolution_cc)

R__LOAD_LIBRARY(MPFXIntro_C)

void mk_MPFIntro() {

  MPFXIntro();
}
