# dijet
CMS dijet analysis

Quick starter guide:
root -l -b -q mk_CondFormats.C
root -l -b -q mk_DijetHistosFill.C\(\"ERA\"\)\,\(\"VERSION\"\)  [or 'python runAllIOVs.py']
['hadd' files by hand if using runAllIOVs.py]
root -l -b -q DijetHistosCombine.C+g
root -l -b -q DijetHistosJER.C+g
root -l -b -q DijetHistosOverlay.C+g

Aim of this package is to consolidate multiple analyses done with JetMET and ZeroBias primary data sets together with the corresponding QCD dijet MC. These include: MC truth JEC, jet veto maps, Jet ID, JER pT-dependent scale factors, dijet eta-dependent JEC (L2Res), multijet pT-dependent JEC (L3Res), PF jet composition (JEC global fit), inclusive jet cross section (physics analysis).

The general core design principles are:
- use (JME)NANO files as input (local data storage)
- use TTree::MakeClass as template (beginner-friendly code, re-freshable header)
- use standalone ROOT (quick prototyping)
- store results in hadd'able format, e.g TH2D, TProfile2D (allow partial re-processing)
- utilize TProfile2D over TH3D whenever possible (fast and small)
- combine triggers and eras in post-processing (superfast iteration)
- store each analysis in separate folder (parallel development, switching on/off)
- provide post-processing tools (e.g. JER extraction, unfolding)
- aspire for linear processing (no loops, e.g. L2Res->JER SF->L2Res again)
- provide paper quality plotting (e.g. DP Notes)

Ultimate goal would be to run mature JetMET analyses automatically and directly after Prompt Reco (JEC4Prompt) for each fill of ~1/fb within 24h, to ensure best possible time-stability and DQM monitoring with robust, physics-relevant standard candles (JEC, JER, jet cross section).