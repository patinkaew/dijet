How to RUN on Hefaistos:
------------------------
(try mosh to not drop connection)
- from local: 'rsync -rutP DijetHistosFill.C DijetHistosFill.h Hefaistos:/media/storage/dijet/'
- source /work/data/rootbinaries/root/bin/thisroot.sh [6.26.10]
  [was: source /work/data/root/bin/thisroot.sh]
- (rm *.d *.so *.pcm)
- root -l -b -q mk_CondFormats.C
- #define GPU in mk_DijetHistosFill.C
=> edit (version, IOV_list) and execute 'python runAllIOVs.py'
[- nohup root -l -b -q mk_DijetHistosFill.C\(\"X\"\) > log.txt & ]
+ runtime about X
=> edit (version, IOV_list) and execute 'python renameAllIOVs.py'

+ tail -f log.txt
+ starting up takes quite a bit (~X sec) due to GetEntries() call
   => code TStopWatch to ignore this startup time, or skip GetEntries

- rsync -rutP files from Hefaistos

To-do:
- add proper MET for MPF => done
- add Crecoil for multijets +more advanced 2D version => done 1D+2D
- add PF composition folder
- add trigger turn-on folder
- add jet veto maps
- add MET filters => done: Flag_METFilter; also need others?
- multijet: improve high pT efficiency (relative jet veto around leading jet?)
  - also figure out O(5-10%) mpfu even at highest pT
- add smeared collection of jets for MC
- option: add monitoring of <rho> for each analysis?

// v18: Implement pT,avp and proper bisector for dijet and multijet following JME-21-001. Keep previous pT,ave and bisector as dijet and multijet axes, respectively. Raise pT of recoil jets to 30 GeV in multijet in hopes of reducing MPFu. Add PF composition for incjet, dijet, multijet. Add 15-40 GeV bins for dijet.

(Implement MC (flat, pThat binned, HT binned))
// v17: Add Jet_jetId>=4 (tightLepVeto) and Flag_METFilters>0 (all?) for multijet, dijet, inclusive jet and jet veto maps

// v16: Fix linking of fixedGridRhoFastjetAll (Rho_*) for Run 2. Modify multijet selection (drop all |eta|>2.5 jets from recoil, not just 1+2), patch pT2/pT1<0.6 to pT2/pTrecoil<0.6 (5% bias!), veto all forward jets, add Crecoil 1D+2D. Implement proper L1L2L3-RC MET for MPF. Clean file structure.

// v15: Add UL2016GH data+MC. Run time about 9h for data (77M/245M=>7.5k)? Added JetA and Rho as new inputs. Multijet still not working.

// v14: fix absetamin thresholds for mt["HLT_PFJetFwdX"] for X>=140. Add missing HLT_PFJet400 trigger. Add Jetveto/heta_pretrg for minitools/drawFwdTrigEff.C.

// v13:	invert reference axes on Multijet so that result is pTlead/pTrecoil (v12 and earlier effectively pTrecoil/pTlead). Add requirement that jet[1] and jet[2] at |eta|<2.5, plus check both for pT>30 and pT<0.7*pT[0]
// Implement also new dijetHistos2 for efficient analysis.
// RunE v12e: 5.5/fb (5.476540005), RunF v12f: 1.4/fb (1.363428299)
// RunCD: 7.6/fb (7.561991506)

// v12: Add p2chf,p2nef,p2nhf and p2chftp,p2chftp,p2nhftp to Jetveto. Update JSON file for 2022BCDEF. Add RunF to mk_*.C listing. Add multijet/h2m2a control.
=> todo: add Dijet folder with TProfile2D's vs eta and pT

// v11: Add doJetveto and Jetveto folder. Changed multijet alpha veto from
//      DELTAPHI to DELTAR (avoid bias on forward PU jets)

// v10: Add HLT_ZeroBias. Split text files to ZeroBias and JetMET.
//      Replace TrigObj with TrigObjJMEAK4. Add TH2D h2m0a and
//      hcosdphi for multijet controls
// Oct 12 full JSON: 9.705969501/fb. RunE only: 2.058990167/fb. (ABCD 7.6/fb)
// (JMENANO on eos: v2p0 4.9 TB, v2p1 6.1 TB)

// v9: Add inclusive jet and multijet folders. Add trigger counts.
// v9mc: Switch off L2L3Res

// v8: add UL2018 samples. Implement redoing JEC, add HLT JES with T&P
// data: 16:44->20:42 (3h58min). TrigObj doubling time?
// mc: ~20:42->01:15 (4h33min). TrigObj doubling time?

// v7: add PF composition vs eta and pT. Add hmjj2 and hmjj213
// data 17:15->19:20 (2h5min). mc 21:15->23:27 (2h12min, died).
// mcflat: 08:04->10:00 (1h56)
// ul18a: 10:24->16:46 (6h22)
// Found 10,793,303 bad events according to new JSON (events cut)
// Found 113,084,154 bad events according to trigger bits (events cut)
// Processed 117 runs, 55,380 luminosity blocks and 160,691,332 events
// Analyzed 47,607,178 events
// ul18mc: 9:14->11:59 (2h45h)
// Processed 1 runs, 19,928 luminosity blocks and 19,928,000 events
// Analyzed 19,928,000 events

// v6bothmc on dataFile_QCDFlats.txt
Loaded 38,809,694 entries
Processed 1 runs, 19,996 luminosity blocks and 38,809,694 events
(19:41->23:56->4h15min)

// v6flatmc on dataFile_FlatQCD.txt (new version)
Loaded 18,904,000 entries
Processed 1 runs, 18,904 luminosity blocks and 18,904,000 events
(15:48->17:47, 1h59min)

// v6: run on dataFiles_RunC.txt with 4.86/fb golden JSON => 2.8954/fb
Loaded 129,672,779 entries
Found 49,796,825 bad events according to new JSON (events cut)
Found 67,883,399 bad events according to trigger bits (events cut)
Processed 82 runs, 16,935 luminosity blocks and 79,875,954 events
Saving these to file rootfiles/jmenano_v6.json for brilcalc
Analyzed 11,992,555 events
(21:18->23:13, 1h 55min)

// v5: 4.86/fb golden JSON
Found 4,833,924 bad events according to new JSON (events cut)
Found 8,891,715 bad events according to trigger bits (events cut)
Processed 31 runs, 5,290 luminosity blocks and 10,787,317 events
Saving these to file rootfiles/jmenano.json for brilcalc
Analyzed 1,895,602 events
(09:06->09:19, 13 mins; 0.450/fb out of 4.86/fb)

// v4: move run, LS id earlier just in case. (TBD: Add METfilters).
// v4: 1.44/fb Golden JSON
Found 4868229 bad events according to new JSON (events cut)
Found 8864174 bad events according to trigger bits (events cut)
Processed 30 runs, 5266 luminosity blocks and 10753012 events
Saving these to file rootfiles/jmenano.json for brilcalc
Analyzed 1888838 events

// v3: add runs (and LS)  listing for Santiago+Pallabi + luminosity
// v3: add PtBins
// TBD: add TrigObj (HLT jet?)
Found 1673601 bad events according to new JSON (events cut)
Found 11450550 bad events according to trigger bits (events cut)
Processed 55 runs, 7004 luminosity blocks and 2497090 events
Analyzed 2497090 events
(17:16->17:36, 20 mins)

// v2_fwd40 tried HLT_PFJetFwd40 (HLT_PFJet40 empty)
Found 1673601 bad events according to new JSON (events cut)
Found 13819691 bad events according to trigger bits (events cut)
Found 10619 bad events not in fwd trigger phase space (events cut)
Analyzed 127949 events

// v2 added genWeight for MC (add runbx plot for data?)
// v2 moved to HLT_DiPFJetAve40 only for data (690k->490k)
// processing times 4 and 12 minutes
Found 1673601 bad events according to new JSON (events cut)
Found 13459281 bad events according to trigger bits (events cut)
Found 0 bad events not in fwd trigger phase space (events cut)
Analyzed 488359 events


// v1 data 2022-08-18 16:51 -> 16:57 (6 minutes)
Found 1673601 bad events according to new JSON (events cut)
Found 13257262 bad events according to trigger bits (events cut)
Found 51040 bad events not in fwd trigger phase space (events cut)
Analyzed 690378 events

// v1 mc 2022-08-16 17:07 -> 17:22 (15 minutes)
Loaded 7893896 entries
Analyzed 7893896 events
