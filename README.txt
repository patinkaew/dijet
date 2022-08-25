How to RUN on Hefaistos:
------------------------
(try mosh to not drop connection)
- from local: 'rsync -rutP DijetHistosFill.C DijetHistosFill.h Hefaistos:/media/storage/dijet/'
- source /work/data/root/bin/thisroot.sh
- (rm *.d *.so *.pcm)
- root -l -b -q mk_CondFormats.C
- #define GPU in mk_DijetHistosFill.C
=> edit (version, IOV_list) and execute 'python runAllIOVs.py'
[- nohup root -l -b -q mk_DijetHistosFill.C > log.txt & ]
+ runtime about X
=> edit (version, IOV_list) and execute 'python renameAllIOVs.py'

+ tail -f log.txt
+ starting up takes quite a bit (~X sec) due to GetEntries() call
   => code TStopWatch to ignore this startup time, or skip GetEntries

- rsync -rutP files from Hefaistos

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
