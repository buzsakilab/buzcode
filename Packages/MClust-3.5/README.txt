MClust version 3.5

Version 3.5.00 4 Feb 2008 
 - ADR
Version 3.5.01 6 Feb 2008 
 - NCST added corrections to WavePC1,2,3
Version 3.5.02 9 Feb 2008 
 - NCST fixed bug in ClearWorkspace re cleaning global variables
 - NCST fixed 32-sample assumption in ShowAverageWaveform
 - NCST added ylim for plotting waveforms global variable
 - ADR created MClustResetGlobalVariables function to bring globals
 into one place 
Version 3.5.03 19 Feb 2008
 - ADR/JCJ added GetName to cluster classes so they can be called correctly
 - ADR renamed wavePC names returned 
 - ADR ClearWorkspace bug fixed
 - ADR Fixed ExportWithColorMerge bug
 - ADR added SetName to access cluster names
 - ADR Fixed CopyCluster
Version 3.5.04 20 Feb 2008
 - ADR ClearWkspace clears undo ptr
 - ADR Cleaned up undo code
 - ADR Added PeakHists
 - ADR Removed extraneous +1 in ShowHistISI
 - ADR Fixed name bug
 - ADR Fixed major SaveClusters when writing files bug
 - ADR Fixed major LoadClusters bug
Version 3.5.05 25 Feb 2008
 - ADR Fixed major export with merge bug
 - ADR Removed exit when exporting from KKwik
Version 3.5.06 26 Feb 2008
 - ADR Fixed several bugs
Version 3.5.07 27 Feb 2008
 - ADR Changed calling methodology of MClustCutter to simplify calls
 - ADR Fixed several bugs
 - ADR Cleaned up ClusterQuality
Version 3.5.08 02 Mar 2008
 - ADR Redid undo mechanism completely
 - ADR Fixed keyboard control in KlustaKwikDecisionWindow
Version 3.5.09 09 Mar 2008
 - ADR/NCST Fixed WaveformCutter undo
 - ADR/NCST cluster name wasn't being carried when convert to mcconvexhull
 - ADR removed those annoying error sounds by making undo/redo windows non-modal
 - ADR fixed undo/redo name bug
 - ADR put filename back in load features text note
 - ADR removed CreateClusterFromDotTFile
Version 3.5.10 15 Mar 2008
 - ADR fixed clear KKC bug
Version 3.5.11 1 Apr 2008
 - ADR fixed Run_KKwik bug
Version 3.5.12 11 Apr 2008
 - ADR mcconvexhull - if no limits, then no points
 - ADR mcconvexhull - changed InPolygon to correctly ID points on the boundaries

VERSION 3.5A.00 19 Apr 2009
 - ADR changed back to only loading in features as they are
       displayed. Surprisingly, this is much faster. 
       
       Time is now a true feature.  Also, Peak and Time return _Peak
       and _Time as names, thus we can no longer assume that names
       match feature files.  Also, since Time only returns a single
       feature, features do not necessarily return 4 features.  Done
       using the concept of "FeatureSources" which is a nF by 2 cell
       array of fd filenames and feature-component-within-the-files.
       This required extensive rewriting of all cluster codes, so 3.5A
       is no longer backwards compatible.

 - ADR added type to add cluster as

 - ADR added delete cluster

 - ADR fixed bug in ConvexHull with points lining up (very rare case
   unless have many many points)

 - ADR fixed bug in InPolygon with points exactly lining up (very rare
   case unless have many many hulls to match)

Version 3.5A.01 26 Apr 2008
 - ADR fixed bug in feature_Peak.  

Version 3.5A.02 1 May 2008
 - ADR fixed bug in Write_fd_file.

Version 3.5A.03 10 May 2008
 - ADR fixed display bug in KlustaKwikCallbacks
 - ADR fixed bug in ParseBatchFile in finding iCV

Version 3.5A.04 29 May 2008
 - ADR fixed FindInClusters bug in KlustaKwikCallbacks (Return Merged)

Version 3.5A.05 28 July 2008
 - ADR added EvalOverlap
 - ADR fixed bug in DeleteLimit when there was no limit to delete