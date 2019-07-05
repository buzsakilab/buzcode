function funcsynapses = addStartEndsToFuncSynapsesFromOtherFuncSynapses(funcsynapses,fs2)

funcsynapses.CnxnStartTimesVsRefSpk = fs2.CnxnStartTimesVsRefSpk;
funcsynapses.CnxnStartBinsVsRefSpk = fs2.CnxnStartBinsVsRefSpk;
funcsynapses.CnxnStartTimesVsCCGStart = fs2.CnxnStartTimesVsCCGStart;
funcsynapses.CnxnStartBinsVsCCGStart = fs2.CnxnStartBinsVsCCGStart;
funcsynapses.CnxnEndTimesVsRefSpk = fs2.CnxnEndTimesVsRefSpk;
funcsynapses.CnxnEndBinsVsRefSpk = fs2.CnxnEndBinsVsRefSpk;
funcsynapses.CnxnEndTimesVsCCGStart = fs2.CnxnEndTimesVsCCGStart;
funcsynapses.CnxnEndBinsVsCCGStart = fs2.CnxnEndBinsVsCCGStart;

funcsynapses.ZeroLag.CnxnStartTimesVsRefSpk = fs2.ZeroLag.CnxnStartTimesVsRefSpk;
funcsynapses.ZeroLag.CnxnStartBinsVsRefSpk = fs2.ZeroLag.CnxnStartBinsVsRefSpk;
funcsynapses.ZeroLag.CnxnStartTimesVsCCGStart = fs2.ZeroLag.CnxnStartTimesVsCCGStart;
funcsynapses.ZeroLag.CnxnStartBinsVsCCGStart = fs2.ZeroLag.CnxnStartBinsVsCCGStart;
funcsynapses.ZeroLag.CnxnEndTimesVsRefSpk = fs2.ZeroLag.CnxnEndTimesVsRefSpk;
funcsynapses.ZeroLag.CnxnEndBinsVsRefSpk = fs2.ZeroLag.CnxnEndBinsVsRefSpk;
funcsynapses.ZeroLag.CnxnEndTimesVsCCGStart = fs2.ZeroLag.CnxnEndTimesVsCCGStart;
funcsynapses.ZeroLag.CnxnEndBinsVsCCGStart = fs2.ZeroLag.CnxnEndBinsVsCCGStart;