DATA FORMATTING STANDARDS

While there are many pre-processing pipelines used within the lab, the code in this repository makes several assumptions.


1) All analyses are done on files formatted in the FMAToolbox standard (res, fet, clu, spk) 
	http://neurosuite.sourceforge.net/formats.html

2) All analyses are done on variables formatted in the FMAToolbox standard
	- spikes  N x 1 cell array, where each element of N is the spike times of a single unit (in seconds)
	- lfp
	- units   N x 2 matrix where column 1 is the shank ID and column 2 is the unit ID

3) If you use a preprocessing pipeline that is not in the /preprocessing/spikefileconversions folder, feel free to make and add one

4) 
