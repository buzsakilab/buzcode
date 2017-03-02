DATA FORMATTING STANDARDS
(under development: read and contribute to the wiki!)
https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards

While there are many pre-processing pipelines used within the lab, the code in this repository makes several assumptions.


1) All analyses are done on files formatted in the FMAToolbox standard (res, fet, clu, spk) 
	http://neurosuite.sourceforge.net/formats.html

2) All analyses are done on variables formatted in the FMAToolbox standard
	- spikes  N x 1 cell array, where each element of N is the spike times of a single unit (in seconds)
	- lfp
	- units   N x 2 matrix where column 1 is the shank ID and column 2 is the unit ID
	- pos	N X M matrix where N is the number of frames captured during a recording, and M is a 			set of columns (5 for LED 11 for Optitrack)
			1-time stamps that match recording time,
			2-X positions,
			3-Y positions,
			4-Z positions, (2nd LED-X positions) 
			5-X rotation,  (2nd LED-Y positions)
			6-Y rotation,
			7-Z rotation,
			8-W rotation,
			9-Error Per Marker,
			10-Frame Count,
			11-Frame Time (Optitrack clock) 
	- events struct that is common to all event types (ripples, spindles, etc);	
		.starts
		.peaks
		.stops
		.peakAmplitudes
		.detectorMethod (full path to function used, example: /buzcode/externalPackages/FMAToolbox/Analyses/FindRipples.m)
		.detectorArgs.inputs (extra inputs given to detectorMethod)

3) If you use a preprocessing pipeline that is not in the /preprocessing/spikefileconversions folder, feel free to make and add one

4) 
