# Sleep-State-Score

File Structure
datasetfolder: the top level folder of your dataset, most likely holds a folder for each of a number of recordings
recordingname: the name for a given recording.  All output files will be of the form /datasetfolder/recordingname/recordingname_file.ext

Necessary Files
/datasetfolder/recordingname/recordingname.lfp
/datasetfolder/recordingname/recordingname.xml

Important Steps:
Identify any "bad" channels.  For example particularly noisy or very low amplitude channels.  Put these in a .txt document: datasetfolder/recordingname/bad_channels.txt.
-If channel noise appears in the middle of the recording, you may need to run the state scoring, identify bad channels, and then run state scoring again.
-bad_channels are 0-indexed, as they are in Neuroscope.

-output channels are also 0-indexed, as in Neuroscope
