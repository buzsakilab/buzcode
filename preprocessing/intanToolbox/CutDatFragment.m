% function CutDatFragment(fname,timeperiod)
% Creating a new dat file from a fragment of an original dat file
% INPUTs
% fname: filename of original dat
% timeperiod: start and stop of file fragment to keep
%
%
% Original by Antonio Ruiz
% Modified by Peter Petersen

fname = 'Peter_MS13_171130_143907';
timeperiod = [0,44*60+9];

Intan_rec_info = read_Intan_RHD2000_file_Peter(pwd);

ch = length(Intan_rec_info.amplifier_channels);
sr = Intan_rec_info.frequency_parameters.amplifier_sample_rate;
a = memmapfile([fname '.dat'],'Format','int16');
aa = a.data;
fid = fopen([fname '_fragment.dat'],'W');
for i = timeperiod(1)+1:timeperiod(2)
     fwrite(fid,aa((i-1)*ch*sr+1:i*ch*sr),'int16');
end
fclose(fid);
clear aa a

disp('Writing aux file')
ch = length(Intan_rec_info.aux_input_channels);
h1 = fopen('auxiliary_fragment.dat','W');
m = memmapfile('auxiliary.dat','Format','uint16','writable',false);
fwrite(h1,m.Data(timeperiod(1)*ch*sr+1:timeperiod(end)*ch*sr),'uint16');
fclose(h1);

disp('Writing analogin file')
ch = length(Intan_rec_info.board_adc_channels);
h2 = fopen('analogin_fragment.dat','W');
m = memmapfile('analogin.dat','Format','uint16','writable',false);
fwrite(h2,m.Data(timeperiod(1)*ch*sr+1:timeperiod(end)*ch*sr),'uint16');
fclose(h2);

disp('Writing digitalin file')
h3 = fopen('digitalin_fragment.dat','W');
m = memmapfile('digitalin.dat','Format','int16','writable',false);
fwrite(h3,m.Data(timeperiod(1)*sr+1:timeperiod(end)*sr),'int16');
fclose(h3);

disp('Writing time file')
h4 = fopen('time_fragment.dat','W');
m = memmapfile('time.dat','Format','int16','writable',false);
fwrite(h4,m.Data(timeperiod(1)*sr+1:timeperiod(end)*sr),'int16');
fclose(h4);

disp('Writing supply file')
h5 = fopen('supply_fragment.dat','W');
m = memmapfile('supply.dat','Format','int16','writable',false);
fwrite(h5,m.Data(timeperiod(1)*sr+1:timeperiod(end)*sr),'int16');
fclose(h5);

