function MAfilt = MAfilt(rg,deeg);

fn = [0.01:0.01:0.98];

cfirpm(50,fn,1./fn);
deegF = filtfilt(20,fn,1./fn)

