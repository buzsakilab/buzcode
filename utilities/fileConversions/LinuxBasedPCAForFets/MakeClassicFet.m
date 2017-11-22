function MakeClassicFet(basename,basePath)

%get basic input params if not there
if ~exist('basename','var')
    [~,basename] = fileparts(cd);
end
if ~exist('dirname','var')
    basePath = cd;
end

cd(basePath)

% if already .fets, move them to a new folder called "PreviousFets"
d = dir(fullfile(basePath,[basename '.fet.*']));
if ~isempty(d);
    mkdir(fullfile(basePath,'PreviousFets'))
    for a = 1:length(d);
        movefile(fullfile(basePath,d(a).name),fullfile(basePath,'PreviousFets',d(a).name));
    end
end

%if ~exist([basename '.fil'],'file')
%    cmd = ['process_mhipass ' basename];
%    system(cmd)
%end

if ~exist([basename '.fil'])
    %xml =  LoadParameters([basename '.xml']);
    xml = bz_getSessionInfo(basePath);
    inname = [basename '.dat'];
    outname = [basename '.fil'];
    numch = num2str(xml.nChannels);
    sampl = num2str(xml.rates.wideband);
    lowband = '800';
    highband = '9500';
    forder = '50';
    gain = '1';
    offset = '0';
    
    firfilter(inname,outname,numch,sampl,lowband,highband,forder,gain,offset);
end    


%get xml
if ~exist('xml','var')
    xml = bz_getSessionInfo(basePath);
end
nGrps = length(xml.SpkGrps);


% look for process_pca_multi or ndm_pca
if ~system('which ndm_pca')
    pca = 'ndm_pca ';
elseif ~system('which process_pca_multi')
    pca = 'process_pca_multi ';
end


for g = 1:nGrps
%     cmd = ['process_pca_multi ' basename ' ' num2str(g)];
    cmd = [pca basename ' ' num2str(g)];
    disp(cmd)
    a =  system(cmd)
end


if a ==0
    cmd = ['rm ' basename '.fil'];
    system(cmd)
end

