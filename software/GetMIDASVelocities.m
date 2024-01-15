function [sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xenu50,rho,enuofrac,enusigvpair,nsteps]=GetMIDASVelocities(frame,velfile)
% [sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xenu50,rho,enuofrac,enusigvpair,nsteps]=GetMIDASVelocities(frame,velfile)
%
% frame = 'IGS14' or 'NA' or one of the other plate frames (e.g. PA, EU, ...)
% see http://geodesy.unr.edu/velocities/midas.readme.txt for a list
%
% if velfile is empty it gets the file from http://geodesy.unr.edu/velocities/
% otherwise it uses velfile from local directory 
% 
% sta         - station 4-character ID
% lat, lon, h - station latitude, longitude and height of station
% lab         - label indicating version of midas used
% t1, tm      - time series first and last epochs
% dT          - duration of time series in years (difference between first and last)
% m           - number of epochs of data, used or not
% ngood       - number of epochs of good data, i.e. used in at least one velocity sample
% ve,vn,vu    - east, north, up mode velocities (m/yr)
% se,sn,su    - east, north, up mode velocity uncertainties (m/yr)
% xenu50      - east, north, up offset at at first epoch (m)
% rho         - mean number of epochs per day (mean temporal density, 1 = no missing data)
% enuofrac    - east, north, up fraction of outliers
% enusigvpair - east, north, up standard deviation velocity pairs
% nsteps      - number of steps assumed, determined from our steps database

if nargin<=1
    velfile=[];
end

if isempty(velfile)
    disp('Reading MIDAS velocities from http://geodesy.unr.edu');
    velfile = ['http://geodesy.unr.edu/velocities/midas.' char(frame) '.txt'];
    [S,stat]=urlread(velfile);
    if stat==1
        C=textscan(S,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    else
        error(['Could not get midas velocity file: ' velfile]);
    end
else
    disp(['Reading MIDAS velocities from ' velfile]);
    fid=fopen(velfile,'r');
    if fid==-1
        error(['Could not read midas velocity file: ' velfile ' from local directory.']);
    end
    C=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);
end

sta=C{1};
lab=C{2};
t1=C{3};
tm=C{4};
dT=C{5};
m=C{6};
ngood=C{7};
numsol=C{8};
ve=1000*C{9};
vn=1000*C{10};
vu=1000*C{11};
se=1000*C{12};
sn=1000*C{13};
su=1000*C{14};
lat = C{25};
lon = C{26};
h = C{27};

xenu50(:,1)=C{15};
xenu50(:,2)=C{16};
xenu50(:,3)=C{17};

enuofrac(:,1)=C{18};
enuofrac(:,2)=C{19};
enuofrac(:,3)=C{20};

enusigvpair(:,1)=C{21};
enusigvpair(:,2)=C{22};
enusigvpair(:,3)=C{23};

nsteps=C{24};

rho = m./(dT*365);






