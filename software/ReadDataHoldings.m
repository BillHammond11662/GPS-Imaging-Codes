function [sta,lat,lon,h,x,y,z,dtbeg,dtend,dtmod,nsol,staorig,nstao,dT] = ReadDataHoldings(dhfile)
% [sta,lat,lon,h,x,y,z,dtbeg,dtend,dtmod,nsol,staorig,nstao,dT] = ReadDataHoldings(dhfile)

if nargin==0
    dhfile=[];
end


if isempty(dhfile)
    disp('Reading DataHoldings.txt from http://geodesy.unr.edu');
    S=webread('http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt');
    C=textscan(S,'%s','Headerlines',1,'Delimiter','\n\r');
else
    disp(['Reading DataHoldings.txt from ' dhfile]);
    fid=fopen(dhfile,'r');
    if fid<0
        disp(['Could not open ' dhfile '.  Quitting.']);
    end
    C=textscan(fid,'%s','Headerlines',1,'Delimiter','\n\r');
    fclose(fid);
end

A=C{1};
N= size(A,1);
sta=cell(N,1);
lat=nan(N,1);
lon=nan(N,1);
h=nan(N,1);
x=nan(N,1);
y=nan(N,1);
z=nan(N,1);
dtbeg=cell(N,1);
dtend=cell(N,1);
dtmod=cell(N,1);
nsol=nan(N,1);
staorig=cell(N,1);
nstao=zeros(N,1);
    
for i=1:N
    B=char(A(i));
    
    len=length(B);
    
    sta{i}=B(1:4);
    lat(i)=str2double(B(6:13));
    lon(i)=str2double(B(15:24));
    h(i)=str2double(B(26:33));
    x(i)=str2double(B(35:47));
    y(i)=str2double(B(49:61));
    z(i)=str2double(B(63:75));
    dtbeg{i}=B(77:86);
    dtend{i}=B(88:97);
    dtmod{i}=B(99:108);
    nsol(i)=str2double(B(110:115));

    if len>115
        still=1;
        slist = cell(0);
        j=1;
        while still
           ib = 116+[((j-1)*4+1):(j*4)]+j-1;
           slist{j}=B(ib);
           j=j+1;
           if len<(ib(end)+5)
               still=0;
           end
        end
        nstao(i)=j-1;
        staorig{i}=slist;
    end
    
end

lon=mod(lon,360);
lon(lon>180)=lon(lon>180)-360;

dT=nan(size(sta));
for i=1:length(sta)
    nbeg = datenum(dtbeg(i),'yyyy-mm-dd');
    nend = datenum(dtend(i),'yyyy-mm-dd');
    dT(i)=(nend-nbeg)/365.25;
end

end
