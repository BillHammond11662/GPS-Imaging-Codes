function [vmsf]=msf(lon,lat,v,sv,ssf,opt,optsphere,optpar,optrnf)
% function [vmsf]=msf(lon,lat,v,sv,ssf,opt,optsphere,optpar,optrnf)
%
% Calculates the spatially filtered velocity field.
%
% lat,lon are coordinates of GPS stations
% v, sv are values at stations and uncertainties
%
% ssf is the spatial structure function.  nx2 where column 1 is distances in degrees,
% column 2 is values, generally decreasing with distance. 
%
% if ssf is empty uses 1/distance weighting
%
% opt=1 is weighted median
% opt=2 is ordinary median filter
% opt=3 is ordinary mean
%
% optsphere = 0 uses regular Delaunay triangulation
% oprsphere = 1 uses Delaunay triangulation on a sphere (slower)
%
% optpar = 1 runs in parallel using matlab parallel computing toolbox
% optpar = 0 runs in sequential mode (if you don't have toolbox)
%
% optrnf = 1 uses robust network filter option of kreemer et al., 2020
% optrnf = 0 does not
%
% Bill Hammond 2023-01-15
% University of Nevada, Reno
% Hammond et al., JGR 2016, doi:10.1002/2016JB013458.

% check ssf
if ~isempty(ssf)
    [~,b]=size(ssf);
    if b~=2
        error('If using ssf, it must be nx2');
    end
    if any(diff(ssf(:,1))<0)
        error('If using ssf, distances must increase monotonically.');
    end
end

% remove zero length baselines by taking median of zero baseline clusters
lon0=lon;
lat0=lat;

still=1;
i=0;
lat2=[];
lon2=[];
v2=[];
sv2=[];
kmap = {};
while still
    i=i+1;
    
    j=find(lon(1)==lon & lat(1)==lat);
    k=find(lon(1)==lon0 & lat(1)==lat0);
    kmap{i,1}=k;  % so v,sv,lon,lat can be reconstructed
  
    lon2=[lon2;lon(j(1))];
    lat2=[lat2;lat(j(1))];
    v2=[v2;median(v(j))];
    sv2=[sv2;min(sv(j))];
    
    lon(j)=[];
    lat(j)=[];
    v(j)=[];
    sv(j)=[];
    if isempty(lon)
        still=0;
    end
end

if optsphere==1
    [x,y,z,~]=latlon2xyz(lat2,lon2,zeros(size(lon2)),[]);
    n=size(lat2,1);
    [~,tri] = sphere_delaunay(n,[x';y';z']);
    tri=tri';
else
    tri=delaunay(lon2,lat2);
end

vmsf2=nan(size(v2));

if optpar==1
    parfor i=1:length(lon2)
        [vmsf2(i)]=calcv(tri,i,lat2,lon2,v2,sv2,ssf,opt,optrnf);
    end
else
    for i=1:length(lon2) 
        [vmsf2(i)]=calcv(tri,i,lat2,lon2,v2,sv2,ssf,opt,optrnf);     
    end   
end

% now map v2 back to v
vmsf=nan(size(v));
for i=1:length(kmap)
    vmsf(kmap{i},1)=vmsf2(i);
end
