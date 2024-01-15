function [VMSF,SVMF]=msf2grid(lon,lat,v,sv,LON,LAT,INNA,ssf,opt,optsphere,optpar,optrnf)
% function [VMSF,SVMF]=msf2grid(lon,lat,v,sv,LON,LAT,INNA,ssf,opt,optsphere,optpar,optrnf)
%
% calculates the weighted median spatial filter of the velocity field
% onto a grid
%
% lon,lat,v,sv are input data list of longitudes, latitudes, values, and
% value uncertainties respectivly 
%
% LON,LAT are matrices with evaluation point longitude and latitude respectively, 
% generated using e.g., "meshgrid"
%
% ssf is the spatial structure function.  nx2 where column 1 is distances,
% column 2 is values, generally decreasing with distance. Distances 0 and
% Inf are implicity assumed to be before begninning and after end of
% vector.
%
% opt=1 is weighted median
% opt=2 is ordinary median filter
% opt=3 is weighted mean
%
% optsphere = 1 uses Delaunay Triangulation on a sphere, otherwise cartesian
%
% optpar = 1 runs in matlab parallel mode
%
% INNA is mask, when value is nan then the grid value is not calculated
%
% VMSF is the output filtered and interpolated rate field. Same dimensions as LON,LAT,INNA
% SVMF is the same dimension as VMSF by 3.  It has three different measures
% of uncertainty in the field.  See msf2pt.m for definitions. 
%
% Bill Hammond 2023-01-15
% University of Nevada, Reno
% Hammond et al., JGR 2016, doi:10.1002/2016JB013458.

% shorten lat/lon list to include only values inside LON,LAT plus a 25% buffer
dd = max([max(max(LON))-min(min(LON))  max(max(LAT))-min(min(LAT))])/4;

iout=find(lon<(min(min(LON))-dd) | lon>(max(max(LON))+dd) | ...
          lat<(min(min(LAT))-dd) | lat>(max(max(LAT))+dd));
lon(iout)=[];
lat(iout)=[];
v(iout)=[];
sv(iout)=[];

[n,m]=size(LON);
VMSF=nan(n,m);
SVMF=nan(n,m,3);

dp=fix(n/20);
steps=dp:dp:n;
it=1;

if optsphere==1
    [x,y,z,~]=latlon2xyz(lat,lon,zeros(size(lon)),[]);
    nn=size(lat,1);
    [~,tri0] = sphere_delaunay(nn,[x';y';z']);
    tri0=tri0';
else
    tri0=delaunay(lon,lat);
end

if optpar==1
       
    for i=1:n
        svmf1=nan(1,m);
        svmf2=nan(1,m);
        svmf3=nan(1,m);
       
        if ~isempty(find(i==steps, 1))
            if it==1
                fmt='%3.0f%%\n';
            else
                fmt='\b\b\b\b%3.0f%%\n';
            end
           fprintf(fmt,100*i/n);
           it=it+1;
        end
        
        parfor j=1:m
            
            if INNA(i,j)==1
                
                [VMSF(i,j),svmf1(j),svmf2(j),~,~,~,~,~,~,~]=...
                    msf2pt(lon,lat,v,sv,tri0,LAT(i,j),LON(i,j),ssf,opt,optsphere,optrnf);
                
            end
        end
        SVMF(i,:,1)=svmf1;
        SVMF(i,:,2)=svmf2;
        SVMF(i,:,3)=svmf3;
    end
    
else
        
    for i=1:n
        
        if ~isempty(find(i==steps, 1))
            if it==1
                fmt='%3.0f%%\n';
            else
                fmt='\b\b\b\b%3.0f%%\n';
            end
           fprintf(fmt,100*i/n);
           it=it+1;
        end
        
        for j=1:m
            
            if INNA(i,j)==1
                
                [VMSF(i,j),SVMF(i,j,1),SVMF(i,j,2),~,~,~,~,~,~,~]=...
                msf2pt(lon,lat,v,sv,tri0,LAT(i,j),LON(i,j),ssf,opt,optsphere,optrnf);

%                 if isnan(VMSF(i,j))
%                     keyboard;
%                 end
                
            end
        end
    end
end






