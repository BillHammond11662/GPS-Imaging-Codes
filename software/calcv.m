function vmsf=calcv(tri,i,lat,lon,v,sv,ssf,opt,optrnf)
%  vmsf=calcv(tri,i,lat,lon,v,sv,ssf,opt,optrnf)
%
%  Inputs:
%  lon,lat  - longitude latitudes of entire network
%  vu, sv   - velocity and uncertainty, same dimensions as lon,lat
%  ssf      - spatial structure function 
%  tri      - Delaunay triangulation of points in lon,lat
%  i        - is index of lon,lat we are working on now
%  opt      - specifies weighted median, ordinary median filter or mean
%  optrnf   - specifies wether to use robust network filter option
%
% Bill Hammond 2023-01-15
% University of Nevada, Reno
% Hammond et al., JGR 2016, doi:10.1002/2016JB013458.

if ~isempty(ssf)
    if ssf(1,1)~=0 && opt==1
        disp('Warning: ssf does not extend to zero distance');
    end
end

vmsf=NaN;

[jt,~] = find(tri==i);
[iw,~] = unique(tri(jt,:));


if optrnf==1  % this is the robust network filter of kreemer et al., 2020
    %find all stations within median distance of selected points.
    [dista,~]=baz(lat(i),lon(i),lat,lon);
    [distb,~]=baz(lat(i),lon(i),lat(iw),lon(iw));
    md = median(distb);
    k=find(dista<=md);
    if ~isempty(setdiff(k,iw))
        iw = [iw;setdiff(k,iw)];
        iw = unique(iw);
    end
end


[dist,~]=baz(lat(i),lon(i),lat(iw),lon(iw));

if isempty(iw)
    return;
end

if opt==1 || opt==3
    w=interp1(ssf(:,1),ssf(:,2),dist,'linear')./sv(iw);
    if all(w==0)
        w=ones(size(w));
    end
    w=w/sum(w);
    if opt==1
        vmsf = weightedMedian(v(iw),w);
    elseif opt==3
        vmsf = sum(v(iw).*w);
    end
elseif opt==2
    vmsf=median(v(iw));
    w = ones(size(iw));
    w = w/length(w);
else
    error('Unrecognized option for opt');
end

if any(isnan(w))
    disp('Problem: we have NaNs in w');
    keyboard;
end

if isnan(vmsf)
    disp('Something went wrong and vmsf=NaN');
    keyboard
end

