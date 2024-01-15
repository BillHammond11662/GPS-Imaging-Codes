function [vupt,supt1,supt2,supt3,dwin,iloc,igps,loni,lati,w,tri]=msf2pt(lon,lat,v,sv,tri0,latpt,lonpt,ssf,opt,optsphere,optrnf)
%  [vupt,supt1,supt2,supt3,dwin,iloc,igps,loni,lati,w,tri]=msf2pt(lon,lat,v,sv,tri0,latpt,lonpt,ssf,opt,optsphere,optrnf)
%
%  INPUTS
%
%  latpt, lonpt  scalar, describe location of evaluation point
%  lat, lon      vectors, describe list of input data locations of
%  tri0          Delaunay triangulation of input lon,lat
%  v and sv      input data of values and their uncertainties
%  ssf           the spatial structure function.  nx2 where column 1 is distances,
%                column 2 is values, generally decreasing with distance. Distances 0 and
%                Inf are implicity assumed to be before begninning and after end of
%                vector.
%  opt           1 = uses weighted median, 2 = uses regular median, 3 = weighted mean
%  
%  optsphere     1 uses Delauny traingulation on a sphere
%
%  optrnf = 1 uses robust network filter option of kreemer et al., 2020
%  optrnf = 0 does not.
%
%  OUTPUTS
%
%  vupt      estimate of vu at latpt, lonpt
%  iloc 	 first local group of lat,lon
%  igps 	 the list of indeces specifing gps stations used for the estimate (indeces of lat,lon)
%  lati,loni coordinates of triangles used in estimate
%  w         is list of weights used for points lati,loni
%  tri       the delaunay triangulation
%  dwin      radius inside which stations are used to estimate local value
%
%  uncertainty computed in three different ways
%    supt1 = based on formal uncertainty of mean
%    supt2 = based on RMS of residual scatter of contributing points
%    supt3 = based on robust standard deviation of residual scatter (mad times 1.4826)
%
% Bill Hammond 2023-01-15
% University of Nevada, Reno
% Hammond et al., JGR 2016, doi:10.1002/2016JB013458.

supt1=nan;
supt2=nan;
dwin =[];
igps =[];
loni=[];
lati=[];
w=[];
tri=[];

% Shorten the list of stations to be level 2 in Delaunay triangulation from eval point.
% but don't do it if optrnf == 1 since that requires preserving all
% stations close to evaluation point, not just the nearest delaunay
% connected
% using optrnf==1 will therefore take longer to run
% but will incorporate more data near the evalation point thus weighting
% local data more


if optrnf~=1
    [dist,~]=baz(lat,lon,latpt,lonpt);
    [~,imin]=min(dist);
    [il1,~]= find(tri0==imin);
    k1=unique(tri0(il1,:));
    k2=[];
    for k=k1'
        [il2,~]=find(tri0==k);
        k2=[k2;unique(tri0(il2,:))];
    end
    iloc=unique(k2);
else
    iloc=1:length(lat);
    iloc=iloc';
end
if length(iloc)<=2
    vupt = NaN;
    return;
end
loni=lon(iloc);
lati=lat(iloc);
vi=v(iloc);
svi=sv(iloc);

%if eval point is exactly one of the lat/lons put it at the end of the list.
isptinlist=find(lonpt==loni & latpt==lati);
if ~isempty(isptinlist)
    loni=[loni;lonpt];
    lati=[lati;latpt];
    vi=[vi;vi(isptinlist)];
    svi=[svi;svi(isptinlist)];
    lati(isptinlist)=[];
    loni(isptinlist)=[];
    vi(isptinlist)=[];
    svi(isptinlist)=[];
    iloc=[iloc;iloc(isptinlist)];
    iloc(isptinlist)=[];
else
    loni=[loni;lonpt];
    lati=[lati;latpt];
end
q=length(loni);

% checks that iloc was reorganized correctly
% if ~isempty(find(loni~=lon(iloc), 1))
%     keyboard;
% end

if optsphere==1
    [x,y,z,~]=latlon2xyz(lati,loni,zeros(size(loni)),[]);
    n=size(lati,1);
    [~,tri] = sphere_delaunay(n,[x';y';z']);
    tri=tri';
else
    tri=delaunay(loni,lati);
end

% hold on;
% triplot(tri,loni,lati);

% get the set of points connected to the evaluation point 
[jt,~] = find(tri==q);
[iw,~] = unique(tri(jt,:));

% this is the robust network filter of kreemer et al., 2020
% find all stations within median distance of selected points.
if optrnf==1  
    [dista,~]=baz(latpt,lonpt,lati,loni);
    [distb,~]=baz(latpt,lonpt,lati(iw),loni(iw));
    md = median(distb);
    k=find(dista<=md);
    if ~isempty(setdiff(k,iw))
        iw = [iw;setdiff(k,iw)];
        iw = unique(iw);
    end
end

% if evaluation point is not a gps station
% remove it from the contributing points since there is
% no gps data there
if isempty(isptinlist)
   iw=setdiff(iw,q); 
end

[n,m]=size(iw);
if min([n;m])~=1
    error('iw is wrong shape');
end
if m>1
    iw=iw';
end

igps=iloc(iw);

% figure(2);
% clf;
% plot(loni,lati,'ko');
% hold on;
% triplot(tri,loni,lati);
% plot(lonpt,latpt,'r*')
    
[dist,~]=baz(latpt,lonpt,lati(iw),loni(iw));

% vmsf=calcv(tri,i,lat,lon,v,sv,ssf,opt,optrnf);

if ~isempty(iw)
    if opt==1 || opt ==3
        if isempty(ssf)
            w=1./dist;
        else
            w=interp1(ssf(:,1),ssf(:,2),dist,'linear');
        end
        w = w./svi(iw);
        if all(w==0)
            w=ones(size(w));
        end
        w = w/sum(w);
        if any(isnan(w))
            disp('Getting NaNs in w');
            keyboard;
        end
        if opt==1
            vupt = weightedMedian(vi(iw),w);
        elseif opt==3
            vupt = sum(vi(iw).*w);
        end
        supt1 = sqrt(sum((w.*svi(iw)).^2));
        supt2 = sqrt(sum((vi(iw)-vupt).^2)/length(iw));
        supt3 = mad(vi(iw)-vupt)*1.4826;
    elseif opt==2
        vupt=median(vi(iw));
        supt1 = sqrt(sum(svi(iw).^2)/sqrt(length(iw)));
        supt2 = sqrt(sum((vi(iw)-vupt).^2)/length(iw));
        supt3 = mad(vi(iw)-vupt)*1.4826;
    else
        error('Unrecognized option');
    end
end



