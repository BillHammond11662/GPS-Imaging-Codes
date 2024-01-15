
clear
close all;
tic;                

optpar = 0;     % Determines whether Matlab will run in parallel mode.
                % Can only use if you have the parallel computing toolbox.
                % If you do not have the parallel toolbox use optpar=0;
                % It can sometimes take longer in parallel mode owing to
                % overhead.
                
velfile=[];     % MIDAS (or any) velocity file
                % setting this to empty set will automatically access the version 
                % of the MIDAS velocity field that is online at the Nevada
                % Geodetic Laboratory website e.g., http://geodesy.unr.edu/velocities/midas.IGS14.txt
                % This file is updated weekly.  If you have downloaded a
                % version of the midas velocity file previously and want
                % to use it, set velfile='/the/path/and/name/of/the/midas/velocity/file.txt'
                % or you can use any velocity file you wish, as long as
                % it is populates the same variable names when it is read
                % in.
                
frame='IGS14';  % At this time should just be 'IGS14' (global frame). It is also OK to use any other of the plate 
                % reference frames listed at http://geodesy.unr.edu under the "MIDAS Velocity Fields".
                % However, stations are only available in these files when
                % they are on the specified plate.
                % So it is easeist to just use 'IGS14' since all stations are 
                % available in this file, and for the vertical component 
                % all plate reference frames have identical time series anyway.
                % However, if using this on horizontal velocities then the plate
                % frame will have to be chosen carefully. 
                % This is only needed if velfile=[];
                
bounds = [-125 -114 32.5 42];    
                % lon,lat bounds of domain of interest
                % format is [lowlon highlon lowlat hightlat]

dTmin=2.5;      % Minimium time series duration used. Omits time series with shorter durations

sigmax=2;       % Maximum velocity uncertainty to accept. Omits velocities with larger uncertainties. 

dvmax= 25;      % For SSF construction, maximum velocity difference allowed. 
                % Can help.  Making it very large (e.g. 1000) causes it to have no effect.

opt = 1;        % opt=1 uses weighted median, opt=2 uses ordinary median, opt=3 uses weighted mean

optsphere=0;    % Uses delauanay triangulation on a sphere (Renka, 1997)
                % spherical delaunay has substantially longer run time,
                % but will provide better results when looking at very large
                % areas with irregular station spacing or near poles. 
                % If you want to use this you can use codes available here:
                % http://people.sc.fsu.edu/~jburkardt/m_src/sphere_delaunay/sphere_delaunay.html
                % get "sphere_delaunay.m". Make sure its folder is in
                % your path and set optsphere=1

optrnf=0;       % optrnf==1 uses Robust Network Filter option of Kreemer et al., 2020
                % incorporates points within median distance of points 
                % connected to evalution point in Delaunay triangulation.
                % This increases weight of local data, especially in cases where
                % GPS network is dense on one side, sparse on the other.
                % Increases run time. 

usemask = 1;    % if you want to mask out some areas, like oceans or areas 
                % with few GPS stations usemask=1 and do some work in the
                % section that makes the mask variable INNA below.  This
                % example constructs INNA to work for Californian and
                % Nevada.

optiterate=0;   % If optiterate == 1 then an iterative resampling of the velocity field
                % and ssf is performed, deriving a median field from many iterations of the
                % imaging of the data. This helps reduce minor outliers associated with the imaging algorithm. 
%                
%
% Bill Hammond 2023-01-15
% University of Nevada, Reno
% Hammond et al., JGR 2016, doi:10.1002/2016JB013458.

%%  Get GPS velocities, selcted those to use, decluster stations
%
% use this if you want to use MIDAS rates estimated by NGL

disp('Getting GPS velocities from NGL MIDAS file...');

[sta,lat,lon,~,~,~,~,dT,~,~,~,~,~,vu,~,~,su,~,~,~,~,~]=...
    GetMIDASVelocities(frame,velfile);

iout=find( su>sigmax | dT<dTmin | lon<bounds(1) | lon>bounds(2) | lat<bounds(3)| lat>bounds(4) | isnan(vu) | isnan(lon) | isnan(lat));

sta(iout)=[];
lon(iout)=[];
lat(iout)=[];
vu(iout)=[];
su(iout)=[];


%% Declustering of GPS stations

disp('Declustering GPS station locations...');

dcol = 0.5;  % distance to define cluster radius in km
iclusters=decluster_sta(lat,lon,dcol,0)';

latd = nan(size(iclusters));
lond=latd;
vud =latd;
sud =latd;
stad=cell(length(iclusters),1);

for i=1:length(iclusters)
   j=iclusters{i};
%    if length(j)>1
%        keyboard;
%    end
   latd(i) = median(lat(j));
   lond(i) = median(lon(j));
   vud(i) = median(vu(j)); 
   sud(i) = norm(su(j));  % negotiable
   stad{i}=[];
   for k=1:length(j)
       if k==1
           stad(i)=sta(j(k));
       else
           stad{i} = [stad{i} '_' char(sta(j(k)))];
       end
   end
end
sta = stad;
lon=lond;
lat=latd;
vu=vud;
su=sud;

N=length(sta);
disp(['After declustering have ' num2str(N) ' stations.']);

% For other areas this will need to be changed or omitted. 
% It gets the data to plot state lines on the figure. 
% However, it is also used to make a mask below, which specifies areas
% that will be ignored by the imaging algorithm.  That is useful e.g., when some of the
% area within the bounds is ocean, and you don't want to consider it.
states = shaperead('usastatehi');
xCA = states(5).X;
yCA = states(5).Y;
xNV = states(28).X;
yNV = states(28).Y;

figure(1);
clf;
plot(xCA,yCA,'k-');
hold on;
plot(xNV,yNV,'k-');
plot(lon,lat,'k.');
set(gca,'dataaspectratio',[1 cosd(38) 1]);
axis(bounds)
title('GPS Station Locations')
drawnow

%%  Make the Spatial Structure function

disp('Making the Spatial Structure Function (SSF)...');
    
if optpar==1
    delete(gcp('nocreate'))
    pool=parpool;
end

% See Hammond et al., 2016 JGR for explanation of the SSF.
% It makes the matrix "ssf" and can be run separately 
% beforehand, as it can take a while to run if there are lots of GPS stations. 
% MakeSSF is robust, and insensitive to outliers, but outliers can be
% eliminated from vu prior to running this if desired.  Using dvmax also
% helps remove large outliers from the dataset during SSF construction.

ssf=MakeSSF(lon,lat,vu,su,0,dvmax); 



%%

MakeColorMap    % creates red-white-blue colormap with white at v=0, 
                % saturating at -3 (blue) and +3 (red) mm/yr.
                
figure(3);
clf;
plot(xCA,yCA,'k-');
hold on;
plot(xNV,yNV,'k-');

sz = 7; 
for i=1:length(lon)
    iclr = fix(m*vu(i) + b);
    iclr(iclr<1)=1;
    iclr(iclr>NC)=NC;
    hp=plot(lon(i),lat(i),'o','markersize',sz,'markerfacecolor',cmap(iclr,:),'markeredgecolor',[.5 .5 .5]);
end
colormap(cmap);
clim([vumin vumax]);
colorbar('vert');
xlabel('Longitude');
ylabel('Latitude');
axis(bounds);
set(gca,'dataaspectratio',[1 median(cosd(lat)) 1]);
title('MIDAS Vertical Rates');


%%  

disp('Performing Preliminary Median Spatial Filtering...');

vumsf=msf(lon,lat,vu,su,ssf,opt,optsphere,optpar,optrnf);

figure(4);
clf;
plot(xCA,yCA,'k-');
hold on;
plot(xNV,yNV,'k-');

for i=1:length(lon)
    iclr = fix(m*vumsf(i) + b);
    iclr(iclr<1)=1;
    iclr(iclr>NC)=NC;
    hp=plot(lon(i),lat(i),'o','markersize',sz,'markerfacecolor',cmap(iclr,:),'markeredgecolor',[.5 .5 .5]);
end
clim([vumin vumax]);
colormap(cmap);
colorbar('vert');
xlabel('Longitude');
ylabel('Latitude');
title('Median Spatial Filtered Vertical Rates');

axis(bounds);
set(gca,'dataaspectratio',[1 median(cosd(lat)) 1]);

figure(5);
clf;
plot(xCA,yCA,'k-');
hold on;
plot(xNV,yNV,'k-');

for i=1:length(lon)
    iclr = fix(m*(vu(i)-vumsf(i)) + b);
    iclr(iclr<1)=1;
    iclr(iclr>NC)=NC;
    hp=plot(lon(i),lat(i),'o','markersize',sz,'markerfacecolor',cmap(iclr,:),'markeredgecolor',[.5 .5 .5]);
end
clim([vumin vumax]);
colormap(cmap);
colorbar('vert');
xlabel('Longitude');
ylabel('Latitude');
title('Speckle Noise');

axis(bounds);
set(gca,'dataaspectratio',[1 median(cosd(lat)) 1]);


%%

disp('GPS Imaging...');

numpts = 100;
int = min([bounds(2)-bounds(1) bounds(4)-bounds(3)])/numpts;
[LON,LAT]=meshgrid(bounds(1):int:bounds(2),bounds(3):int:bounds(4));

% Here is the part where the mask "INNA" is made.  It is a matrix the same
% size as LON and LAT, but only has 0 or 1.  Areas with 0 are not imaged so
% can save on computation time. 
if usemask == 1
    INNACA=inpolygon(LON,LAT,xCA(1:end-1),yCA(1:end-1));
    INNANV=inpolygon(LON,LAT,xNV(1:end-1),yNV(1:end-1));
    INNA=(INNACA | INNANV);
else
    INNA=ones(size(LON));
end
               
[VU,SVU]=msf2grid(lon,lat,vumsf,su,LON,LAT,INNA,ssf,opt,optsphere,optpar,optrnf);

figure(6);
clf;
pcolor(LON,LAT,VU);
colormap(cmap);
shading interp;
hold on;
plot(xCA,yCA,'k-');
plot(xNV,yNV,'k-');
plot(lon,lat,'k.');
clim([vumin vumax]);
colorbar('vert');
axis(bounds);
xlabel('Longitude');
ylabel('Latitude');
set(gca,'dataaspectratio',[1 median(cosd(lat)) 1]);


figure(7);
clf;
pcolor(LON,LAT,SmoothScene(VU,1));
colormap(cmap);
shading interp;
hold on;
plot(xCA,yCA,'k-');
plot(xNV,yNV,'k-');
plot(lon,lat,'k.');
clim([vumin vumax]);
colorbar('vert');
axis(bounds);
xlabel('Longitude');
ylabel('Latitude');
set(gca,'dataaspectratio',[1 median(cosd(lat)) 1]);

toc

%%  Iterate to reduce artifacts

if optiterate == 1

    disp('Alternative: Iterating to reduce artifacts')

    nits = 10; % or 100 or 1000?
    n=length(vu);
    ns = size(ssf,1);
    s_ssf = 0.2*ones(ns,1);  % 0.2 is assume here. This was not determined from the data and is negotiable.

    [ni,mi]=size(LON);
    VUcube = nan(ni,mi,nits);

    % In each iteration resamples the velocity field within as assumed Gaussian
    % distribute whose standard deviation is defined by the velocity uncertainties.
    % Also resamples the ssf assuming that it is uncertain by the value given
    % by s_ssf above.  In this case its 20%. The velocity field is then median
    % spatiall filtered and imaged in each iteration.  The final field is the
    % median value at each pixel from all the iterations.

    for i=1:nits

        disp(['   Iteration: ' num2str(i)]);

        vui = vu + su.*randn(n,1);

        ssfi = ssf;
        ssfi(:,2) = ssf(:,2) + s_ssf.*randn(ns,1);
        ssfi(ssfi(:,2)>1,2)=1;
        ssfi(ssfi(:,2)<=0,2)=0;

        vumsfi=msf(lon,lat,vui,su,ssfi,opt,optsphere,optpar,optrnf);

        [VUI,~]=msf2grid(lon,lat,vumsfi,su,LON,LAT,INNA,ssfi,opt,optsphere,optpar,optrnf);

        VUcube(:,:,i)=VUI;

    end

    VUmed = median(VUcube,3);

    figure(8);
    clf;
    pcolor(LON,LAT,VUmed);
    colormap(cmap);
    shading interp;
    hold on;
    plot(xCA,yCA,'k-');
    plot(xNV,yNV,'k-');
    plot(lon,lat,'k.');
    clim([vumin vumax]);
    colorbar('vert');
    axis(bounds);
    xlabel('Longitude');
    ylabel('Latitude');
    set(gca,'dataaspectratio',[1 median(cosd(lat)) 1]);

end

%%

save('AfterGPSImaging');

toc