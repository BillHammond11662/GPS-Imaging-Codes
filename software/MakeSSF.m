function ssf=MakeSSF(lon,lat,vu,su,optpar,dvmax)


N=length(vu);
dvu = nan(N,N);
R=dvu;
I=nan(N,N);
J=nan(N,N);

disp('Finding distances and dvu of many pairs...');

if optpar
    parfor i=1:N
        if (mod(i,100)==0)
            disp(['i=' sprintf('%.0f',i) ' of ' sprintf('%.0f',N)]);
        end
        for j=1:N
            if i>j
                I(i,j)=i;
                J(i,j)=j;
                dvu(i,j)=vu(i)-vu(j);
                [R(i,j),~] = baz(lat(i),lon(i),lat(j),lon(j));
            end
        end
    end
else
    for i=1:N
        if (mod(i,100)==0)
            disp(['i=' sprintf('%.0f',i) ' of ' sprintf('%.0f',N)]);
        end
        for j=1:N
            if i>j
                I(i,j)=i;
                J(i,j)=j;
                dvu(i,j)=vu(i)-vu(j);
                [R(i,j),~] = baz(lat(i),lon(i),lat(j),lon(j));
            end
        end
    end
end

%% 

% vumsf = msf(lon,lat,vu,su,[],2,0,0,0);

q = find(~isnan(dvu) & abs(dvu)<dvmax);

% These intervals can be adjusted to the data
%xlogbin = -2:.25:1.25;
%xlogbin = [-2 -1:.25:1.25];
xlogbin = [-2 -1:.25:1.25];

xmed=nan(length(xlogbin)-1,1);
vup50=nan(length(xlogbin)-1,1);

Rq=R(q);
dvuq=dvu(q);

vup50(1) = median(su);
xmed(1) = 10^((xlogbin(2)+xlogbin(1))/2);

for i=2:(length(xmed))
    
    x1=10^xlogbin(i);
    x2=10^xlogbin(i+1);
    xmed(i) = 10^((xlogbin(i+1)+xlogbin(i))/2);
    j=find(Rq>=x1 & Rq<=x2);
    
    vup50(i) = max([vup50(1:(i-1));median(abs( dvuq(j)-median(dvuq(j))))]);
        
end

% equations 1 and 2 of Hammond et al., 2016
ssf = 1./vup50;
ssf = ssf/max(ssf);
ssf(end)=0;

% smooths it down to 0 at end (has little effect since only connected points
% in triangulation are weighted by ssf)
ssf(end-1)=mean(ssf([end-2 end]));
ssf(end-2)=mean(ssf([end-3 end-1]));

% extend domain to 0 and 180 just to make sure no points fall outside
xmed = [0;xmed;180];
ssf = [1;ssf;0];
ssf = [xmed ssf];

   
%%

figure(2);
clf;

subplot(311);
histogram(log10(R(q)),xlogbin);
set(gca,'xlim',[xlogbin(1) xlogbin(end)]);
xt=get(gca,'xticklabel');
xtl=cell(size(xt));
for  i=1:length(xt)
   xtl{i}=['10' '^{' char(xt{i}) '}'];
end
hold on;
grid on;
set(gca,'xticklabel',xtl');
title('Histogram of Baseline Distances');
xlabel('Distance Between Station Pairs (degrees)');

subplot(312);
semilogx(R(q),abs(dvu(q)),'k.');
grid on;
xlabel('Distance Between Station Pairs (degrees)');
ylabel('Vu difference (mm/yr)');
set(gca,'xlim',[10^xlogbin(1) 10^xlogbin(end)]);
set(gca,'ylim',[.1 dvmax]);
hold on;

subplot(313)
hs=semilogx(ssf(:,1),ssf(:,2),'k-');
set(hs,'Linewidth',2)
grid on;
set(gca,'xlim',[10^xlogbin(1) 10^xlogbin(end)]);
set(gca,'ylim',[0 1]);
xlabel('Distance Between Station Pairs (degrees)');
title('Spatial Structure Function');
