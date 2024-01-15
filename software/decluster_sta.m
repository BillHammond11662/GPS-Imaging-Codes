function iclusters=decluster_sta(lat,lon,dcol,optverb)
% iclusters=decluster_sta(lat,lon,dcol,optverb)
%
% dcol is in km

if any(isnan(lon)) || any(isnan(lat))
   error('No Nans allowed in lat or lon'); 
end

inotdone = [1:length(lat)]';
idone=[];

done=0;
ccnt = 0;
while ~done
    ccnt = ccnt+1;
    
    [dist,~] = distance(lat,lon,lat(1),lon(1));
    jnear=find(dist<km2deg(dcol));
    
    iclusters{ccnt}=inotdone(jnear);
    idone=[idone;inotdone(jnear)];
    inotdone(jnear)=[];
    lat(jnear)=[];
    lon(jnear)=[];
    
    if optverb==1
        disp([num2str(ccnt) ' : ' num2str(length(jnear)) ': not done ' ...
            num2str(length(inotdone)) ': done ' num2str(length(idone)) ': ' num2str(length(lat))]);
    end
    %     for p=1:length(jnear)
    %         disp(['   :' num2str(jnear(p))]);
    %     end
    
    if isempty(inotdone)
        done=1;
    end
end

