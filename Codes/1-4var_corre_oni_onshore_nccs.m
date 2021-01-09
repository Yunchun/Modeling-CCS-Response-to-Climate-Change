cd('D:\19580112-20070112')
%read ONI from file
oni = textread('oni.txt');
%remove the first column(year index)
oni = (oni(:,2:13))';
%reshape matrix into a vector
oni= oni(:);

%Remove decadal variations of anomalies using Lanczos Filter
oni = lanczosfilter(oni,1,1/121,[],'high');

%oni_avg - oni averaged over 3 months(jan-apr)
oni_avg = zeros(49,1)
for i = 0:(length(oni_avg)-1)
    %get onis from jan to apr
    jan_apr = oni((12*i+1):(12*i+4));
    oni_avg(i+1,1) = nanmean(jan_apr);
end

long = ncread('CCS_7k_0-360_fred_grd.nc','lon_rho');
lat = ncread('CCS_7k_0-360_fred_grd.nc','lat_rho');
%lon - the indices of longitude near onshore northern CCS
%lat - the indices of latitude near onshore northern CCS
lon = find(233 <= long & long <= 237);
la = find(38 <= lat & lat <= 42);

%index - indices that are shared in lon and lat
index = zeros(1,1);
for l = 1:length(lon)
    for p = 1:length(la)
        if lon(l,1)==la(p,1)
            index1 = la(p,1);
            index = cat(1,index,index1);
        end
    end
end
index = index(2:length(index),1);

%row - row indices of lon and lat of onshore northern CCS
%col - column indices of lon and lat of onshore northern CCS
col = zeros(length(index),1);
row = zeros(length(index),1);
for i = 1:length(index)
    if rem(index(i),182) == 0
        col(i,1)= fix(index(i)/182);
        row(i,1) = 182;
    else
        col(i,1)= fix(index(i)/182)+1;
        row(i,1) = rem(index(i),182);
    end
end

%remove the repeating elements in row and col
row = unique(row);
col = unique(col);

% z_r - the 3D matrix of depth at every grid point
N = 50;
s_rho = ncread('CCS_7k_0-360_fred_grd.nc','s_rho');
hc = ncread('CCS_7k_0-360_fred_grd.nc','hc');
h = ncread('CCS_7k_0-360_fred_grd.nc','h');
Cs_r = ncread('CCS_7k_0-360_fred_grd.nc','Cs_r');
 
zeta = ncread('1959-01-12.nc','zeta'); % free surface elevation (zet ; 
zeta=squeeze(zeta(:,:,10));

for  k = 1:N;
    z0 = (hc.*s_rho(k) + h * Cs_r(k))./(hc + h);
    z_r(k,:,:) = zeta + (zeta + h).*z0; 
end

%depth_table - the table where each element is the actual depth at grid
%points near onshore northern CCS
depth_table = zeros(50,length(row),length(col));
for i = 1:length(row)
    for j = 1:length(col)
        depth_table(:,i,j) = squeeze(z_r(:,row(i,1),col(j,1)));
    end
end

%depth_ind - a column vector, each value is the depth layer at which
%actual depth is closest to 100m below the sea surface

depth_ind = zeros(length(row),length(col));
for i = 1:length(row)
    for j = 1:length(col)
        x = 1;
        for l = 2:50
            %compare the differences between the actual depths and -100m
            if abs(depth_table(l,i,j)+100) < abs(depth_table(x,i,j)+100)
                x = l;
            end
        end
        depth_ind(i,j) = x;
    end
end


%month_avg - col:each grid point, row: each month in each year
%each value is variable averaged over onshore northern ccs area 
%and over 100m of each grid point
month_avg= zeros(49,length(row),length(col));
projectdir = 'D:\19580112-20070112';
dinfo = dir( fullfile(projectdir, '*12.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

%calculate Pzooplankton averaged over 100m from 1-4
for K = 1 : num_files
    this_file = filenames{K};
    %read pzoo value from each file
    pzoo = ncread(this_file, 'Pzooplankton');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(pzoo(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_temp - pzoo at depths above -100m of one grid cell
                grid_pzoo = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_pzoo),1);
                %sum - the sum of pzoo*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_pzoo(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of pzoo anomalies with jan-apr's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end
    
new_long = long(row,col);
new_lat = lat(row,col);

figure(1)
h=suptitle({'correlation of anomalies with ONI over Northern CCS','lag=3 months(Jan-Apr)'});
set(h,'FontSize',10,'FontWeight','normal')
subplot(2,3,1)
%set up the projection of m_map
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m Pzoo');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar;
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on


%Calculate Microzoo values averaged over 100m from 1-4
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read microzoo value from each file
    micro = ncread(this_file, 'microzooplankton');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(micro(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_micro - microzoo at depths above -100m of one grid cell
                grid_micro = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_micro),1);
                %sum - the sum of microzoo*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_micro(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's microzoo anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,2)
%set up the projection of m_map
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m Microzoo');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar;
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on


%Calculate Mesozoo values averaged over 100m from 1-4
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read mesozoo value averaged over 100m from 1-4
    meso = ncread(this_file, 'mesozooplankton');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(meso(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_meso - mesozoo at depths above -100m of one grid cell
                grid_meso = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_meso),1);
                %sum - the sum of mesozoo*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_meso(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's meso anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,3)
%set up the projection of m_map
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m Mesozoo');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on



%Calculate nanophyto values averaged over 100m from 1-4
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read nanophyto value from each file
    nano = ncread(this_file, 'nanophytoplankton');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(nano(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_nano - nanophyto at depths above -100m of one grid cell
                grid_nano = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_nano),1);
                %sum - the sum of nanophyto*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_nano(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's nanophyto anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,4)
%set up the projection of m_map
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m Nanophy');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on


%Calculate diatom values averaged over 100m from 1-4
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read diatom value from each file
    diatom = ncread(this_file, 'diatom');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(diatom(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_dia - diatom at depths above -100m of one grid cell
                grid_dia = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_dia),1);
                %sum - the sum of diatom*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_dia(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's diatom anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,6)
%set up the projection of m_map
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m Diatom');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold off




%Calculate temperature values averaged over 100m from 1-4
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read the temperature value from each file
    temp = ncread(this_file, 'temp');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(temp(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_temp - temperature at depths above -100m of one grid cell
                grid_temp = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_temp),1);
                %sum - the sum of T*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_temp(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's st anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

figure(2)
h = suptitle({'correlation of anomalies with ONI over Northern CCS','lag=3 months(Jan-Apr)'});
set(h,'FontSize',12,'FontWeight','normal')
subplot(2,3,1)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m Sea Temp');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on


%Calculate salinity averaged over 100m from 1-4
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read salt values from each file
    salt = ncread(this_file, 'salt');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(salt(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_salt - salinity at depths above -100m of one grid cell
                grid_salt = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_salt),1);
                %sum - the sum of salinity*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_salt(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's salinity anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,4)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m Salt');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on


%Calculate u velocity averaged over 100m from 1-4 
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read u velocity value from each file
    u = ncread(this_file, 'u');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(u(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_u - u velocity at depths above -100m of one grid cell
                grid_u = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_u),1);
                %sum - the sum of u velocity*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_u(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's u velocity anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,2)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m u velocity');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on


%Calculate v velocity values averaged over 100m from 1-4
month_avg= zeros(49,length(row),length(col));

cd('D:\19580112-20070112')
for K = 1 : num_files
    this_file = filenames{K};
    %read v velocity value from each file
    v = ncread(this_file, 'v');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(v(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(i,j):50,i,j);
                %grid_v - v velocity at depths above -100m of one grid cell
                grid_v = grid((depth_ind(i,j)+2):50,1);
                %interval - dz
                interval = zeros(length(grid_v),1);
                %sum - the sum of v velocity*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_v(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(size(depth_ind));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of 1-4's v velocity anomalies with 1-4's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,5)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('0-100m v velocity');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar 
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on


%Average NO3 over 25-100m
depth_ind = zeros(2,length(row),length(col));
for i = 1:length(row)
    for j = 1:length(col)
        x = 1;
        y = 1;
        for l = 2:50
            %compare the differences between the actual depths and -100m
            if abs(depth_table(l,i,j)+100) < abs(depth_table(x,i,j)+100)
                x = l;
            %compare the differences between the actual depths and -100m
            elseif abs(depth_table(l,i,j)+25) < abs(depth_table(y,i,j)+25)
                y = l;
            end   
        end
        depth_ind(1,i,j) = x;
        depth_ind(2,i,j) = y;
    end
end


month_avg= zeros(49,length(row),length(col));


%calculate NO3 averaged over 25-100m from 1-4
for K = 1 : num_files
    this_file = filenames{K};
    %read NO3 value from each file
    NO3 = ncread(this_file, 'NO3');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(NO3(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(1,i,j):depth_ind(2,i,j),i,j);
                %grid_NO3 - NO3 at depths above -100m of one grid cell
                grid_NO3 = grid((depth_ind(1,i,j)+2):(depth_ind(2,i,j)),1);
                %interval - dz
                interval = zeros(length(grid_NO3),1);
                %sum - the sum of NO3*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_NO3(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(length(row),length(col));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of NO3 anomalies with jan-apr's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,3)
%set up the projection of m_map
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('25-100m NO3');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar;
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold on



month_avg= zeros(49,length(row),length(col));


%calculate SiOH4 averaged over 25-100m from 1-4
for K = 1 : num_files
    this_file = filenames{K};
    %read SiOH4 value from each file
    SiOH4 = ncread(this_file, 'SiOH4');
    for i = 1:length(row)
        for j = 1:length(col)
            avg = 0;
            for m = 1:4
                grid = squeeze(SiOH4(row(i),col(j),:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = depth_table(depth_ind(1,i,j):depth_ind(2,i,j),i,j);
                %grid_si - SiOH4 at depths above -100m of one grid cell
                grid_si = grid((depth_ind(1,i,j)+2):(depth_ind(2,i,j)),1);
                %interval - dz
                interval = zeros(length(grid_si),1);
                %sum - the sum of SiOH4*depth interval
                sum = 0;
                for l = 1:length(interval)
                    interval(l,1) = (depth(l+2,1)-depth(l,1))/2 ;
                    sum = sum + interval(l,1)*grid_si(l,1);
                end
                avg = avg + sum/nansum(interval);
            end
            month_avg(K,i,j) = nanmean(avg);
        end
    end
end

cor = zeros(length(row),length(col));
for i = 1:length(row)
    for j = 1:length(col)
        mean = nanmean(month_avg(:,i,j));
        anom = month_avg(:,i,j) - mean;
        %cor - the correlation of SiOH4 anomalies with jan-apr's oni
        cor(i,j) = corr(oni_avg,anom);
    end
end

subplot(2,3,6)
%set up the projection of m_map
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-127 -123],'lat',[38 42]);
my_map = m_pcolor(new_long-360,new_lat,cor);
xlabel('25-100m SiOH4');
set(gca,'FontSize',8);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.7 0.7]);
colorbar;
cd('D:\m_map1.4\m_map')
m_coast;
m_grid;
hold off