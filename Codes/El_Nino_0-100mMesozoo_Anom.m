cd('D:\19580112-20070112')
% Read longitude, latitude, and temperature from files
long = ncread('CCS_7k_0-360_fred_grd.nc1','lon_rho');
lat = ncread('CCS_7k_0-360_fred_grd.nc1','lat_rho');s

% z_r - the 3D matrix of depth at every grid point
N = 50;
s_rho = ncread('CCS_7k_0-360_fred_grd.nc1','s_rho');
hc = ncread('CCS_7k_0-360_fred_grd.nc1','hc');
h = ncread('CCS_7k_0-360_fred_grd.nc1','h');
Cs_r = ncread('CCS_7k_0-360_fred_grd.nc1','Cs_r');
 
zeta = ncread('1959-01-12.nc','zeta'); % free surface elevation (zet ; 
zeta=squeeze(zeta(:,:,10));

for  k = 1:N;
    z0 = (hc.*s_rho(k) + h * Cs_r(k))./(hc + h);
    z_r(:,:,k) = zeta + (zeta + h).*z0; 
end

%depth_ind - a column vector, each value is the depth layer at which
%actual depth is closest to 100m below the sea surface

depth_ind = zeros(size(long));
for i = 1:182
    for j = 1:482
        x = 1;
        for l = 2:50
            %compare the differences between the actual depths and -100m
            if abs(z_r(i,j,l)+100) < abs(z_r(i,j,x)+100)
                x = l;
            end
        end
        depth_ind(i,j) = x;
    end
end

%Average El Nino Variable Anomalies over 100m 
month_avg= zeros(182,482,588);
projectdir = 'D:\19580112-20070112';
dinfo = dir( fullfile(projectdir, '*.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );

for K = 1 : num_files
    this_file = filenames{K};
    %read mesozoo value from each file
    mesozoo = ncread(this_file, 'mesozooplankton');
    for i = 1:182
        for j = 1:482
            for m = 1:12
                grid = squeeze(mesozoo(i,j,:,m));
                %depth - actual depths that are above -100m of one grid cell
                depth = squeeze(z_r(i,j,depth_ind(i,j):50));
                total_depth = depth(length(depth))-depth(1);
                %grid_meso - mesozoo at depths above -100m of one grid cell
                grid_meso = grid((depth_ind(i,j)+2):50,1);
                %sum - the sum of mesozoo*depth interval
                sum = 0;
                for l = 1:length(grid_meso)
                    sum = sum + ((depth(l+2,1)-depth(l,1))/2)*grid_meso(l,1);
                end
                if abs(total_depth-100)>10
                   month_avg(i,j,:) = NaN;
                end
                month_avg(i,j,12*(K-1)+m) =  sum/total_depth;
            end
        end
    end
end

%meso_mean - monthly mean mesozoo over 100m
meso_mean = zeros(182,482,588);
for i = 1:12
    sum = zeros(size(long));
    for j = 1:49
        sum = sum + squeeze(month_avg(:,:,12*(j-1)+i));
    end
    meso_mean(:,:,i) = sum/49;
    for m = 2:49
        meso_mean(:,:,12*(m-1)+i) = meso_mean(:,:,i);
    end
end

meso_anom = zeros(182,482,12,49);
for i = 1:182
    for j = 1:482
        anom = squeeze(month_avg(i,j,:)- meso_mean(i,j,:));
        anom = lanczosfilter(anom,1,1/121,[],'high');
        for m = 1:49
            meso_anom(i,j,:,m) = anom((12*(m-1)+1):(12*m));
        end
    end
end

cd('D:\19580112-20070112\var_Anoms_CCS')
load('mesozoo_0-100m_anom.mat')
El_pre_peak = [5;7;10;14;24;28;29;33;36;39;44];
El_post_peak = [6;8;11;15;25;29;30;34;37;40;45];

El_Anom = zeros(182,482,12,11);
for m = 1:11
    El_Anom(:,:,9:12,m) = meso_anom(:,:,9:12,El_pre_peak(m));
    El_Anom(:,:,1:8,m) = meso_anom(:,:,1:8,El_post_peak(m));
end


%average monthly mean mesozoo anomalies in 11 El Nino years
final_Anom = zeros(182,482,12);
for i = 1:182
    for j = 1:482
        for m = 1:12
            final_Anom(i,j,m) = nansum(El_Anom(i,j,m,:))/11;
        end
    end
end

figure(1)
suptitle('El Niño anomalies of 0-100m Mesozooplankton')
subplot(2,4,1)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,9)));
xlabel('Sep');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold on


subplot(2,4,2)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,10)));
xlabel('Oct');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');;
m_grid('box','on','linestyle','none','xtick',4);
hold on
    
subplot(2,4,3)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,11)));
xlabel('Nov');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');;
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,4)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,12)));
xlabel('Dec');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');;
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,5)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,1)));
xlabel('Jan');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');;
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,6)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,2)));
xlabel('Feb');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');;
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,7)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,3)));
xlabel('Mar');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map');
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,8)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,4)));
xlabel('Apr');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold off
h  = colorbar('Position',[0.92 0.2 0.0126370011898221 0.672734475483289]);
h.FontSize = 12;
h.Label.FontSize = 10;

h.Label.String = 'mmol N_2 m^-^3' ;


figure(2)
suptitle('El Niño anomalies of 0-100m Mesozooplankton(Post-peak year)')
subplot(2,4,1)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,1)));
xlabel('Jan');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');;
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,2)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,2)));
xlabel('Feb');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');;
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,3)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,3)));
xlabel('Mar');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map');
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,4)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,4)));
xlabel('Apr');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,5)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,5)));
xlabel('May');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold on
    
subplot(2,4,6)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,6)));
xlabel('June');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,7)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,7)));
xlabel('July');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold on

subplot(2,4,8)
cd('D:\m_map1.4\m_map')
m_proj('lambert','long',[-142 -110],'lat',[18 52]);
my_map = m_pcolor(long-360,lat,squeeze(final_Anom(:,:,8)));
xlabel('Aug');
set(gca,'FontSize',9);
cd('D:\kakearney-cptcmap-pkg-845bf83\kakearney-cptcmap-pkg-845bf83\cptcmap')
cptcmap('temp_19lev.cpt');
caxis([-0.03 0.03])
cd('D:\m_map1.4\m_map')
m_coast('patch',[0.4 0.4 0.3],'edgecolor','black');
m_grid('box','on','linestyle','none','xtick',4);
hold off

h  = colorbar('Position',[0.92 0.2 0.0126370011898221 0.672734475483289]);
h.FontSize = 12;
h.Label.FontSize = 10;

h.Label.String = 'mmol N_2 m^-^3' ;
 