Ln_pre_peak = [12;13;15;17;25;26;30;37;40;41];
Ln_post_peak = [13;14;16;18;26;27;31;38;41;42];

Ln_Anom = zeros(182,482,12,10);
for m = 1:10
    Ln_Anom(:,:,9:12,m) = meso_anom(:,:,9:12,Ln_pre_peak(m));
    Ln_Anom(:,:,1:8,m) = meso_anom(:,:,1:8,Ln_post_peak(m));
end
    
final_Anom = zeros(182,482,12);
for i = 1:182
    for j = 1:482
        for m = 1:12
            final_Anom(i,j,m) = nansum(Ln_Anom(i,j,m,:))/10;
        end
    end
end

figure(1)
suptitle('La Niña anomalies of 0-100m Mesozooplankton')
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
suptitle('La Niña anomalies of 0-100m Mesozooplankton(Post-peak year)')
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
 