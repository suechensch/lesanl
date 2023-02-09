
clear
close all

dpath ='../data/';
fpath ='../figures/';
sname=[dpath 'tke_cld_all.mat']
load(sname)
exps={'L1','L1sst','L2','L2sst'};
np=length(exps);

figure;
set(gca,'fontsize',15)
tmp=squeeze(nanmean(tke_cloud,2));
tmp2=squeeze(nanmean(tke_clear,2));
tmp3=squeeze(nanmean(tke_all,2));
tmp4=squeeze(nanmean(tke_sc,2));
mycolor={'b','m',[0.4940 0.1840 0.5560],[0.8500 0.3250 0.0980]};


  set(gca,'fontsize',15)
  hold all
  for i=1:np
    plot(tmp4(:,i),Z,'color',mycolor{i},'linewidth',2)
    x=[6,30];
    y=[8,30];
  end
  ylim([0 1500])
  legend('L1','L1sst','L2','L2sst','location','southeast')
% create smaller axes in top right, and plot on it
  axes('Position',[.7 .7 .2 .2])
  box on
  hold all
  for i=1:np
  plot(tmp4(:,i),Z,'color',mycolor{i},'linewidth',2)
  x=[6,30];
  y=[8,30];
  end
  ylim([0 50])

xlabel('TKE (m^2/s^2)')
ylabel('Height (m)')
fname=[fpath 'fig7_sc_tke_cloud_profile.png']
print ('-dpng',fname)

