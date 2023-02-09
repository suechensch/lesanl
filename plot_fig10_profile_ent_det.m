
clear all
close all
%addpath '/users/chen/Matlab/functions';
addpath '/p/home/chens/Matlab/functions';

exps={'camp2xSep27','sstSep27','camp2xSep28-t2','sstSep28'};
subn='new8';
%exps={'sstSep27'};
kpath=['/p/work1/chens/les/data/' exps{1} '/'];
dtg='2019092721';
nest=1;
%outpath=['/chen4/cloud/matfiles/'];
yy=str2num(dtg(1:4));
mm=str2num(dtg(5:6));
dd=str2num(dtg(7:8));
hh=str2num(dtg(9:10));
%opath=['/p/cwfs/chens/' exps{1} '/matfiles/'];
[data_grid]=read_datahd_ff(dtg, kpath);
ma=data_grid.nest.nx;
na=data_grid.nest.ny;
kka=data_grid.nz;
mo=num2str(ma,'%04d');
no=num2str(na,'%04d');
Z=flipud(data_grid.sigm);
k2=20; % zk mean vertical array
mylabel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
%exps={'sstSep28','camp2xSep28-t2'};
nexp=length(exps);
for np=1:nexp
%for np=[1]

%outpath=['/p/cwfs/chens/cld_les/' exps{np} '/matfiles/'];
  outpath=['/p/cwfs/chens/' exps{np} '/matfiles/new7/'];
  sfile=[outpath exps{np} '_cldtb_ent_det_' subn '.mat'];
  load(sfile)

  [nt,nk]=size(h_mentr_new);

%  hr_entr(:,:,np)=h_mentr;
%  hr_detr(:,:,np)=h_mdetr;
%  hr_entr(:,:,np)=h_mentr_new;
%  hr_detr(:,:,np)=h_mdetr_new;
  ick=0;
iavg=0;
if iavg==1
% average to 100 m
% average to 40 m
%  for k=1:10:kka
  for k=1:4:kka
    en_std(k,np)=nanstd(h_smentr(k,:));
    de_std(k,np)=nanstd(h_smdetr(k,:));
    ik1=(k-1)+1;
    if ik1 > kka
     break;
    else
    ik2=ik1+9; 
    ick=ick+1;
    zk(ick)=Z(ik1);    
    hr_entr(ick,:,np)=nanmean(h_smentr(ik1:ik2,:),1);
    hr_detr(ick,:,np)=nanmean(h_smdetr(ik1:ik2,:),1);
    end
  end

else
  hr_entr(:,:,np)=h_smentr;
  hr_detr(:,:,np)=h_smdetr;
end %iavg
 
% mean dmass/dz
   hr_dmdz(:,:,np)=dmassdz2;
   hr_mass(:,:,np)=mass2;

end

zk=Z;
k2=80;
fc=1000; % cover to km-1
mcurve=fc*nanmean(hr_entr(1:k2,:,1),2);
mcurve2=fc*nanmean(hr_entr(1:k2,:,2),2);
mcurve3=fc*nanmean(hr_entr(1:k2,:,3),2);
mcurve4=fc*nanmean(hr_entr(1:k2,:,4),2);
%curve1 = mcurve-sqrt(fc*en_std(1:k2,1));
%curve2 = mcurve+sqrt(fc*en_std(1:k2,1));
%xy2 = [zk(1:k2)/fc, fliplr(zk(1:k2)/fc)];
%inBetween = [curve1, curve2];
figure;
  subplot(1,2,1)
hold all
%TT=tiledlayout(2,1);
%TT.Padding = 'compact';
%TT.TileSpacing = 'compact';

plot(fc*nanmean(hr_entr(1:k2,:,1),2),zk(1:k2)/fc,'b','linewidth', 2)
plot(fc*nanmean(hr_entr(1:k2,:,2),2),zk(1:k2)/fc,'r','linewidth', 2)
plot(fc*nanmean(hr_entr(1:k2,:,3),2),zk(1:k2)/fc,'k','linewidth', 2)
plot(fc*nanmean(hr_entr(1:k2,:,4),2),zk(1:k2)/fc,'m','linewidth', 2)
%plot(mcurve,zk(1:k2)/fc,'b','linewidth', 2)
%plot(mcurve2,zk(1:k2)/fc,'r','linewidth', 2)
%plot(mcurve3,zk(1:k2)/fc,'k','linewidth', 2)
%plot(mcurve4,zk(1:k2)/fc,'m','linewidth', 2)
%xxyy2=fill(curve1,curve2, [43 54 255]/255); alpha(xxyy2,0.5);
ylabel('Height (km)')
xlabel('Entrainment Rate (km^-^1)')
  title(mylabel{1})
  ax = gca;
  ax.TitleHorizontalAlignment = 'left';
box on;
legend('L1','L1sst','L2','L2sst','location','southeast')
%legend('L1','L1sst','L2','L2sst','location','northeast')
%xlim([0 3e-2])
ylim([0.6 1.8])
oname=[outpath  'mean_entr_' subn '_profile.png'];
%print ('-dpng',oname)

%figure;
%  ax = nexttile;
  subplot(1,2,2)
hold all
mcurve=fc*nanmean(hr_detr(1:k2,:,1),2);
mcurve2=fc*nanmean(hr_detr(1:k2,:,2),2);
mcurve3=fc*nanmean(hr_detr(1:k2,:,3),2);
mcurve4=fc*nanmean(hr_detr(1:k2,:,4),2);
plot(fc*nanmean(hr_detr(1:k2,:,1),2),zk(1:k2)/fc,'b','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,2),2),zk(1:k2)/fc,'r','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,3),2),zk(1:k2)/fc,'k','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,4),2),zk(1:k2)/fc,'m','linewidth', 2)
%plot(mcurve,zk(1:k2)/fc,'b','linewidth', 2)
%plot(mcurve2,zk(1:k2)/fc,'r','linewidth', 2)
%plot(mcurve3,zk(1:k2)/fc,'k','linewidth', 2)
%plot(mcurve4,zk(1:k2)/fc,'m','linewidth', 2)
%xxyy2=fill(curve1,curve2, [43 54 255]/255); alpha(xxyy2,0.5);
ylabel('Height (km)')
xlabel('Detrainment Rate (km^-^1)')
ylim([0.6 1.8])
  title(mylabel{2})
%  ax = gca;
  ax.TitleHorizontalAlignment = 'left';
box on;
ylim([0.6 1.8])
%legend('L1','L1sst','L2','location','southeast')
%legend('L1','L1sst','L2','L2sst','location','northeast')
%xlim([0 3e-2])
oname=[outpath  'mean_detr_' subn '_profile.png'];
print ('-dpng',oname)

figure;
hold all;

colormap(jet)
set(gca,'fontsize',15);
plot(fc*nanmean(hr_entr(1:k2,:,1),2),zk(1:k2)/fc,'b','linewidth', 2)
plot(fc*nanmean(hr_entr(1:k2,:,2),2),zk(1:k2)/fc,'r','linewidth', 2)
plot(fc*nanmean(hr_entr(1:k2,:,3),2),zk(1:k2)/fc,'k','linewidth', 2)
plot(fc*nanmean(hr_entr(1:k2,:,4),2),zk(1:k2)/fc,'m','linewidth', 2)

%plot(fc*nanmean(hr_entr(1:k2,1:10,1),2),zk(1:k2)/fc,'color',[255 204 204]/255,'linewidth', 1)
%plot(fc*nanmean(hr_entr(1:k2,1:10,2),2),zk(1:k2)/fc,'color',[204 255 255]/255,'linewidth', 1)
%plot(fc*nanmean(hr_entr(1:k2,1:10,3),2),zk(1:k2)/fc,'color',[204 255 255]/255,'linewidth', 1)

plot(fc*nanmean(hr_detr(1:k2,:,1),2),zk(1:k2)/fc,'-.b','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,2),2),zk(1:k2)/fc,'-.r','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,3),2),zk(1:k2)/fc,'-.k','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,4),2),zk(1:k2)/fc,'-.m','linewidth', 2)

%plot(fc*nanmean(hr_detr(1:k2,1:10,1),2),zk(1:k2)/fc,'color',[255 204 204]/255,'linewidth',1,'linestyle','-.')
%plot(fc*nanmean(hr_detr(1:k2,1:10,2),2),zk(1:k2)/fc,'color',[204 255 255]/255,'linewidth', 1,'linestyle','-.')
%plot(fc*nanmean(hr_detr(1:k2,1:10,3),2),zk(1:k2)/fc,'color',[204 255 255]/255,'linewidth', 1,'linestyle','-.')

ylabel('Height (km)')
xlabel('Entrainment/Detrainment Rate (km^-^1)')
box on;
ylim([0.5 1.8])
legend('L1','L1sst','L2','L2sst','location','northeast')
%legend('L1','L1sst','L2','L2sst','location','southeast')
%xlim([0 3e-2])
%ylim([0.6 2])
oname=[outpath  'mean_entr-detr_' subn '_profile.png'];
print ('-dpng',oname)

k3=70;

figure;
hold all;
colormap(jet)
set(gca,'fontsize',15);
plot(nanmean(massdz2(1:k3,:),2),Z(1:k3)/fc,'b','linewidth', 2)
ylim([0.6 1.8])
title('1/mass * dmass/dz');
ylabel('Height (km)')
oname=[outpath  'massdz_' subn '_profile.png'];
print ('-dpng',oname)

figure;
  subplot(1,2,1)
hold all;
colormap(jet)
set(gca,'fontsize',15);
plot(nanmean(hr_mass(1:k3,:,1),2),Z(1:k3)/fc,'b','linewidth', 2)
plot(nanmean(hr_mass(1:k3,:,2),2),Z(1:k3)/fc,'r','linewidth', 2)
plot(nanmean(hr_mass(1:k3,:,3),2),Z(1:k3)/fc,'k','linewidth', 2)
plot(nanmean(hr_mass(1:k3,:,4),2),Z(1:k3)/fc,'m','linewidth', 2)
legend('L1','L1sst','L2','L2sst','location','northeast')
ylim([0.6 1.8])
xlabel('mass (kg m^-^2 s^-^1)');
ylabel('Height (km)')
  title('(a)')
  ax = gca;
  ax.TitleHorizontalAlignment = 'left';
%oname=[outpath  'mass_' subn '_profile.png'];
%print ('-dpng',oname)

%figure;
  subplot(1,2,2)
hold all;
colormap(jet)
set(gca,'fontsize',15);
%plot(nanmean(dmassdz(1:k3,:),2),Z(1:k3)/fc,'b','linewidth', 2)
plot(nanmean(hr_dmdz(1:k3,:,1),2),Z(1:k3)/fc,'b','linewidth', 2)
plot(nanmean(hr_dmdz(1:k3,:,2),2),Z(1:k3)/fc,'r','linewidth', 2)
plot(nanmean(hr_dmdz(1:k3,:,3),2),Z(1:k3)/fc,'k','linewidth', 2)
plot(nanmean(hr_dmdz(1:k3,:,4),2),Z(1:k3)/fc,'m','linewidth', 2)
ylim([0.6 1.8])
legend('L1','L1sst','L2','L2sst','location','northwest')
xlabel('dmass/dz (kg m^-^3 s^-^1)');
ylabel('Height (km)')
  title('(b)')
  ax = gca;
  ax.TitleHorizontalAlignment = 'left';
oname=[outpath  'dmassdz_' subn '_profile.png'];
print ('-dpng',oname)
figname=[outpath  'dmassdz_' subn '_profile.fig'];
savefig(figname)

ifplot=0;

if ifplot==1

figure;
hold all;

colormap(jet)
set(gca,'fontsize',15);
plot(fc*nanmean(hr_detr(1:k2,:,1),2),Z(1:k2)/fc,'b','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,2),2),Z(1:k2)/fc,'r','linewidth', 2)
plot(fc*nanmean(hr_detr(1:k2,:,3),2),Z(1:k2)/fc,'k','linewidth', 2)
ylim([0.6 2])
ylabel('Height (km)')
xlabel('Detrainment rate (km^-^1)')
box on;
%legend('L1','L1sst','L2','location','northeast')
legend('L1','L1sst','L2','L2sst','location','southeast')
%xlim([0 3e-2])
oname=[outpath  'mean_detr_new5_profile.png'];
print ('-dpng',oname)

end
