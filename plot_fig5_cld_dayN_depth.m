
clear all
close all
expn={'L1','L1sst','L2','L2sst'};

ifload=1; % load saved matfile
nbin=10;

if ifload == 1
figure;
colormap(jet)
subplot(2,1,1)
hold all
for np=[1,3]

exp=expn{np};
kpath=['../data/' exp '/'];
outpath=['../figures/'];

lfile=[kpath exp '_cld_depth.mat'];
  load(lfile)
  if np == 1
    h1d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor','b','Linewidth',1);
    h1n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor','m','Linewidth',1);
    h1n.BinWidth=250;
    h1d.BinWidth=250;
  elseif np == 2
    h2d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor','b','Linewidth',1);
    h2n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','coun','FaceColor','none','EdgeColor','m','Linewidth',1);
    h2n.BinWidth=250;
    h2d.BinWidth=250;
  elseif np == 3
    h3d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor',[0.4940 0.1840 0.5560],'Linewidth',1);
    h3n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','countdensity','FaceColor','none','EdgeColor',[0.8500 0.3250 0.0980],'Linewidth',1);
    h3n.BinWidth=250;
    h3d.BinWidth=250;
  elseif np == 4
    h4d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','countdensity','FaceColor','none','EdgeColor',[0.8500 0.3250 0.0980],'Linewidth',1);
    h4n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor',[0.4940 0.1840 0.5560],'Linewidth',1);
    h4n.BinWidth=250;
    h4d.BinWidth=250;
  end

end % nexp
set(gca,'fontsize',15);
xlabel('cloud depth (m)')
ylabel('Count')
if strcmp(exp,'sstSep28')==1
  legend('L2 day','L2 night','L2sst day','L2sst night')
else
  legend('L1 day','L1 night','L1sst day','L1sst night')
end
title('(a)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

xlabel('Cloud depth (m)');


colormap(jet)
subplot(2,1,2)
hold all
for np=[2,4]

exp=expn{np};
kpath=['../data/' exp '/'];
outpath=['../figures/'];

lfile=[kpath exp '_cld_depth.mat'];
  load(lfile)
  if np == 1
    h1d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor','b','Linewidth',1);
    h1n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor','m','Linewidth',1);
    h1n.BinWidth=250;
    h1d.BinWidth=250;
  elseif np == 2
    h2d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor','b','Linewidth',1);
    h2n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','coun','FaceColor','none','EdgeColor','m','Linewidth',1);
    h2n.BinWidth=250;
    h2d.BinWidth=250;
  elseif np == 3
    h3d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor',[0.4940 0.1840 0.5560],'Linewidth',1);
    h3n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','countdensity','FaceColor','none','EdgeColor',[0.8500 0.3250 0.0980],'Linewidth',1);
    h3n.BinWidth=250;
    h3d.BinWidth=250;
  elseif np == 4
    h4d=histogram(cld_depth(index_dayt,:,:),nbin,'Normalization','count','FaceColor','none','EdgeColor',[0.4940 0.1840 0.5560],'Linewidth',1);
    h4n=histogram(cld_depth(index_night,:,:),nbin,'Normalization','countdensity','FaceColor','none','EdgeColor',[0.8500 0.3250 0.0980],'Linewidth',1);
    h4n.BinWidth=250;
    h4d.BinWidth=250;
  end

end % nexp
set(gca,'fontsize',15);
xlabel('cloud depth (m)')
ylabel('Count')
  legend('L2 day','L2 night','L2sst day','L2sst night')
title('(b)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

xlabel('Cloud depth (m)');

fname=[outpath exp '_cld_depth_all_hist.png'];
print ('-dpng',fname)

  
else
GMT2LT=8;
sunrise_lt=6;
sunset_lt=18; % model forecast hour + initial hour + GMT to LT
nest=1;
k1=1;
k2=199;
i1=50;
j1=50;
j2=220;
i2=220;

%% start process
t1=4;
t2=21;  % hour
tint=5/60; % hour
nt=int8((t2-t1)/tint+1);
ifplot=0;   %=1 plot cld base 
cld_cut=0.0001; % cloud mixing ratio threshold to sample the cloudy region 
cld_cut_e=1e-9; % cloud mixing ratio threshold to sample the environment region 



for np=1:length(expn)
ic=0;
exp=expn{np};
kpath=['../data/' exp '/'];
outpath=['../figures/'];
if strcmp(exp,'L2') == 1
   dtg='2019092723';
elseif strcmp(exp,'L1') == 1
   dtg='2019092721';
elseif strcmp(exp,'L1sst') == 1
   dtg='2019092721';
elseif strcmp(exp,'L2sst') == 1
   dtg='2019092723';
end

yy=str2num(dtg(1:4));
mm=str2num(dtg(5:6));
dd=str2num(dtg(7:8));
hh=str2num(dtg(9:10));
[data_grid]=read_datahd_ff(dtg, kpath);
ma=data_grid.nest.nx;
na=data_grid.nest.ny;
kka=data_grid.nz;
mo=num2str(ma,'%04d');
no=num2str(na,'%04d');
Z=flipud(data_grid.sigm);
z2=num2str(data_grid.sigm(1),'%06d');
ztop=data_grid.ztop;

kka=data_grid.nz;

fld(1:kka,1:nt)=0.0;
datev(1:nt)=0;
cld_base(1:na,1:ma)=0.;
cld_top(1:na,1:ma)=0.;

for t=t1:tint:t2
  ic=ic+1;
  [fy,fm,fd,fh,fmin,fs]=datevec(datenum(yy,mm,dd,0+t,0,0));
  if t >=24
    fh=fh+24;
  end
  str2=[num2str(fh,'%04d') num2str(fmin,'%02d') num2str(fs,'%02d') ]
  datev(ic)=datenum(fy,fm,fd,fh,fmin,fs);
  fname=['cldmix_sig_' z2 '_000005_' num2str(nest,'%1d') 'a' mo 'x' no '_' dtg '_' str2 '_fcstfld'];
  if exist([kpath fname]) == 0
   continue
  end

    cld=readff(kpath,fname,kka);
  for k=k1:k2
   indexc=find(cld(i1:i2,j1:j2,k)>=cld_cut_e);
   if isempty(indexc)
    cld_fract(np,ic,k)=0;
   else
    cld_fract(np,ic,k)=length(indexc)/((i2-i1+1)*(j2-j1+1));
   end
  end
   index_e=find(cld<=cld_cut_e);
   for i=1:na
     for j=1:ma
      index=find(cld(i,j,k1:k2)>=cld_cut);
      ilen=length(index);
      if ilen == 0
        cld_base(i,j)=nan;
        ktop(ic,i,j)=nan;
      else
        cld_base(i,j)=squeeze((index(1)));
        cld_top(i,j)=squeeze(Z(index(ilen)));
      end
     end
   end 
  
cld_depth(ic,:,:)=cld_top-cld_base;
cldt_top(ic,:,:)=cld_top;
cldt_base(ic,:,:)=cld_base;

end % end t
local_time=mod([t1:1:t2]+GMT2LT+hh,24);
local_time(find(local_time == 0))= 24;
index_dayt=find(local_time>=sunrise_lt & local_time< sunset_lt);
index_night=find(local_time<sunrise_lt | local_time>= sunset_lt);
index_dn=find(local_time==sunrise_lt | local_time==sunset_lt);

for ii=1:length(local_time)
  local_str{ii}=[num2str(local_time(ii),'%02d') 'LT'];
end

for ii=1:length(index_dayt)
  dayt_str{ii}=[num2str(local_time(index_dayt(ii)),'%02d') 'LT'];
end

for ii=1:length(index_night)
  night_str{ii}=[num2str(local_time(index_night(ii)),'%02d') 'LT'];
end

figure;
colormap(jet)
hold all
set(gca,'fontsize',15);
xlabel('cloud depth (m)')
ylabel('Number of clouds')
h=histogram(cld_depth(index_dayt,:,:))
h2=histogram(cld_depth(index_night,:,:))
fewerbins(h)
fewerbins(h2)
legend('Day','Night')
xlabel('Cloud depth (m)');
fname=[outpath  'fig5_cld_depth_day_night_hist.png'];
print ('-dpng',fname)

lfile=[outpath exp '_cld_depth.mat'];
save(lfile,'cld_depth','index_dayt','index_night','cldt_top','cldt_base','cld_fract')


for i=1:length(expn)
   Z2d(i,:)=Z;
end

figure;
colormap(jet)
hold all
set(gca,'fontsize',15);
xlabel('cloud fraction')
ylabel('Height (m)')

for i=1:length(expn)
  if i==1
   plot(squeeze(nanmean(cld_fract(i,:,k1:k2),2)),Z(k1:k2),'b','linewidth',2)
  else
   plot(squeeze(nanmean(cld_fract(i,:,k1:k2),2)),Z(k1:k2),'k','linewidth',2)
  end

end
xlim([0 0.05])
ylim([0 3500])
legend('L1','L2')
fname=[outpath exp '_cld_fraction.png'];
print ('-dpng',fname)

figure;
colormap(jet)
hold all
set(gca,'fontsize',15);
ylabel('cloud fraction')

for i=1:length(expn)
  if i==1
   plot(datev,squeeze(nanmean(cld_fract(i,:,k1:k2),3)),'b','linewidth',2)
  else
   plot(datev,squeeze(nanmean(cld_fract(i,:,k1:k2),3)),'k','linewidth',2)
  end

end

  xlabel('Local Time (h)');
  datetick('x','hh','keepticks')
legend('L1','L2')
fname=[outpath exp '_ts_cld_fraction.png'];
print ('-dpng',fname)

end % ifload

end % no of experiment
