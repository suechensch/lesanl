%
clear all
close all

exps={'L1','L1sst','L2','L2sst'};
opath=['../data/'];
subn='fig8';
npt2=23;
npt=21;
npt3=19;
nexp=length(exps);
mcld_en_base(1:npt2)=nan;
mcld_en_top(1:3,1:npt2)=nan;
mcld_de_base(1:3,1:npt2)=nan;
mcld_de_top(1:3,1:npt2)=nan;
tmp(1:nexp,1:30)=nan;
tt=[1:npt];
ebase(1:4,1:npt2)=nan;
dbase(1:4,1:npt2)=nan;
etop(1:4,1:npt2)=nan;
dtop(1:4,1:npt2)=nan;
datehv(1:4,1:npt2)=nan;

lt=4; % model data starting forecast hour
add_lt=8+lt;
for np=1:nexp
  expn=exps{np};

if strcmp(expn,'L2') == 1 | strcmp(expn,'L2sst') == 1
  for t=1:npt
    datehv(np,t)=datenum(datevec(datenum(2019,09,27,23+(t-1)+add_lt,0,0)));
  end
else
   npt=21;
  for t=1:npt
    datehv(np,t)=datenum(datevec(datenum(2019,09,27,21+t-1+add_lt,0,0)));
  end
end
   ntime=npt*12;
sfile=[opath exps{np} '_mse_cld_entr_.mat'];
load(sfile)

ic2=0;
[nt2,nx,ny]=size(cld_en_base);
en_base_2d=reshape(cld_en_base,[nt2,nx*ny]);
de_base_2d=reshape(cld_de_base,[nt2,nx*ny]);
en_top_2d=reshape(cld_en_top,[nt2,nx*ny]);
de_top_2d=reshape(cld_de_top,[nt2,nx*ny]);
for i=1:nt2
  Y_en(i,:)=quantile(en_base_2d(i,:),[0.25,0.975]);
  Y_de(i,:)=quantile(de_base_2d(i,:),[0.25,0.975]);
  Yt_en(i,:)=quantile(en_top_2d(i,:),[0.25,0.975]);
  Yt_de(i,:)=quantile(de_top_2d(i,:),[0.25,0.975]);
end

for i=1:12:ntime
  ic2=ic2+1;
  if ic2 > npt
   break;
  end
  ii=i+11;
  if ii > ntime | ii > nt2
   break
  end
  hr_Y_en(ic2,:,np)=nanmean(Y_en(i:ii,:),1);
  hr_Y_de(ic2,:,np)=nanmean(Y_de(i:ii,:),1);
  hr_Yt_en(ic2,:,np)=nanmean(Yt_en(i:ii,:),1);
  hr_Yt_de(ic2,:,np)=nanmean(Yt_de(i:ii,:),1);
  hr_en_base(ic2,j1:j2,i1:i2)=nanmean(cld_en_base(i:ii,j1:j2,i1:i2),1);  
  hr_de_base(ic2,j1:j2,i1:i2)=nanmean(cld_de_base(i:ii,j1:j2,i1:i2),1);  
  hr_en_top(ic2,j1:j2,i1:i2)=nanmean(cld_en_top(i:ii,j1:j2,i1:i2),1);  
  hr_de_top(ic2,j1:j2,i1:i2)=nanmean(cld_de_top(i:ii,j1:j2,i1:i2),1);  
end % end i

  length(cld_en_base)
  ic3=size(hr_en_base,1);
  mcld_en_base(np,1:ic3)=nanmean(nanmean(hr_en_base,2),3);
  mcld_de_base(np,1:ic3)=nanmean(nanmean(hr_de_base,2),3);
  mcld_en_top(np,1:ic3)=nanmean(nanmean(hr_en_top,2),3);
  mcld_de_top(np,1:ic3)=nanmean(nanmean(hr_de_top,2),3);
  ebase=mcld_en_base;
  dbase=mcld_de_base;
  etop= mcld_en_top;
  dtop=mcld_de_top;
ifshift = 0;
if ifshift == 0
  ebase=mcld_en_base;
  dbase=mcld_de_base;
  etop= mcld_en_top;
  dtop=mcld_de_top;
else
if strcmp(exps{np},'L1') == 1
  ebase=mcld_en_base;
  dbase=mcld_de_base;
  etop= mcld_en_top;
  dtop=mcld_de_top;
elseif strcmp(exps{np},'L1sst') == 1
  ebase=mcld_en_base;
  dbase=mcld_de_base;
  etop= mcld_en_top;
  dtop=mcld_de_top;
% time shif the array by 2 h
elseif strcmp(exps{np},'L2') == 1
  ebase(3,1:npt-2)=mcld_en_base(3,3:npt);
  etop(3,1:npt-2)=mcld_en_top(3,3:npt);
  dbase(3,1:npt-2)=mcld_de_base(3,3:npt);
  dtop(3,1:npt-2)=mcld_de_top(3,3:npt);
elseif strcmp(exps{np},'L2sst') == 1
  ebase(4,1:npt-2)=mcld_en_base(4,3:npt);
  etop(4,1:npt-2)=mcld_en_top(4,3:npt);
  dbase(4,1:npt-2)=mcld_de_base(4,3:npt);
  dtop(4,1:npt-2)=mcld_de_top(4,3:npt);
end   
  
end % ifshift

end % np


%% plot clod base and top ent/det
% conversion factor to km-1
% time shif the array by 2 h
fc=1000;
figure;
colormap(jet)
  hold all;
  set(gca,'fontsize',15)
  plot(datehv(1,1:npt3),fc*ebase(1,1:npt3),'linewidth', 2)
  plot(datehv(2,1:npt3),fc*ebase(2,1:npt3),'r','linewidth', 2)
  datetick('x','hh')
  plot(datehv(3,1:npt3),fc*ebase(3,1:npt3),'k','linewidth', 2)
  plot(datehv(4,1:npt3),fc*ebase(4,1:npt3),'m','linewidth', 2)
  legend('L1','L1sst','L2','L2sst','location','northeast')
  ylabel('Entrainment (km^-^1)')
  xlabel('Local Time (h)')
axis tight
box on;
title('Cloud Base')
oname=[opath exps{1} '_cld_base_ent_' subn '.png'];
print ('-dpng',oname)


figure;
colormap(jet)
  hold all;
  set(gca,'fontsize',15)
  plot(datehv(1,1:npt3),fc*dbase(1,1:npt3),'linewidth', 2)
  plot(datehv(2,1:npt3),fc*dbase(2,1:npt3),'r','linewidth', 2)
  plot(datehv(3,1:npt3),fc*dbase(3,1:npt3),'k','linewidth', 2)
  plot(datehv(4,1:npt3),fc*dbase(4,1:npt3),'m','linewidth', 2)
  datetick('x','hh')
  legend('L1','L1sst','L2','L2sst','location','northeast')
  ylabel('Detrainment (km^-^1)')
  xlabel('Local Time (h)')
axis tight
box on;
title('Cloud Base')
oname=[opath exps{1} '_base_cld_det_' subn '.png'];
print ('-dpng',oname)

figure;
colormap(jet)
  hold all;
  set(gca,'fontsize',15)
  plot(datehv(1,1:npt3),fc*etop(1,1:npt3),'linewidth', 2)
  plot(datehv(2,1:npt3),fc*etop(2,1:npt3),'r','linewidth', 2)
  plot(datehv(3,1:npt3),fc*etop(3,1:npt3),'k','linewidth', 2)
  plot(datehv(4,1:npt3),fc*etop(4,1:npt3),'m','linewidth', 2)
  datetick('x','hh')
  ylabel('Entrainment (km^-^1)')
  xlabel('Local Time (h)')
axis tight
box on;
  legend('L1','L1sst','L2','L2sst','location','northeast')
title('Cloud Top')
oname=[opath exps{1} '_cld_top_ent_' subn '.png'];
print ('-dpng',oname)

%% detrainment
figure;
colormap(jet)
  hold all;
  set(gca,'fontsize',15)
  plot(datehv(1,1:npt3),fc*dtop(1,1:npt3),'linewidth', 2)
  plot(datehv(2,1:npt3),fc*dtop(2,1:npt3),'r','linewidth', 2)
  plot(datehv(3,1:npt3),fc*dtop(3,1:npt3),'k','linewidth', 2)
  plot(datehv(4,1:npt3),fc*dtop(4,1:npt3),'m','linewidth', 2)
  datetick('x','hh')
  ylabel('Detrainment (km^-^1)')
  xlabel('Local Time (h)')
  legend('L1','L1sst','L2','L2sst','location','northeast')
axis tight
box on;
title('Cloud Top')
oname=[opath exps{1} '_cld_top_det_' subn '.png'];
print ('-dpng',oname)

