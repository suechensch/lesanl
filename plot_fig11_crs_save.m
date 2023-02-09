% matlab script to:
% 1. read COAMPS atmos sig filed
% 2. interpolate to the obs location
% 3. plot time-height cross section of the fields
clear
close all

addpath '/p/home/chens/Matlab/functions';
addpath '/p/home/chens/coamps/utility/matlab';
%exp={'cSep27_22'};
expn={'camp2xSep27'};
%expn={'L1'};
%expn={'camp2xSep27','sstSep27','camp2xSep28-t2','sstSep28'};
np=length(expn);
%ii1=160;
%ii2=170;
ii1=150;
ii2=180;
%texp=[7,7,5,5];
texp=[7];

for i=1:np
   outpath=['/p/cwfs/chens/' expn{np} '/matfiles/'];
   kpath{i}=['/p/work1/chens/les/data/' expn{i} '/'];
%dtg='2019092722';
%dtg='2019092721';
if strcmp(expn{i},'camp2xSep28-t2') == 1
   dtg='2019092723';
elseif strcmp(expn{i},'camp2xSep27') == 1 
   dtg='2019092721';
elseif strcmp(expn{i},'sstSep27') == 1
   dtg='2019092721';
elseif strcmp(expn{i},'sstcSep27') == 1
   dtg='2019092722';
elseif strcmp(expn{i},'cSep27_22') == 1
   dtg='2019092722';
elseif strcmp(expn{i},'sstcSep27') == 1
   dtg='2019092722';
elseif strcmp(expn{i},'sstSep28') == 1
   dtg='2019092723';
end
%dtg={'2018030612'};
nest=1;
nesta=num2str(nest,'%01d');
yy=str2num(dtg(1:4));
mm=str2num(dtg(5:6));
dd=str2num(dtg(7:8));
hh=str2num(dtg(9:10));
[data_grid]=read_datahd_ff(dtg, kpath{i});
m=data_grid.nest.nx(nest);
n=data_grid.nest.ny(nest);
%m=180;
%n=150;
ma=num2str(m,'%04d');
na=num2str(n,'%04d');
kka=data_grid.nz;
Z=flipud(data_grid.sigm);
sig1=num2str(data_grid.sigm(1),'%06d');
sig2=num2str(data_grid.sigm(kka),'%06d');
sigw1=num2str(data_grid.sigw(1),'%06d');
sigw2=num2str(data_grid.sigw(kka+1),'%06d');
ztop=data_grid.ztop;
j1=1; j2=n;
i1=1; i2=m;
%ii1=105; % 24h
%ii2=115;

fname=['longit_sfc_000000_000000_' nesta 'a' ma 'x' na '_' dtg '_00000000_fcstfld'];
%lon=readff(kpath{i},fname,1);

fname=['latitu_sfc_000000_000000_' nesta 'a' ma 'x' na '_' dtg '_00000000_fcstfld'];
%lat=readff(kpath{i},fname,1);

fname=['lndsea_sfc_000000_000000_' nesta 'a' ma 'x' na '_' dtg '_00000000_fcstfld'];
%mask=readff(kpath{i},fname,1);
%landpoints=find(mask==1);
%oceanpoints=find(mask==0);
%lakepoints=find(mask==-1);


%plevels=[1000:-10:100];
plevels=Z*0.001;
k1=1;
k2=100;
k3=65;
lm=length(plevels);   

%obs.lat=[0.0];
%obs.lon=[110.0];
%nobs=length(obs.lon);
obsi=125;
obsj=125;

[x2d,lvl2d]=meshgrid(i1:i2,1:kka);
for ik=i1:i2
  for k=1:lm
    x2d(ik,k)=ik;
    lvl2d(ik,k)=plevels(k);
 end
end


%t1=28;
%t1=18+25/60;
%t2=30;  % hour
%t2=28;  % hour
%t2=5/60;
tint=10/60; % hour
%nt=int8((t2-t1)/tint+1);
nt=1;

%% compute julian date

varname={'pcph20','pottmp','cldmix','uuwind','vvwind','cldmix','turbke','htgrte','relhum','wwwind'};
vname={'tpw','pott','qc','u','v','tmix','tke','htgrte','relhum','w'};
cmax=[100,380,0.0002,8,8,20*10^-3,10^-3,1,100,2];
cmin=[0,320,0,-8,-8,0,0,0,0,-1];
nvar=size(vname,2)
%nvar=1;
cld_cut=0.0001; % cloud mixing ratio threshold to sample the cloudy region 
cld_cut_alto=0.0005; % cloud mixing ratio threshold to sample the cloudy region 

ic=0;
%  fname=[varname{3} '_sig_' sig1 '_' sig2 '_' nesta 'a' ma 'x' na '_' dtg{nn} '_00000000_fcstfld']
%       fld0=readff(kpath{nn},fname,kka);
for n_day=1:nt
%    inc=(n_day-1)*time_int;
%    [yy,mm,dd,hh,min,sec]=datevec(datenum(beg_yy,beg_mm,beg_dd,beg_hh+inc,0,0));
%for t=[6,12,18,24,30]
t1=texp(i);
%for t=t1:tint:t1+1
for t=t1:tint:t1
  ic=ic+1;
%    [fy,fm,fd,fh,fmin,fs]=datevec(datenum(yy,mm,dd,hh+t,0,0));
 
%str_cdtg=[num2str(yy,'%4d') num2str(mm,'%02d') num2str(dd,'%02d') num2str(hh,'%02d')]
%str1=num2str(fd,'%02d');

str=num2str(t,'%04d');

%for nv=4:4  
for nv=3:3 
  [fy,fm,fd,fh,fmin,fs]=datevec(datenum(yy,mm,dd,0+t,0,0));
  if t >=24
    fh=fh+24;
  end
  str=[num2str(fh,'%04d') num2str(fmin,'%02d') num2str(fs,'%02d') ]
  datev(ic)=datenum(fy,fm,fd,fh,fmin,fs);
%for nv=1:nvar  
%  if strcmp(varname{nv},'wwwind') == 1
%    fname=[varname{nv} '_sig_' sigw1 '_' sigw2 '_' nesta 'a' ma 'x' na '_' dtg '_' str '_fcstfld']
%    fldw=readff(kpath{nn},fname,kka+1);
%    fld_up=fldw(:,:,2:kka+1);
%    fld_dn=fldw(:,:,1:kka);
%    fld=0.5*(fld_up+fld_dn);
%  else
    fname=[varname{nv} '_sig_' sig1 '_' sig2 '_' nesta 'a' ma 'x' na '_' dtg '_' str '_fcstfld']
    qc=readff(kpath{i},fname,kka);
%    fld(find(fld<cld_cut))=nan;
%  end

% for obsj=50:200
 for obsj=133
       m2d=squeeze(qc(obsj,:,:));
  fname=['wwwind_sig_' num2str(ztop,'%06d') '_000000_' nesta 'a' ma 'x' na '_' dtg '_' str '_fcstfld'];
       fldw=readff(kpath{i},fname,kka+1);
       fld_up=fldw(:,:,2:kka+1);
       fld_dn=fldw(:,:,1:kka);
       wwind=0.5*(fld_up+fld_dn);
       m2d1=squeeze(wwind(obsj,:,:));

     fname=['turbke_sig_' sig1 '_' sig2 '_' nesta 'a' ma 'x' na '_' dtg '_' str '_fcstfld']
     tke=readff(kpath{i},fname,kka);
       m2d2=squeeze(tke(obsj,:,:));

     fname=['uuwind_sig_' sig1 '_' sig2 '_' nesta 'a' ma 'x' na '_' dtg '_' str '_fcstfld']
%     uu=readff(kpath{i},fname,kka);
%       m2d3=squeeze(uu(obsj,:,:));

       figure
       pcolor(x2d(:,k1:k2),lvl2d(:,k1:k2),m2d1(:,k1:k2));
       shading interp;
       caxis([cmin(10) cmax(10)])
       colorbar;
       hold all
%       [C,h]=contour(x2d(:,k1:k2),lvl2d(:,k1:k2),m2d1(:,k1:k2),'b');
%       clabel(C,h,[0.01],'color','b')
%       [C2,h2]=contour(x2d(:,k1:k2),lvl2d(:,k1:k2),m2d1(:,k1:k2),[-1:1:0],'--b');
%       clabel(C2,h2)
%       hold on;
       contour(x2d(:,k1:k3),lvl2d(:,k1:k3),m2d(:,k1:k3),[0.0001:0.0002:0.0015],'showtext','on','color','k','linewidth',1);
%       [C3,h3]=contour(x2d(:,k1:k3),lvl2d(:,k1:k3),m2d(:,k1:k3),[cld_cut],'k','linewidth',2);
%       clabel(C3,h3,[cld_cut_alto],'color','k')
%       hold on;
%       contour(x2d(:,k3:k2),lvl2d(:,k3:k2),m2d(:,k3:k2),'showtext','on','color','r','linewidth',1);
%       [C4,h4]=contour(x2d(:,k3:k2),lvl2d(:,k3:k2),m2d(:,k3:k2),[cld_cut cld_cut],'r','linewidth',2);
%       clabel(C3,h3,[cld_cut],'color','r')
%       caxis([-8 8])
%       quiver(x2d(:,k1:k3),lvl2d(:,k1:k3),m2d3(:,k1:k3),m2d1(:,k1:k3),1,'color','k','LineWidth',1);
       titlestring=[vname{nv} ', ' str 'h'];
       title(titlestring);
       xlabel('Distance (m)');
       ylabel('Height (km)');
       fname=[outpath expn{i} '_' vname{nv} '_cs_' str '.png'];
       print('-dpng',fname)
%%

       figure
       pcolor(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d1(ii1:ii2,k1:k2));
       shading interp;
       caxis([cmin(10) cmax(10)])
       colorbar;
       hold all
       contour(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d(ii1:ii2,k1:k2),[0.0001:0.005:0.01],'k');
       [C,h]=contour(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d2(ii1:ii2,k1:k2),[0:0.02:0.1],'r');
%       clabel(C,h,[0.0],'color','k')
       [C2,h2]=contour(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d1(ii1:ii2,k1:k2),[-0.1:0.05:0],'--w');
%       clabel(C2,h2,[-0.3],'color','k')
%       hold on;
%       contour(x2d(:,k1:k3),lvl2d(:,k1:k3),m2d1(:,k1:k3),'showtext','on','color','k','linewidth',1);
%       [C3,h3]=contour(x2d(:,k1:k3),lvl2d(:,k1:k3),m2d(:,k1:k3),[cld_cut_alto cld_cut_alto],'k','linewidth',2);
%       clabel(C3,h3,[cld_cut_alto],'color','k')
%       hold on;
%       contour(x2d(:,k3:k2),lvl2d(:,k3:k2),m2d(:,k3:k2),'showtext','on','color','r','linewidth',1);
%       [C4,h4]=contour(x2d(:,k3:k2),lvl2d(:,k3:k2),m2d(:,k3:k2),[cld_cut cld_cut],'r','linewidth',2);
%       clabel(C3,h3,[cld_cut],'color','r')
%       caxis([-8 8])
inc=10;
%       quiver(x2d(:,k1:inc:k2),lvl2d(:,k1:inc:k2),m2d3(:,k1:inc:k2),m2d1(:,k1:inc:k2),1,'color','k','LineWidth',1);
%       titlestring=[vname{nv} ', ' str 'h'];
%       title(titlestring);
%       xlabel('Distance (m)');
       ylabel('Height (km)');
       fname=[outpath expn{i} '_' vname{10} '_' num2str(obsj,'%03d') '_cs_' str '.png'];
       print('-dpng',fname)

end % obsj
end % nv
sname=[outpath 'L1_qc_w_tke.mat'];
save(sname,'qc','wwind','tke','str','data_grid','dtg');
end % t
end % n_day

end % end np

    


