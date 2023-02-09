% matlab script to:
% 1. read COAMPS atmos sig filed
% 2. interpolate to the obs location
% 3. plot time-height cross section of the fields
clear
close all

expn={'L1'};
np=length(expn);
ii1=150;
ii2=180;
texp=[7];

i=1;
nest=1;
outpath=['../figures/'];
kpath{i}=['../data/' expn{i} '/'];
sname=[kpath{i} expn{i} '_qc_w_tke.mat'];
load(sname,'qc','wwind','tke','str','data_grid','dtg');
yy=str2num(dtg(1:4));
mm=str2num(dtg(5:6));
dd=str2num(dtg(7:8));
hh=str2num(dtg(9:10));
m=data_grid.nest.nx(nest);
n=data_grid.nest.ny(nest);
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
plevels=Z*0.001;
k1=1;
k2=100;
k3=65;
lm=length(plevels);   

obsi=125;
obsj=125;

[x2d,lvl2d]=meshgrid(i1:i2,1:kka);
for ik=i1:i2
  for k=1:lm
    x2d(ik,k)=ik;
    lvl2d(ik,k)=plevels(k);
 end
end


nt=1;


vname={'tpw','pott','qc','u','v','tmix','tke','htgrte','relhum','w'};
cmax=[100,380,0.0002,8,8,20*10^-3,10^-3,1,100,2];
cmin=[0,320,0,-8,-8,0,0,0,0,-1];
nvar=size(vname,2)
cld_cut=0.0001; % cloud mixing ratio threshold to sample the cloudy region 
nv=3;

ic=0;
t1=texp(i);
for t=t1
  ic=ic+1;

str=num2str(t,'%04d');

  [fy,fm,fd,fh,fmin,fs]=datevec(datenum(yy,mm,dd,0+t,0,0));
  datev(ic)=datenum(fy,fm,fd,fh,fmin,fs);

 for obsj=133
       m2d=squeeze(qc(obsj,:,:));
       m2d1=squeeze(wwind(obsj,:,:));
       m2d2=squeeze(tke(obsj,:,:));

%%

       figure
       pcolor(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d1(ii1:ii2,k1:k2));
       shading interp;
       caxis([cmin(10) cmax(10)])
       colorbar;
       hold all
       contour(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d(ii1:ii2,k1:k2),[0.0001:0.005:0.01],'k');
       [C,h]=contour(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d2(ii1:ii2,k1:k2),[0:0.02:0.1],'r');
       [C2,h2]=contour(x2d(ii1:ii2,k1:k2),lvl2d(ii1:ii2,k1:k2),m2d1(ii1:ii2,k1:k2),[-0.1:0.05:0],'--w');
       ylabel('Height (km)');
       fname=[outpath expn{i} '_' vname{10} '_' num2str(obsj,'%03d') '_cs_' str '.png'];
       print('-dpng',fname)

end % obsj
end % t


    


