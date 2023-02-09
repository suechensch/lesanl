%clear all
close all

exps={'L1','L1sst','L2','L2sst'};
subn={'fig3','fig4'};

outpath='../figures/';
data_path='../data/';
sname=[data_path 'cld_iso.mat'];
load(sname)

i1=1;i2=ma;
j1=1;j2=na;
yslice=220;
k1=1;
k2=199;
k12=k2-k1+1;

[lono,lato]=meshgrid(1:ma,1:na);
[slon,slat,sz]=meshgrid(i1:i2,j1:j2,k1:k2);

for j=1:j2-j1+1
  for i=1:i2-i1+1
for k=1:k12
  kk=k1+k-1;
  jj=j1+j-1;
  ii=i1+i-1;
  slon(j,i,k)=lono(jj,ii);
  slat(j,i,k)=lato(jj,ii);
  sz(j,i,k)=Z(kk);
end
end
end


mylabel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
isoval=0.0001;  % isosurface value

bp=[1,3];
ep=[2,4];
for nf=1:2
  figure
  ic=0;

for t=1:2
 for np=bp(nf):ep(nf)
    ic=ic+1;
    subplot(2,2,ic)
    hold all;
    fld=squeeze(cld(np,t,j1:j2,i1:i2,k1:k2));
% display isovalue of 1 g/kg
    colormap(flipud(gray(100)))
    s =isosurface(slon,slat,sz,fld,isoval,slon);
    p=patch(s);
    set(p,'FaceColor',[0.97 0.97 0.97]);
    set(p,...
            'EdgeColor','none', ...
            'AmbientStrength',.2,...
            'SpecularStrength',.7,...
            'DiffuseStrength',.4); 
    caxis([0 5])
    colorbar('southoutside')
    isonormals(slon,slat,sz,fld,p)

    camlight('right'); 
    lighting gouraud;
    whitebg('white')
    set(gcf,'color','white')
    grid on
    grid minor

    if ic < 3
      titlestring=[mylabel{ic} ' 6 am (' exps{np} ': ' num2str(exp_tt(np,t),'%02d') ' h  FCST)'];
    else
      titlestring=[mylabel{ic} ' 6 pm (' exps{np} ': ' num2str(exp_tt(np,t),'%02d') ' h  FCST)'];
    end
    title(titlestring,'fontsize',12)
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';

    xlabel('X');
    ylabel('Y');
    zlabel('Z (km)');
    fig = gcf;
    fig.InvertHardcopy = 'off';

    view(45,25); axis tight
    zlim([0 2.5])

 end % np
end % end t
    fname=[outpath subn{nf} '_cld_iso.png'];
    saveas(gcf,fname)
    print('-dpng',fname)
end % end nf: number of figures

