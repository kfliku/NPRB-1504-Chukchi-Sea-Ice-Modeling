% construct Central Channel grid.

filename = '';
Lm = ;
Mm = ;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid parameters.

L = Lm+1;
M = Mm+1;
Lp= L +1;
Mp= M +1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct grid, and calculate grid metrics.

% Construct dx_r and dy_r.  Integrate to get x_r and y_r.
dx_r = 1000*ones(1,Lp);    % Grid resolution in x
dy_r = 1000*ones(1,Mp);    % Grid resolution in y

disp(['Minimum dx = ',num2str(min(dx_r))]);
disp(['Minimum dy = ',num2str(min(dy_r))]);
[pm,pn]=meshgrid(1./dx_r,1./dy_r);

% Calculate x_r and y_r
dx=[dx_r(1)./2 0.5.*(dx_r(1:end-1)+dx_r(2:end))];
dy=[dy_r(1)./2 0.5.*(dy_r(1:end-1)+dy_r(2:end))];
[x_r,y_r]=meshgrid(cumsum(dx),cumsum(dy));

% Shift grid so x_u(:,1)=0 and y_v(1,:)=0.
y_r = y_r - y_r(1,1) - (y_r(2,1)-y_r(1,1))/2;
x_r = x_r - x_r(1,1) - (x_r(1,2)-x_r(1,1))/2;

% Calculate dmde and dndx.
dndx = zeros(Mp,Lp);
dmde = zeros(Mp,Lp);

dndx(2:M,2:L) = (1./pn(2:M,3:Lp) - 1./(pn(2:M,1:Lm)))/2;
dmde(2:M,2:L) = (1./pm(3:Mp,2:L) - 1./(pm(1:Mm,2:L)))/2;

% Calculate x_u, etc.
x_u = (x_r(:,1:L) + x_r(:,2:Lp))/2;
y_u = (y_r(:,1:L) + y_r(:,2:Lp))/2;

x_v = (x_r(1:M,:) + x_r(2:Mp,:))/2;
y_v = (y_r(1:M,:) + y_r(2:Mp,:))/2;

x_p = (x_r(1:M,1:L) + x_r(2:Mp,2:Lp))/2;
y_p = (y_r(1:M,1:L) + y_r(2:Mp,2:Lp))/2;

el = y_u(end,1);
xl = x_v(1,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct grid topography.

h = ones(Mp,Lp)*45;
h_c = 350; k_c = 650;
a = 90; b = 200;

% Hyperboloids Channel
for j=150:700
  for i=51:350
    xw = -a*sec(atan((j-k_c)./b))+h_c;
    h(j,i) = 47.5-3.5*tanh((xw-i)/25);
  end
  for i=351:800
    xe = a*sec(atan((j-k_c)./b))+h_c;
    h(j,i) = 47.5-3.5*tanh((i-xe)/25);
  end
end
% South Extent
xi = 1:651;
ys = ceil(0.0014*(xi-350).^2);
ys(ys<2)=2;
for i = 1:651
  dh = (h(150,i)-45)/(150-ys(i));
  for j = ys(i):150
    h(j,i) = h(j-1,i)+dh;
  end
end

h(abs(h-50)<1e-2)=50;
h(h<45)=45; 
h(isnan(h))=45;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create masking grids, Coreolis grid, and angle.

% Masking at RHO-points.
rmask=ones(size(x_r));

for i=2:Lp,
  for j=1:Mp,
    umask(j,i-1)=rmask(j,i)*rmask(j,i-1);
  end,
end,
for i=1:Lp,
  for j=2:Mp,
    vmask(j-1,i)=rmask(j,i)*rmask(j-1,i);
  end,
end,
for i=2:Lp,
  for j=2:Mp,
    pmask(j-1,i-1)=rmask(j,i)*rmask(j,i-1)*rmask(j-1,i)*rmask(j-1,i-1);
  end,
end,

% Coriolis parameter.
f=1.4e-4 ;

% Angle of the grid is zero.
angle=zeros(size(x_r));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create NetCDF file.

% Open file
nc=netcdf(filename,'clobber');
nc.Description = 'Central-Channel Grid';
nc.Author = 'Kofan Lu';
nc.Created = datestr(now);
nc.type = 'Central-Channel GRID file';

% Dimensions
nc('xi_rho')=Lp;
nc('xi_u')  =L;
nc('xi_v')  =Lp;
nc('xi_psi')=L;

nc('eta_rho')=Mp;
nc('eta_u')  =Mp;
nc('eta_v')  =M;
nc('eta_psi')=M;

nc('one')=1;

% Create variables
dims = {'eta_rho'; 'xi_rho'};
nc{'x_rho'}= ncdouble(dims);
nc{'x_rho'}(:,:)=x_r;

dims = {'eta_psi'; 'xi_psi'};
nc{'x_psi'}= ncdouble(dims);
nc{'x_psi'}(:,:)=x_p;

dims = {'eta_u'; 'xi_u'};
nc{'x_u'}= ncdouble(dims);
nc{'x_u'}(:,:)=x_u;

dims = {'eta_v'; 'xi_v'};
nc{'x_v'}= ncdouble(dims);
nc{'x_v'}(:,:)=x_v;

dims = {'eta_rho'; 'xi_rho'};
nc{'y_rho'}= ncdouble(dims);
nc{'y_rho'}(:,:)=y_r;

dims = {'eta_psi'; 'xi_psi'};
nc{'y_psi'}= ncdouble(dims);
nc{'y_psi'}(:,:)=y_p;

dims = {'eta_u'; 'xi_u'};
nc{'y_u'}= ncdouble(dims);
nc{'y_u'}(:,:)=y_u;

dims = {'eta_v'; 'xi_v'};
nc{'y_v'}= ncdouble(dims);
nc{'y_v'}(:,:)=y_v;

dims = {'eta_rho'; 'xi_rho'};
nc{'pm'}= ncdouble(dims);
nc{'pm'}(:,:)=pm;

dims = {'eta_rho'; 'xi_rho'};
nc{'pn'}= ncdouble(dims);
nc{'pn'}(:,:)=pn;

dims = {'eta_rho'; 'xi_rho'};
nc{'dmde'}= ncdouble(dims);
nc{'dmde'}(:,:)=dmde;

dims = {'eta_rho'; 'xi_rho'};
nc{'dndx'}= ncdouble(dims);
nc{'dndx'}(:,:)=dndx;

dims = {'eta_rho'; 'xi_rho'};
nc{'angle'}= ncdouble(dims);
nc{'angle'}(:,:)=angle;

dims = {'eta_rho'; 'xi_rho'};
nc{'mask_rho'}= ncdouble(dims);
nc{'mask_rho'}(:,:)=rmask;

dims = {'eta_psi'; 'xi_psi'};
nc{'mask_psi'}= ncdouble(dims);
nc{'mask_psi'}(:,:)=pmask;

dims = {'eta_u'; 'xi_u'};
nc{'mask_u'}= ncdouble(dims);
nc{'mask_u'}(:,:)=umask;

dims = {'eta_v'; 'xi_v'};
nc{'mask_v'}= ncdouble(dims);
nc{'mask_v'}(:,:)=vmask;

dims = {'eta_rho'; 'xi_rho'};
nc{'h'}= ncdouble(dims);
nc{'h'}(:,:)=h;

dims = {'eta_rho'; 'xi_rho'};
nc{'f'}= ncdouble(dims);
nc{'f'}(:,:)=f;

dims = {'one'};
nc{'el'} = ncdouble(dims);
nc{'el'}(:) = el;

dims = {'one'};
nc{'xl'} = ncdouble(dims);
nc{'xl'}(:) = xl;

dims = {'one'};
nc{'spherical'} = ncchar(dims);
nc{'spherical'}(:) = 'F';

close(nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
plot1=1;

if plot1==1
figure;clf
rmask(find(rmask==0)) = nan;
%pcolor (x_r(2:20:end,2:20:end),y_r(2:20:end,2:20:end),...
%        h(2:20:end,2:20:end).*rmask(2:20:end,2:20:end))
pcolor(x_r,y_r,h.*rmask);
%set(gca,'dataaspectratio',[1 1 1])
axis equal;
hc=colorbar('horiz');
hcl=get(hc,'ylabel');
set(hcl,'string','Depth (m)');
xlabel ('xi distance (m)');
ylabel ('eta distance (m)');
title ('Topography');


return

figure;clf
plot (y_r(2:end-1,:),h(2:end-1,:),'-o');
end

figure;clf
subplot(2,1,1);
plot (x_r./1000,dx_r./1000);hold on;
plot (x_r./1000,dx_r./1000,'o');
xlabel('Along channel distance (km)');
ylabel('Along channel grid spacing (km)');
subplot(2,1,2);
plot (y_r./1000,dy_r./1000);hold on;
plot (y_r./1000,dy_r./1000,'o');
ylim([0 0.5]);
xlabel('Cross channel distance (km)');
ylabel('Cross channel grid spacing (km)');


