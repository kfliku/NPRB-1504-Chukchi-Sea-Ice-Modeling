%1) Enter name of netcdf forcing file to be created.
%   If it already exists it will be overwritten!!.
    forc_file='';

%2) Enter times of river forcings data, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
    river_time=[0:3600:120*86400];             % model time (days)
    Qr_flow(1:length(river_time)) = 600000;    % inflow discharge (m^3/s); 
    C_river = 0                                % sediment concentration in river
    C_dye   = 0                                % dye concentration in river
%
    num_river_times=length(river_time);        % do not change this.

%3) Enter number of vertical sigma levels in model.
%   This will be same value as entered in mod_param.F
    N = 40;

%4) Enter the values of theta_s, theta_b, and Tcline from your *.in file.
    theta_s = 1.0
    theta_b = 0.5
    Tcline =  5.0
    Vtransform = 2
    Vstretching = 3
    hc = 5
    Tcline = hc

%5) Enter value of h, Lm, and Mm.
%   This info can come from a grid file or user supplied here.
%   
%   Are you entering a grid file name (1 = yes, 0 = no)? 
    get_grid = 1;    %<--- put a 1 or 0 here
                                                                 
    if (get_grid)
      grid_file=''    %<-enter name of grid here
%
% Get some grid info, do not change this.
% 
      nc=netcdf(grid_file);
      h=nc{'h'}(:);
      [MP,LP]=size(h);
      close(nc);
      Lm=LP-2;
      Mm=MP-2;
%
    else
      Lm=100;       %<--- else put size of grid here, from mod_param.F
      Mm=20;        %<--- else put size of grid here, from mod_param.F
      LP = Lm+2;    %don't change this.
      MP = Mm+2;    %don't change this.

      % enter depth, same as in ana_grid
      for j=1:MP
        for i=1:LP
          h(j,i)=18-16*(Mm-j)/(Mm-1);
        end
      end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc some grid stuff here - do not change this.
% You go on to step 6.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   L  = Lm+1;
   M  = Mm+1;
   xi_psi  = L;
   xi_rho  = LP;
   xi_u    = L;
   xi_v    = LP;
   eta_psi = M;
   eta_rho = MP;
   eta_u   = MP;
   eta_v   = M;

%
% Don't change this either.  This is from set_scoord.F
% This info is calculated so you know the vertical spacings
% of the grid at startup.
% Go to step 6 now.
%
   [sc_r,Cs_r] = stretching(Vstretching,theta_s,theta_b,hc,N,0,0);
   [sc_w,Cs_w] = stretching(Vstretching,theta_s,theta_b,hc,N,1,0);

%
% Assume zeta starts at 0.
%
%    for j=1:eta_rho
%      for i=1:xi_rho
%        zeta(1,j,i) = 0;
%      end
%    end
%
%  Calc z at rho and w points.
%  Don't change this.
%
%   [z_r] = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,1,h',squeeze(zeta(1,:,:))');
%   [z_w] = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,5,h',squeeze(zeta(1,:,:))');

%   z_r = permute(z_r,[3 2 1]);
%   z_w = permute(z_w,[3 2 1]);
%   Hz = diff(z_w,1,1);
   
if 1
    %6) Enter number of rivers.
    % I am assuming that this is actually the number of point sources, not the
    % number of distinct rivers
    num_rivers=100;

    %7) Initialize river location and direction parameters.
    %   Currently, the direction can be along XI-direction
    %   (river-direction = 0) or along ETA-direction (river_direction > 0).
    %   The mass sources are
    %   located at U- or V-points so the grid locations should range from
    %   1 =< river_Xposition =< L  and  1 =< river_Eposition =< M
    %
    %  'river flag: 1=temp, 2=salt, 3=temp+salt, 4=temp+salt+sed, 5=temp+salt+sed+bio

    river_Xposition=[251:251+num_rivers-1]';     % num_rivers values
    river_Eposition=ones([1 num_rivers])';      % num_rivers values
    river_direction=ones([1 num_rivers])';      % num_rivers values
    river_flag=3*ones([1 num_rivers])';           % num_rivers values
else    
    % Do it automatically (assumes ONE river coming in from the East side)
    disp('Finding river location automatically')
    mask_rho = nc_varget(grid_file,'mask_rho');
    riv_nxi = LP-1; % add river one gridpoint in from the East boundary
    river_Eposition = find(mask_rho(:,LP-1)==1)-1;
    num_rivers = length(river_Eposition);
    river_Xposition = riv_nxi*ones(num_rivers,1);
    river_direction = 0 * ones(num_rivers,1);
                                                          
    river_flag = 4 * ones(num_rivers,1);
end

%8) Initialize river shape.
    for i=1:num_rivers
      for k=1:N
        river_Vshape(k,i)=1/N;
      end
    end

%9) Initialize river flow.
    ramp_u=0*24;              % start ramp UP at RAMP_UP hours
    ramp_time=1.0*24;         % ramp from 0 to 1 over RAMP_TIME hours
    ramp_time_2 = 1.0*24;
    ramp_d=365*24;              % start ramp DOWN at RAMP_DOWN hours
    %fac=min( (0.5*(tanh((river_time/3600-ramp_u)/(ramp_time/2))+1)), ...
    %         (1-(0.5*(tanh((river_time/3600-ramp_d)/(ramp_time_2/2))+1))) );
    fac=(tanh((river_time/3600-ramp_u)/(ramp_time/5)));
    fac2 = ones(size(fac));
    for i=1:num_rivers
        river_transport(:,i)= fac.*Qr_flow/num_rivers;
%        river_transport(:,i)= -fac2.*Qr_flow/num_rivers;
    end
%10) Time series of river temp and salt.
    river_temp = ones(num_river_times,N,num_rivers)*(-1.8);
    river_salt = ones(num_river_times,N,num_rivers)*33;
%   for time=1:num_river_times
%     for k=1:N
%       for i=1:num_rivers
%         river_temp(time,k,i)=5;
%         river_salt(time,k,i)=30;
%       end
%     end
%   end

                                           
%11) Enter number of mud sediments (NCS) and number of sand sediments (NNS).
%   These values should be the same as in mod_param.F
    NCS = 0;   %number of cohesive sed classes
    NNS = 0;   %number of non-cohesive sed classes
    NPT = 0;   %number of inactive passive tracer
%       
% calc sed parameters. Do not alter.
%   
   NAT=2;  %assume temp + salt are active
   NST = NCS + NNS;     % total number of sed tracers.
   NT = NAT+NST+NPT;        % total number of tracers.
%   NT = NAT+NST;        % total number of tracers.
    
%12) Sediment class properties (in order, mud first then sand).
%  These values should coincide with your sediment.in file.
  mud_Srho=ones(1,NCS)*2650;        %kg m-3, NCS values
  mud_Sd50=[0.03 0.03]/1000;        %m,      NCS values
  mud_Wsed=[0.1 0.1]/1000;         %m s-1,  NCS values
  mud_tau_ce=[0.05 0.05];           %N m-2,  NCS values
  mud_Erate=[5 5]*1e-5;             %kg m-2 s-1, NCS values
  sand_Srho=ones(1,NNS)*2650;       %kg m-3, NNS values
  sand_Sd50=[1.0]/1000;             %m,      NNS values
  sand_Wsed=[1.0]/1000;             %m s-1,  NNS values
  sand_tau_ce=[0.07];               %N m-2,  NNS values
  sand_Erate=[1]*1e-5;              %kg m-2 s-1, NNS values
%         
% make some combined arrays.  Do not alter.
%       
  Srho=  [mud_Srho,sand_Srho];
  Sd50=  [mud_Sd50,sand_Sd50];
  Wsed=  [mud_Wsed,sand_Wsed];

  tau_ce=[mud_tau_ce,sand_tau_ce];
  Erate= [mud_Erate,sand_Erate];


%13) Time series of river mud and sand.
%
% mud.
if NST ~= 0
  display('Initializing river sediments.')
%
  b = 1;        % from sediment rating curve c = a * Q^b (see Kao & Milliman 2008)
  C_river_temp = C_river * (fac.*Qr_flow/max(Qr_flow)).^b;
%  C_river_temp = C_river * (fac2.*Qr_flow/max(Qr_flow)).^b;
  for idmud=1:NCS
    count=['0',num2str(idmud)];
    count=count(end-1:end);
    for time=1:num_river_times
      for k=1:N
        for i=1:num_rivers
          river_mud_temp(time,k,i) = C_river_temp(time);   % mud conc in river
        end
      end
    end
    eval(['river_mud_',count,' = river_mud_temp;'])
    clear river_mud_temp;
  end
%
% sand.
%
  for isand=1:NNS
    count=['0',num2str(isand)];
    count=count(end-1:end);
    for time=1:num_river_times
                                  
      for k=1:N
        for i=1:num_rivers
          river_sand_temp(time,k,i) = 0;   % sand conc in river
        end
      end
    end
    eval(['river_sand_',count,' = river_sand_temp;'])
    clear river_sand_temp; 
  end
end
%14) Time series of river dye.
% 
if NPT ~= 0
  display('Initializing river inactive passive tracers.')
  C_dye_temp = river_time * 0;
  C_dye_temp(find(river_time/86400 >= 1.5)) = C_dye;
  for idpt=1:NPT
    count=['0',num2str(idpt)];
    count=count(end-1:end);
    for time=1:num_river_times
      for k=1:N
        for i=1:num_rivers
          river_dye_temp(time,k,i) = C_dye_temp(time);   % dye conc in river
        end
      end 
    end
    eval(['river_dye_',count,' = river_dye_temp;'])
    clear river_dye_temp;
  end
end
    
%%%%% plotting some initial fields
    
                                        
%%%%% plotting some initial fields

figure;
%subplot(2,2,1)
plot(river_time/86400,sum(river_transport,2)); grid on;
xlabel('time'); title('discharge');
%subplot(2,2,2)
%plot(river_time/86400,squeeze(river_salt(:,:,5))); grid on;
%xlabel('time'); title('salinity');
%subplot(2,2,3)
%plot(river_time/86400,squeeze(river_mud_01(:,:,5))); grid on;
%xlabel('time'); title('sed. conc');
%subplot(2,2,4)
%plot(river_time/86400,squeeze(river_dye_01(:,:,5))); grid on;
%xlabel('time'); title('dye conc');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END of USER INPUT                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create init file
nc_forc=netcdf(forc_file,'clobber');
if isempty(nc_forc), return, end

%% Global attributes:

disp(' ## Defining Global Attributes...')
nc_forc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);
nc_forc.type = ncchar('Initialization file from create_roms_init.m');

%% Dimensions:

disp(' ## Defining Dimensions...')
 
nc_forc('xi_psi') = L;
nc_forc('xi_rho') = LP;
nc_forc('xi_u') = L;
nc_forc('xi_v') = LP;

nc_forc('eta_psi') = M;
nc_forc('eta_rho') = MP;
nc_forc('eta_u') = MP;
nc_forc('eta_v') = M;

nc_forc('s_rho') = N;
nc_forc('s_w') = N+1;
nc_forc('tracer') = NT;

nc_forc('one') = 1;
nc_forc('two') = 2;
nc_forc('river')=num_rivers;
nc_forc('river_time')=num_river_times;
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

nc_forc{'theta_b'} = ncdouble; %% 1 element.
nc_forc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc_forc{'theta_b'}.units = ncchar('1');
 
nc_forc{'theta_s'} = ncdouble; %% 1 element.
nc_forc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc_forc{'theta_s'}.units = ncchar('1');

nc_forc{'Tcline'} = ncdouble; %% 1 element.
nc_forc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc_forc{'Tcline'}.units = ncchar('meter');

nc_forc{'hc'} = ncdouble; %% 1 element.
nc_forc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc_forc{'hc'}.units = ncchar('meter');

nc_forc{'Cs_r'} = ncdouble('s_rho');
nc_forc{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc_forc{'Cs_r'}.units = ncchar('1');
nc_forc{'Cs_r'}.valid_min = ncdouble(-1);
nc_forc{'Cs_r'}.valid_max = ncdouble(0);
nc_forc{'Cs_r'}.field = ncchar('Cs_r, scalar');

nc_forc{'Cs_w'} = ncdouble('s_w');
nc_forc{'Cs_w'}.long_name = ncchar('S-coordinate stretching curves at W-points');
nc_forc{'Cs_w'}.units = ncchar('1');
nc_forc{'Cs_w'}.valid_min = ncdouble(-1);
nc_forc{'Cs_w'}.valid_max = ncdouble(0);
nc_forc{'Cs_w'}.field = ncchar('Cs_w, scalar');

nc_forc{'sc_r'} = ncdouble('s_rho');
nc_forc{'sc_r'}.long_name = ncchar('S-coordinate at RHO-points');
nc_forc{'sc_r'}.units = ncchar('1');
nc_forc{'sc_r'}.valid_min = ncdouble(-1);
nc_forc{'sc_r'}.valid_max = ncdouble(0);
nc_forc{'sc_r'}.field = ncchar('sc_r, scalar');

nc_forc{'sc_w'} = ncdouble('s_w');
nc_forc{'sc_w'}.long_name = ncchar('S-coordinate at W-points');
nc_forc{'sc_w'}.units = ncchar('1');

nc_forc{'sc_w'}.units = ncchar('1');
nc_forc{'sc_w'}.valid_min = ncdouble(-1);
nc_forc{'sc_w'}.valid_max = ncdouble(0);
nc_forc{'sc_w'}.field = ncchar('sc_w, scalar');

nc_forc{'river'} = ncdouble('river');
nc_forc{'river'}.long_name = ncchar('river_runoff identification number');
nc_forc{'river'}.units = ncchar('nondimensional');
nc_forc{'river'}.field = ncchar('num_rivers, scalar');

nc_forc{'river_time'} = ncdouble('river_time');
nc_forc{'river_time'}.long_name = ncchar('river_time');
nc_forc{'river_time'}.units = ncchar('seconds');
nc_forc{'river_time'}.field = ncchar('river_time, scalar, series');

nc_forc{'river_Xposition'} = ncdouble('river');
nc_forc{'river_Xposition'}.long_name = ncchar('river runoff  XI-positions at RHO-points');
nc_forc{'river_Xposition'}.units = ncchar('scalar');
nc_forc{'river_Xposition'}.field = ncchar('river runoff XI position, scalar, series');

nc_forc{'river_Eposition'} = ncdouble('river');
nc_forc{'river_Eposition'}.long_name = ncchar('river runoff ETA-positions at RHO-points');
nc_forc{'river_Eposition'}.units = ncchar('scalar');
nc_forc{'river_Eposition'}.field = ncchar('river runoff ETA position, scalar, series');

nc_forc{'river_direction'} = ncdouble('river');
nc_forc{'river_direction'}.long_name = ncchar('river runoff direction, XI=0, ETA>0');
nc_forc{'river_direction'}.units = ncchar('scalar');
nc_forc{'river_direction'}.field = ncchar('river runoff direction, scalar, series');

nc_forc{'river_Vshape'} = ncdouble('s_rho','river');
nc_forc{'river_Vshape'}.long_name = ncchar('river runoff mass transport vertical profile');
nc_forc{'river_Vshape'}.units = ncchar('scalar');
                                                          
nc_forc{'river_Vshape'}.field = ncchar('river runoff vertical profile, scalar, series');

nc_forc{'river_transport'} = ncdouble('river_time','river');
nc_forc{'river_transport'}.long_name = ncchar('river runoff mass transport');
nc_forc{'river_transport'}.units = ncchar('meter^3/s');
nc_forc{'river_transport'}.field = ncchar('river runoff mass transport, scalar, series');

nc_forc{'river_flag'} = ncdouble('river');
nc_forc{'river_flag'}.long_name = ncchar('river flag, 1=temp, 2=salt, 3=temp+salt, 4=temp+salt+sed, 5=temp+salt+sed+bio');
nc_forc{'river_flag'}.units = ncchar('nondimensional');
nc_forc{'river_flag'}.field = ncchar('river flag, scalar, series');

nc_forc{'river_temp'} = ncdouble('river_time','s_rho','river');
nc_forc{'river_temp'}.long_name = ncchar('river runoff potential temperature');
nc_forc{'river_temp'}.units = ncchar('Celsius');
nc_forc{'river_temp'}.field = ncchar('river temperature, scalar, series');

nc_forc{'river_salt'} = ncdouble('river_time','s_rho','river');
nc_forc{'river_salt'}.long_name = ncchar('river runoff salinity');
nc_forc{'river_salt'}.units = ncchar('PSU');
nc_forc{'river_salt'}.field = ncchar('river salinity, scalar, series');

if NPT ~= 0
for idpt=1:NPT
  count=['00',num2str(idpt)];
  count=count(end-1:end);
  eval(['nc_forc{''river_dye_',count,'''} = ncdouble(''river_time'', ''s_rho'', ''river'');'])
  eval(['nc_forc{''river_dye_',count,'''}.long_name = ncchar(''river runoff dye concentration, dye ',count,''');'])
  eval(['nc_forc{''river_dye_',count,'''}.units = ncchar(''kilogram meter-3'');'])
  eval(['nc_forc{''river_dye_',count,'''}.time = ncchar(''river_time'');'])
  eval(['nc_forc{''river_dye_',count,'''}.field = ncchar(''river runoff dye_',count,', scalar, series'');'])
end
end

for mm=1:NCS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_forc{''river_mud_',count,'''} = ncdouble(''river_time'', ''s_rho'', ''river'');'])
  eval(['nc_forc{''river_mud_',count,'''}.long_name = ncchar(''river runoff suspended sediment concentration, size class ',count,''');'])
  eval(['nc_forc{''river_mud_',count,'''}.units = ncchar(''kilogram meter-3'');'])
  eval(['nc_forc{''river_mud_',count,'''}.time = ncchar(''river_time'');'])
  eval(['nc_forc{''river_mud_',count,'''}.field = ncchar(''river runoff mud_',count,', scalar, series'');'])
end
for mm=1:NNS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_forc{''river_sand_',count,'''} = ncdouble(''river_time'', ''s_rho'', ''river'');'])
  eval(['nc_forc{''river_sand_',count,'''}.long_name = ncchar(''river runoff suspended sediment concentration, size class ',count,''');'])
  eval(['nc_forc{''river_sand_',count,'''}.units = ncchar(''kilogram meter-3'');'])
  eval(['nc_forc{''river_sand_',count,'''}.time = ncchar(''river_time'');'])
  eval(['nc_forc{''river_sand_',count,'''}.field = ncchar(''river runoff sand_',count,', scalar, series'');'])
end
  
  
%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')
  
nc_forc{'theta_s'}(:) = theta_s;
nc_forc{'theta_b'}(:) = theta_b;
nc_forc{'Tcline'}(:)  = Tcline;
nc_forc{'Cs_r'}(:) = Cs_r
nc_forc{'Cs_w'}(:) = Cs_w;
nc_forc{'sc_w'}(:) = sc_w;
nc_forc{'sc_r'}(:) = sc_r;
nc_forc{'hc'}(:) = hc;

nc_forc{'river_time'}(:) = river_time;
nc_forc{'river_Xposition'}(:) = river_Xposition;
nc_forc{'river_Eposition'}(:) = river_Eposition;
nc_forc{'river_direction'}(:) = river_direction;
nc_forc{'river_Vshape'}(:) = river_Vshape;
nc_forc{'river_transport'}(:) = river_transport;
nc_forc{'river_flag'}(:) = river_flag;
nc_forc{'river_temp'}(:) = river_temp;
nc_forc{'river_salt'}(:) = river_salt;

if NPT ~= 0
for idpt=1:NPT
  count=['00',num2str(idpt)];
  count=count(end-1:end);
  eval(['nc_forc{''river_dye_',count,'''}(:)     = river_dye_',count,';'])           % dye conc in river
end
end

for mm=1:NCS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_forc{''river_mud_',count,'''}(:)     = river_mud_',count,';'])           %sed conc in water column
end
for mm=1:NNS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_forc{''river_sand_',count,'''}(:)     = river_sand_',count,';'])           %sed conc in water column
end


%close file
close(nc_forc)

