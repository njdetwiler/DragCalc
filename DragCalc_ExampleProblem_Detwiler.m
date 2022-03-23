%% Drag Calculator Project

% Nicholas Detwiler
% 61009785
% to be submitted 11 MAR 2022

%% description:
% this matlab script will find the total drag acting on an aircraft with
% given dimensions and speed. the dimensions, unless otherwise specified,
% will be listed below in the variables section. this code may or may not
% be adapted to allow user input of parameters; but that is not in the
% assignment and would just be for fun/flex
% unlike someone who is smart, I have decided to just put the entirety of
% the needed code in one file and use ananoymous functions for the looping
% required by the assignment

clear;clc;   % clear things

%% variables/constants
% necessary constants
y = 1.4;        % sp heat ratio
R = 1716.5;     % gas constant (freedom units)
ro_STP = 2.377e-3;  % density @ sea level 

% aircraft geometry
W = 98e3;       % weight lbs
% wing
b = 93.2;       % ft
S_ref = 1000;     % ft^2
tc_w= 0.106;    % thickness to cord ratio 0.180
A_w = 24.5;     % sweep angle (deg) 11.0
TR_w= 0.2;      % wing taper ratio
crt_w = 17.8;   % ft
EWAR = 0.83;    % exposed wing area ratio
AR = (b^2)/S_ref; % aspect ratio

% fuselage
L_f = 107;      % fuselage length 117
D_f = 11.5;     % fuselage diameter 10.0
Swet_f = 3280;  %0.8 * pi * D_f * L_f;   % wetted area of fuselage

% vertical tail
Sexp_vt = 161;  % ft^2
tc_vt= 0.09;    % v tail taper ratio
A_vt =0;        % v tail sweep angle (deg)
TR_vt=0.8;      % v tail taper ratio
crt_vt = 15.5;  % ft

% horizontal tail
Sexp_ht = 261;  % ft^2
tc_ht= 0.09;    % h tail taper ratio
A_ht =31.6;     % h tail sweep angle (deg)
TR_ht=0.35;     % h tail taper ratio
crt_ht = 11.1;  % ft

% pylons
Swet_p = 117;
tc_p = 0.06;
A_p = 0;
TR_p = 1;
c_p = 16.2;

% nacelles
Swet_n = 455;
fineness_n = 5.0;
L_n = 16.8;

% flow characteristics
V = 765;        % ft/s
rho = 8.754e-4;  % lb/ft^3
T  = 400;       % deg Rankine
mu = 3.025e-7;  % lb*s/ft^2
a = sqrt(y * R * T); % speed of sound
M = V / a;      % mach
V_ind = V * sqrt(rho / ro_STP); % indicated flight speed
nu = mu / rho;  % kinematic viscosity

%% functions
% anonymous functions in this section will be used for looping through many
% different velocities in the plotting portion of the assignment as well as
% just reducing how much I need to type in general

% Mach number
Ma1 = @(u1) u1./a;

% Friction Coefficient (used curve fitting tool w R^2 of 0.9993)
Cf1 = @(re1) 0.2355.*re1.^(-0.2893)+0.001027;

% Reynolds Number for given length and velocity
RN1 = @(v1,l1) (v1.*l1)./nu;

% z coefficient for k
z1 = @(m1,A1) ((2-m1.^2).*cosd(A1)) ./ sqrt(1 - (m1.*(cosd(A1).^2)));

% form factor k
k1 = @(tc2,m2,A2) 1 + z1(m2,A2).*tc2 + 100*(tc2.^4);

% fuselage correction factor s
s_fus1 = @(db1) 2.993.*db1.^3 + -3.389.*db1.^2 + 0.1836.*db1 + 0.9968;

% Parasitic Drag coefficient for given C_f, k, S_wet, S_ref
CDp1 = @(Cf1,k1,Sw1,Sr1) (Cf1.*k1.*Sw1)./Sr1;

% Lift Coefficient with given velocity V1, planform area S1
CL1 = @(V1,S1) (W./( (.5*rho*(V1.^2)) * S1) );

% Induced Drag Coefficient with a given C_L,AR,e
CDi1 = @(cl1,ar1,e1) (cl1.^2) ./ (pi * AR .* e1);

% dynamic pressure from given v
q1 = @(v2) 0.5*rho .* (v2.^2);


%% drag calculations

% Part 1: Total Parasitic Drag Coefficient
    % Wing
    Sexp_w = EWAR * S_ref; Swet_w = 2 * 1.02 * Sexp_w; % exposed wing area

    ctp_w = TR_w * crt_w;   % chord of wing tip

    c_MACexp = crt_w - (crt_w - ctp_w) * ((D_f/2) / (b / 2));
    MACexp_w = (2/3) * (crt_w + ctp_w - ...
        ((crt_w * ctp_w) / (crt_w + ctp_w)));   % MAC of exp wing
    Re_w = RN1(V,MACexp_w); Cf_w = Cf1(Re_w);   % reynolds # of wing
    k_w = k1(tc_w,M,A_w);   % form factor

    CDp_w = k_w * Cf_w * Swet_w / S_ref;    % parasitic drag coeff
    f_w = k_w * Cf_w * Swet_w;                    % equiv flat plate area f

    % Fuselage
    Re_f = RN1(V,L_f); Cf_f = Cf1(Re_f);    % reynolds # of fuselage
    width_to_span = D_f / b; corr_fus = s_fus1(width_to_span);  % fuse corr
    fineness_f = L_f / D_f;
    k_f = (61.14*fineness_f -117.8) / ...
        (fineness_f^2 + 46.26*fineness_f + -110.4); % R^2 = 0.9999 fit

    CDp_f = k_f * Cf_f * Swet_f / S_ref;
    f_f = k_f * Cf_f * Swet_f;

    % Horz Tail
    Swet_ht = 2 * 1.02 * Sexp_ht;
    ctp_ht = TR_ht * crt_ht; 
    MACexp_ht = (2/3)*(crt_ht + ctp_ht - ...
        ((crt_ht * ctp_ht)/(crt_ht + ctp_ht)));
    Re_ht = RN1(V,MACexp_ht); Cf_ht = Cf1(Re_ht);
    k_ht = k1(tc_ht,M,A_ht);
    CDp_ht = k_ht * Cf_ht * Swet_ht / S_ref;
    f_ht = k_ht * Cf_ht * Swet_ht;

    % Vert Tail
    Swet_vt = 2 * 1.02 * Sexp_vt;
    ctp_vt = TR_vt * crt_vt; 
    MACexp_vt = (2/3)*(crt_vt + ctp_vt - ...
        ((crt_vt * ctp_vt)/(crt_vt + ctp_vt)));
    Re_vt = RN1(V,MACexp_vt); Cf_vt = Cf1(Re_vt);
    k_vt = k1(tc_vt,M,A_vt);
    CDp_vt = k_vt * Cf_vt * Swet_vt / S_ref;
    f_vt = k_vt * Cf_vt * Swet_vt;

    % Pylons
    k_p = k1(tc_p,M,A_p);
    Re_p = RN1(V,c_p);
    Cf_p = Cf1(Re_p);
    CDp_p = k_p * Cf_p * Swet_p / S_ref;
    f_p = k_p * Cf_p * Swet_p;

    % Nacelles
    k_n = (61.14*fineness_n -117.8) / ...
        (fineness_n^2 + 46.26*fineness_n + -110.4); % R^2 = 0.9999 fit
    Re_n = RN1(V,L_n);
    Cf_n = Cf1(Re_n);
    CDp_n = k_n * Cf_n * Swet_n / S_ref;
    f_n = k_n * Cf_n * Swet_n;

    % TOTAL
    f_tot = 1.1 * (f_w + f_f + f_ht + f_vt + f_p + f_n);
    CDp_tot = f_tot / S_ref;

% Part 2: Induced Drag
    % find lift coefficient
    CL = CL1(V,S_ref);
    % find efficiency factor e
    % CDp_tot ~~ 0.0200 so use e_0200 fitted function (R^2 = 0.9999)
    e = 1.255e-05*AR^3 + 0.0001553*+AR^2 + -0.02069*AR + 0.965;
    % find induced drag coefficient
    CDi = CDi1(CL,AR,e);

% Part 3: Total Incompressible Drag Coefficient and Drag, L/D
CD_tot = CDp_tot + CDi;
D_tot = CD_tot * q1(V) * S_ref;
LD_ratio = W/D_tot;

%% plotting

% velocity vs drag


% velocity vs L/D


