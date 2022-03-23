%% Drag Calculator Project

% Nicholas Detwiler
% 61009785
% to be submitted 13 MAR 2022

%% description:
% this matlab script will find the total drag acting on an aircraft with
% given dimensions and speed. the dimensions, unless otherwise specified
% will be listed below in the variables section. this code may or may not
% be adapted to allow user input of parameters; but that is not in the
% assignment and would just be for fun/flex

clear;clc;   % clears workspace and command window


%% variables/constants
% necessary constants
y = 1.4;        % sp heat ratio
R = 1716.5;     % gas constant (freedom units)
ro_STP = 2.377e-3;  % density @ sea level 

% aircraft geometry
W = 98e3;       % weight lbs
% wing
b = 93.2;       % ft
S_ref = 1000;   % ft^2
tc_w= 0.180;    % thickness to chord ratio 0.180
A_w = 11.0;     % sweep angle (deg) 11.0
TR_w= 0.2;      % wing taper ratio
crt_w = 17.8;   % ft
EWAR = 0.83;    % exposed wing area ratio
AR = (b^2)/S_ref; % aspect ratio

% fuselage
L_f = 117;      % fuselage length 117
D_f = 10.0;     % fuselage diameter 10.0
Swet_f = 0.8 * pi * D_f * L_f;   % wetted area of fuselage

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
Swet_p = 117;   % ft^2
tc_p = 0.06;    % pylons thickness/chord
A_p = 0;        % sweepback angle (0)
TR_p = 1;       % taper ratio (no taper)
c_p = 16.2;     % chord length of pylon

% nacelles
Swet_n = 455;   % total wetted area of nacelles
fineness_n = 5.0;   % fineness ratio 
L_n = 16.8;     % characteristic length 

% flow characteristics
V = 765;        % ft/s
rho = 8.754e-4;  % lb/ft^3
T  = 400;       % deg Rankine
mu = 3.025e-7;  % lb*s/ft^2
a = sqrt(y * R * T); % speed of sound
M = V / a;      % mach
V_ind = V * sqrt(rho / ro_STP); % indicated flight speed
nu = mu / rho;  % kinematic viscosity

% optimization and plotting
V_range = 230:880;    % velocity range over which to determine V_opt
LDymin=5; LDymax=20;    % y scaling values for L/D graph
Dymin=0; Dymax = 20000; % for total drag graph
xmin=200; xmax=900;     % x scaling values


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
z1 = @(m1,A1) ((2-m1.^2).*cosd(A1)) ./ sqrt(1 - ((m1.^2).*(cosd(A1).^2)));

% form factor k
k1 = @(tc2,m2,A2) (1 + (z1(m2,A2).*tc2) + 100*(tc2.^4));

% fuselage correction factor s (not used)
s_fus1 = @(db1) 2.993.*db1.^3 + -3.389.*db1.^2 + 0.1836.*db1 + 0.9968;

% Parasitic Drag coefficient for given C_f, k, S_wet, S_ref
CDp1 = @(Cf1,k1,Sw1,Sr1) (Cf1.*k1.*Sw1)./Sr1;

% Lift Coefficient with given velocity V1, planform area S1
CL1 = @(V1,S1) (W./( (.5*rho*(V1.^2)) .* S1) );

% CL but w input of dynamic presure q
CLq = @(q1,S1) (W ./ (q1 .* S1));
% Induced Drag Coefficient with a given C_L,AR,e
CDi1 = @(cl1,e1) (cl1.^2) ./ (pi * AR .* e1);

% dynamic pressure from given v
q1 = @(v2) 0.5*rho .* (v2.^2);


%% drag calculations

% Part 1: Total Parasitic Drag Coefficient
    % Wing
    Sexp_w = EWAR * S_ref; Swet_w = 2 * 1.02 * Sexp_w; % exposed wing area

    ctp_w = TR_w * crt_w;   % chord of wing tip

    c_exp = crt_w - (crt_w - ctp_w) * ((D_f/2) / (b / 2));
    MACexp_w = (2/3) * (c_exp + ctp_w - ...
        ((c_exp * ctp_w) / (c_exp + ctp_w)));   % MAC of exp wing
    Re_w = RN1(V,MACexp_w); Cf_w = Cf1(Re_w);   % reynolds # of wing
    k_w = k1(tc_w,M,A_w);   % form factor

    CDp_w = k_w * Cf_w * Swet_w / S_ref;    % parasitic drag coeff
    f_w = k_w * Cf_w * Swet_w;              % equiv flat plate area f

    % Fuselage
    Re_f = RN1(V,L_f); Cf_f = Cf1(Re_f);    % reynolds # of fuselage
    %width_to_span = D_f / b;corr_fus = s_fus1(width_to_span);  % fuse corr
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
        %e = 1.255e-05*AR^3 + 0.0001553*+AR^2 + -0.02069*AR + 0.965;
    % CDp_tot ~~ 0.0150 so use e_0150 fitted function (R^2 = 1.0)
    e = 2.049e-05*AR^3 + -0.0001136*AR^2 + -0.01502*AR + 0.967;


    % find induced drag coefficient
    CDi = CDi1(CL,e);

% Part 3: Total Incompressible Drag Coefficient and Drag, L/D
CD_tot = CDp_tot + CDi;
D_tot = CD_tot * q1(V) * S_ref;
LD_ratio = W/D_tot;

fprintf(['unique values: \nt/c_w = %1.3f \tA_w = %2.1f deg ' ...
    '\nL_f = %3.0f ft \tD_f = %2.0f ft \nSwet_f = %5.2f' ...
    '\n\nAt V=765 ft/s, the aircraft will have: \nCDp_tot = %2.4f' ...
    '\nequivalent flat plate area f_tot = %2.3f ft^2' ...
    '\ntotal drag D_tot = %5.1f lb' ...
    '\nand final lift-drag ratio = %2.2f\n\n'], ...
    tc_w, A_w, L_f, D_f, Swet_f, ...
    CDp_tot, f_tot, D_tot, LD_ratio)


%% optimization and plotting

% make vector of dynamic pressure q that varies over V_range
q_range = q1(V_range);
% get vector of CLs from q_range
CL_range = CLq(q_range,S_ref);

% loop thru the various components to find friction coeffs, form factors
% and drag coeffs
% this couldve been done w a nested loop w 6 indeces instead but eh

for i=1:length(V_range)

    % mach at this v
    M_r = Ma1(V_range(i));

    % use mach to get form factor k
    k_w_r = k1(tc_w, M_r, A_w);
    k_ht_r= k1(tc_ht,M_r,A_ht);
    k_vt_r= k1(tc_vt,M_r,A_vt);
    k_p_r = k1(tc_p, M_r, A_w);
    k_f_r = k_f;
    k_n_r = k_n;

    % get all them Reynolds numbers and fric coeffs
    Rew_r = RN1(V_range(i),MACexp_w);   Cfw_r = Cf1(Rew_r); 
    Ref_r = RN1(V_range(i),L_f);        Cff_r = Cf1(Ref_r);
    Reh_r = RN1(V_range(i),MACexp_ht);  Cfh_r = Cf1(Reh_r);
    Rev_r = RN1(V_range(i),MACexp_vt);  Cfv_r = Cf1(Rev_r);
    Rep_r = RN1(V_range(i),c_p);        Cfp_r = Cf1(Rep_r);
    Ren_r = RN1(V_range(i),L_n);        Cfn_r = Cf1(Ren_r);
    
    % get all them para drag coeffs
    CDpw_r = k_w_r * Cfw_r * Swet_w / S_ref;
    CDpf_r = k_f_r * Cff_r * Swet_f / S_ref;
    CDph_r = k_ht_r* Cfh_r * Swet_ht/ S_ref;
    CDpv_r = k_vt_r* Cfv_r * Swet_vt/ S_ref;
    CDpp_r = k_p_r * Cfp_r * Swet_p / S_ref;
    CDpn_r = k_n_r * Cfn_r * Swet_n / S_ref;

    % sum all para drag coeffs
    CDptot_r(i) = 1.1 * (CDpw_r + CDpf_r + CDph_r + CDpv_r + ...
        CDpp_r + CDpn_r);

    % if para drag coeff is closer to 0.0200, use e_0200, else use e_0150
    if (CDptot_r(i) - 0.0175 > 0)
        e_r = 1.255e-05*AR^3 + 0.0001553*+AR^2 + -0.02069*AR + 0.965;
    else
        e_r = 2.049e-05*AR^3 + -0.0001136*AR^2 + -0.01502*AR + 0.967;
    end

    % use e_r to determine the induced drag coeff
    CDi_r(i) = CL_range(i)^2 / (pi * AR * e_r);

    % total drag coeff
    CDtot_r(i) = CDptot_r(i) + CDi_r(i);

    % lift-drag ratio for current index
    LDratio_r(i) = CL_range(i) / CDtot_r(i); 

end

Dtot_r = CDtot_r .* q_range .* S_ref;
Di_r = CDi_r .* q_range .* S_ref;
Dptot_r = CDptot_r .* q_range .* S_ref;

LDmax = max(LDratio_r);
index_opt = find(LDratio_r == LDmax);
V_opt = V_range(index_opt);
V_opt_butinknots = V_opt * 0.5925;

fprintf(['The unique aircraft with the above characteristics \nwill ' ...
    'operate at maximum aerodynamic efficiency of %2.2f \nat or ' ...
    'around a velocity of: %3.0f ft/s or %3.1f knots \n'], ...
    LDmax,V_opt,V_opt_butinknots)

% plot them mfs
subplot(1,2,1)
hold on

plot (V_range,Dptot_r,'m',V_range,Di_r,'c',V_range,Dtot_r,'k')
xlim([xmin,xmax])
ylim([Dymin,Dymax])
grid on
xlabel('V (ft/s)')
ylabel('Drag Force (lbs)')
title ('Drag vs Velocity')
legend('D_p','D_i','D')

hold off


subplot(1,2,2)
hold on

plot(V_range,LDratio_r,'r')
xlim([xmin,xmax])
ylim([LDymin,LDymax])
grid on
xlabel('V (ft/s)')
ylabel('L/D')
title('Lift-Drag Ratio vs Velocity')
%text(690,17.3,'piecewise e correction','FontSize',6)

hold off

