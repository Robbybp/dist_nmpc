function xprime=cola2mod(t,X,U) 
%
% cola2mod - Model of two identical columns in series, both based on Column
%            A, with three components.  
%
% colamod - This is a nonlinear model of a distillation column with
%           NT-1 theoretical stages including a reboiler (stage 1) plus a
%           total condenser ("stage" NT). The liquid flow dynamics are
%           modelled by a simple linear relationship.
%           Model assumptions: Two components (binary separation); constant
%           relative volatility; no vapor holdup; one feed and two products;
%           constant molar flows (same vapor flow on all stages); 
%           total condenser
%
%           The model is based on column A in Skogestad and Postlethwaite
%           (1996). The model has 82 states.
%
% Inputs:    t    - time in [min].
%            X    - State, the first 41 states are compositions of light
%                   component A with reboiler/bottom stage as X(1) and 
%                   condenser as X(41). State X(42)is holdup in reboiler/
%                   bottom stage and X(82) is hold-up in condenser. 
%            U(1) - reflux L,
%            U(2) - boilup V,
%            U(3) - top or distillate product flow D,
%            U(4) - bottom product flow B,
%            U(5) - feed rate F,
%            U(6) - feed composition, zF.
%            U(7) - feed liquid fraction, qF.
%
% Outputs:   xprime - vector with time derivative of all the states 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------
% The following data need to be changed for a new column.
% These data are for "colmn A".
% Number of stages (including reboiler and total condenser: 
    NT=41; 
% Location of feed stage (stages are counted from the bottom):
    NF=21;
% Relative volatility
    alpha=[2;1.5];
% Nominal liquid holdups

% Data for linearized liquid flow dynamics (does not apply to reboiler and condenser):
    Kuf = 21.65032 ; % constant above feed
 Kbf = 29.65032 ; % constant below feed
 Muw = 0.25     ; % liquid holdup under weir (kmol)
% End data which need to be changed
%------------------------------------------------------------
%X = [
% Splitting the states
% X ( columns, [xa ,xb,m] , tray ) 
x=zeros(2,2,NT);
x(1,1,:)=X(1:NT);                          % Liquid composition from btm to top
x(1,2,:)=X(NT+1:2*NT);  
x(2,1,:)=X(2*NT+1:3*NT);  
x(2,2,:)=X(3*NT+1:4*NT);  
M=zeros(2,NT) ; 
M(1,:)= X(4*NT+1:5*NT)   ;                % Liquid hold up from btm to top
M(2,:) = X(5*NT+1:6*NT);
% U(column, control) 
% Inputs and disturbances
c=1:2 ; 
LT(1) = U(1);                            % Reflux
LT(2) = U(2); 
VB(1) = U(3);                            % Boilup
VB(2) = U(4);  
D(1)  = U(5);                            % Distillate
D(2)  = U(6); 
B(1)  = U(7);                            % Bottoms
B(2)  = U(8); 

%Feed to first column
F(1)  = 1.41;                            % Feedrate
zF(1,1) = .4;                            % Feed composition
zF(1,2) = .2;
qF(1) = 1;                            % Feed liquid fraction
%Feed to second column
F(2) = B(1) ; 
zF(2,1) = x(1,1,1) ; 
zF(2,2)= x(1,2,1) ; 
qF(2) = 1;
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% THE MODEL
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
c=1:2; 
j=1:2;
% Vapor-liquid equilibria
i=1:NT-1;    y(1,1,i)=(alpha(1).*x(1,1,i))./(1+(alpha(1)-1).*x(1,1,i)+(alpha(2)-1).*x(1,2,i));
y(1,2,i)=(alpha(2).*x(1,2,i))./(1+(alpha(1)-1).*x(1,1,i)+(alpha(2)-1).*x(1,2,i));
 y(2,1,i)=(alpha(1).*x(2,1,i))./(1+(alpha(1)-1).*x(2,1,i)+(alpha(2)-1).*x(2,2,i));
y(2,2,i)=(alpha(2).*x(2,2,i))./(1+(alpha(1)-1).*x(2,1,i)+(alpha(2)-1).*x(2,2,i));

% Vapor Flows assuming constant molar flows
i=1:NT-1;    V(1,i)=VB(1).*ones(1,NT-1);
 V(2,i)=VB(2).*ones(1,NT-1);
i=NF:NT-1;   V(1,i)=V(1,i) + (1-qF(1))*F(1);
 V(2,i)=V(2,i) + (1-qF(2))*F(2);

% Liquid flows assuming linearized tray hydraulics with time constant taul
% Also includes coefficient lambda for effect of vapor flow ("K2-effect").
i=2:NF;
L(1,i) = Kbf*((M(1,i) - Muw+ ((M(1,i) - Muw).^2+10.^(-8)).^(.5))/2 ).^1.5;
L(2,i) = Kbf*((M(2,i) - Muw+ ((M(2,i) - Muw).^2+10.^(-8)).^(.5))/2).^1.5;
i=NF+1:NT-1;
L(1,i) = Kuf*( (M(1,i) - Muw+((M(1,i) - Muw).^2+10.^(-8)).^(.5))/2    ).^1.5;
L(2,i) = Kuf*( (M(2,i) - Muw+((M(2,i) - Muw).^2+10.^(-8)).^.5)/2    ).^1.5;
L(1,NT)=LT(1);
L(2,NT)=LT(2);

% Time derivatives from  material balances for 
% 1) total holdup and 2) component holdup

% Column
i=2:NT-1;
dMdt(1,i) = L(1,i+1)         - L(1,i)       + V(1,i-1)         - V(1,i);
dMdt(2,i) = L(2,i+1)         - L(2,i)       + V(2,i-1)         - V(2,i);
dMxdt(1,1,i)= squeeze(L(1,i+1)).*squeeze(x(1,1,i+1))' - squeeze(L(1,i)).*squeeze(x(1,1,i))' + squeeze(V(1,i-1)).*squeeze(y(1,1,i-1))' - squeeze(V(1,i)).*squeeze(y(1,1,i))';
dMxdt(1,2,i)= squeeze(L(1,i+1)).*squeeze(x(1,2,i+1))' - squeeze(L(1,i)).*squeeze(x(1,2,i))' + squeeze(V(1,i-1)).*squeeze(y(1,2,i-1))' - squeeze(V(1,i)).*squeeze(y(1,2,i))';
dMxdt(2,1,i)= squeeze(L(2,i+1)).*squeeze(x(2,1,i+1))' - squeeze(L(2,i)).*squeeze(x(2,1,i))' + squeeze(V(2,i-1)).*squeeze(y(2,1,i-1))' - squeeze(V(2,i)).*squeeze(y(2,1,i))';
dMxdt(2,2,i)= squeeze(L(2,i+1)).*squeeze(x(2,2,i+1))' - squeeze(L(2,i)).*squeeze(x(2,2,i))' + squeeze(V(2,i-1)).*squeeze(y(2,2,i-1))' - squeeze(V(2,i)).*squeeze(y(2,2,i))';


% Correction for feed at the feed stage
% The feed is assumed to be mixed into the feed stage
dMdt(1,NF) = dMdt(1,NF)  + F(1);
dMxdt(1,1,NF)= dMxdt(1,1,NF) + F(1)*zF(1,1);
dMxdt(1,2,NF)= dMxdt(1,2,NF) + F(1)*zF(1,2);
dMdt(2,NF) = dMdt(2,NF)  + F(2);
dMxdt(2,1,NF)= dMxdt(2,1,NF) + F(2)*zF(2,1);
dMxdt(2,2,NF)= dMxdt(2,2,NF) + F(2)*zF(2,2);

% Reboiler (assumed to be an equilibrium stage)
dMdt(1,1) = L(1,2)      - V(1,1)      - B(1);
dMxdt(1,1,1)= L(1,2)*x(1,1,2) - V(1,1)*y(1,1,1) - B(1)*x(1,1,1);
dMxdt(1,2,1)= L(1,2)*x(1,2,2) - V(1,1)*y(1,2,1) - B(1)*x(1,2,1);
dMdt(2,1) = L(2,2)      - V(2,1)      - B(2);
dMxdt(2,1,1)= L(2,2)*x(2,1,2) - V(2,1)*y(2,1,1) - B(2)*x(2,1,1);
dMxdt(2,2,1)= L(2,2)*x(2,2,2) - V(2,1)*y(2,2,1) - B(2)*x(2,2,1);


% Total condenser (no equilibrium stage)
dMdt(1,NT) = V(1,NT-1)         - LT(1)       - D(1);
dMxdt(1,1,NT)= V(1,NT-1)*y(1,1,NT-1) - LT(1)*x(1,1,NT) - D(1)*x(1,1,NT);
dMxdt(1,2,NT)= V(1,NT-1)*y(1,2,NT-1) - LT(1)*x(1,2,NT) - D(1)*x(1,2,NT);
dMdt(2,NT) = V(2,NT-1)         - LT(2)       - D(2);
dMxdt(2,1,NT)= V(2,NT-1)*y(2,1,NT-1) - LT(2)*x(2,1,NT) - D(2)*x(2,1,NT);
dMxdt(2,2,NT)= V(2,NT-1)*y(2,2,NT-1) - LT(2)*x(2,2,NT) - D(2)*x(2,2,NT);

% Compute the derivative for the mole fractions from d(Mx) = x dM + M dx
i=1:NT;
dxdt(1,1,i) = (squeeze(dMxdt(1,1,i)) - squeeze(x(1,1,i)).*squeeze(dMdt(1,i))')./M(1,i)';
dxdt(2,1,i) = (squeeze(dMxdt(2,1,i)) - squeeze(x(2,1,i)).*squeeze(dMdt(2,i))')./M(2,i)';
dxdt(1,2,i) = (squeeze(dMxdt(1,2,i)) - squeeze(x(1,2,i)).*squeeze(dMdt(1,i))')./M(1,i)';
dxdt(2,2,i) = (squeeze(dMxdt(2,2,i)) - squeeze(x(2,2,i)).*squeeze(dMdt(2,i))')./M(2,i)';

% Output
xprime=[ squeeze(dxdt(1,1,:)) ; squeeze(dxdt(1,2,:)); squeeze(dxdt(2,1,:)); squeeze(dxdt(2,2,:)); squeeze(dMdt(1,:))'; squeeze(dMdt(2,:))' ];
