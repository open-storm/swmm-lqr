sys = ss(A,-G,C,D,T); 
[rsys,U] = minreal(sys,[],0);

Ahat = rsys.a;
Bhat = rsys.b;
Chat = rsys.c;
Dhat = rsys.d;

Qtil = U*Q*U';
Qhat = Qtil( 1:length(Ahat),1:length(Ahat) );
Rhat = R;
%Qhat = (100 * rho * Rhat) + .01; % Hardcoding this for total pond control to attempt to force no levels go below desired setpoint

n_mv = size(Ahat,1); % number of manipulated variables

[Khat1,S1] = dlqr(Ahat,Bhat,Qhat,Rhat);
V    = eye(n_mv);
S    = dare(Ahat,Bhat,V'*Qhat*V,Rhat);
Tmat =-(Ahat'-Ahat'*S*Bhat*((R+Bhat'*S*Bhat)\Bhat')-eye(n_mv))\V'*Qhat;
E    = (R+Bhat'*S*Bhat)\Bhat'*Tmat;
Khat = (R+Bhat'*S*Bhat)\Bhat'*S*Ahat;
%disp([norm(Khat1-Khat) norm(S1-S)])
%u=-Khat*x+E*x_ref;

%%
%{
Atil = U*A*U'; %Atil = U*A*inv(U); Atil = U*A/U;
Btil = -T/As(1)*U*G;
Ktil = K*U';

% Confirm Ahat is contained within Atil: True
isequal(Atil(1:length(Ahat),1:length(Ahat)), Ahat)

% Confirm Bhat is the first indices within Btil: True
isequal(Btil(1:length(Ahat),:), Bhat)

% Check if Khat is the same as the first indices of K*U': apparently, not
% quite... (on the order of 1e-12 difference)
isequal(Ktil(:,1:length(Ahat)), Khat)

all(all( abs(Ktil(:,1:length(Ahat)) - Khat) < 1e-8 ))
%}