function [xw, control_sequence] = reachable_test(A,B,C_h,D,G,runoff,pct_full_init, pct_full_final)

%%%
% precip [=] in/hr
% runoff [=] cfs
% xw     [=] m

%%
% Define initial values and variables
%size(runoff)
%size(B)

% *3.048 converts 1 ft to 10 meters
cfs2cms = 0.02831685;

k  = size(runoff,1);

if numel(pct_full_init) == 1
    x0 = pct_full_init  .*3.048*sum(C_h,2); % 3.048*sum(C_h,2);
else
    x0 = sum(C_h,2);
    x0( find(sum(C_h,2)) ) = pct_full_init*3.048;
end

if numel(pct_full_final) == 1
    xf = pct_full_final .*3.048*sum(C_h,2);
else
    xf = sum(C_h,2);
    xf( find(sum(C_h,2)) ) = pct_full_final*3.048;
end

xw = zeros(size(B,1),1);
Gr = 0 * A;

u  = zeros(size(G,2),k);

% Compute offset from disturbance

for j = 1:k
        xw = xw + A^( k-j )*B*runoff(j,:)'*cfs2cms;
end

xw = xw + (A^k)*x0;

if nargout == 1
    return
end
%disp(xw)

% Calculate desired end state for the control sequence

ksi = xf - xw;


% Compute the Grammian

for j = 1:k
    Gr = Gr + A^(k-j)*G*G'*(A^(k-j))';
end

% Gr is singular...
Gr_inv = pinv(Gr);


%

for j = 1:k
        u(:,j) = G'* (A^( k-j ))'*Gr_inv*ksi;
end

%figure;plot(u')

%disp(u(:,200))

control_sequence = u(:,200);