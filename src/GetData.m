%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get asset return data
%
% INPUT
%   T1: Estimation window size
%   smonth: out-of-sample start month
%   emonth: out-of-sample end month
% OUTPUT
%   r: (T1+T2)xN asset returns (in excess of risk-free rate)
%   rf: (T1+T2)x1 risk-free returns
%   r2: T2xN out-of-sample asset returns (in excess of risk-free rate)
%   rf2: T2x1 out-of-sample risk-free returns
%   * T1, T2: estimation window size and out-of-sample period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,rf,r2,rf2] = GetData(T1, smonth, emonth)

% Risk-free rate from FF Factors
data1 = csvread('FF_Factors.csv',1,0);
sidx1 = find(data1(:,1)==smonth)-T1;
eidx1 = find(data1(:,1)==emonth);
data1 = data1(sidx1:eidx1,:);

% 25 FF factor portfolio returns
RISK_FREE = data1(:,6);
data2 = csvread('25_FF.csv',1,0);
sidx2 = find(data2(:,1)==smonth)-T1;
eidx2 = find(data2(:,1)==emonth);
ASSET_RET = data2(sidx2:eidx2,2:end);

% total return -> excess return
ASSET_RET = ASSET_RET-repmat(RISK_FREE, 1, size(ASSET_RET,2));
r = 0.01*ASSET_RET; % asset returns
rf = 0.01*RISK_FREE; % risk-free rates
r2 = 0.01*ASSET_RET(T1+1:end,:); % out-of-sample asset returns
rf2 = 0.01*RISK_FREE(T1+1:end,:); % out-of-sample risk-free rates

end
