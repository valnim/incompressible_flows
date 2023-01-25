% rho = konst Navier Stokes 
% Kontinuity: delui/delxi = 0 or delu/delx + delv/dely = 0 
% momentum i: delui/dt + del(ui*uj)/delxj = -1/rho * delp/delxi + del
% momentum i: delui/dt + del(ui*uj)/delxj = -1/rho * delp/delxi 
% + deltauij/delxj + fbi
% tij = 2 ny Sij = ny*(delui/delxj + deluj/delxi)

% 2 Step Method until steady state solution (Projekt / Correction)
% Segregated Scheme:
% Start from known time level n
% We know: uin, pn
% Projection Step: advance solution to new time point n + 1
% Solve Navier STokes using pressure from previous level n
% (ui* - uin)/delta t = - 1/rho * delpn/delxi - advection n + diff n + fbn

% ui* = uin + delta t * F (RHS)
% 1st order in time: rhs only computed on level n
% 2nd order in time: (Adams Bashford): F = 3/2(RHS)n - 1/2 (RHS)n-1

% Second Step: Enforce Continuity
% (uin+1 - uin)/delta t = - 1/rho * delp'/delxi
% pn+1 = pn + p'
% deluin+1/deltt != 0
% Derivate 2nd equation by del/delxi
% deluin+1/delxi - delui*/delxi = - delta t / rho * del^2p'/delxi/delxi
% Simplify Poisson Equation:
% del^2p'/delxi^2 = delui*/delxi * rho/delta t

% Solve poisson for p' --> solve uin+1

% 1. solve 1 euqation for ui*
% 2. solve poisson for p'
% 3. solve 2 equation for unin+1 
% 4. new pressure pn+1 = pn + p'
% Advance in time until steady state

% Discretization of F(RHS): 2nd order staggered grid
% Staggered grid to avoid checkerboarding
% Pressure is defined in cell middle
% Velocities are defined in cell intefaces
% Cartesian coordinates xy uv
% delu/delx + delv/dely = 0
% delu/delt + deluu/delx + deluv/dely = - 1/rho * delp/delx + 
%           + ny (del^2u/delx^2 + del^2u/dely^2) + fbx
% delv/delt + delvv/dely + deluv/delx = - 1/rho * delp/dely + 
%           + ny (del^2u/delx^2 + del^2u/dely^2) + fby








