clearvars; clc; close all;

% Temperature Definition
T_H = 40;
T_C = 10;
T_M = (T_H + T_C) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialziation of Grid
ni = 50;           % Number of Cells in X Direction
nj = 50;           % Number of Cells in Y Direction
imax = ni + 2;      % Number of Array Elements in X Direction
jmax = nj + 2;      % Number of Array Elements in Y Direction

deltax = 1 / ni;
deltay = 2 / nj;

% Initialization of Fluid parameter arrays
unm1 = zeros(imax, jmax);
vnm1 = zeros(imax, jmax);
pnm1 = ones(imax, jmax);
Tnm1 = zeros(imax, jmax);

un = zeros(imax, jmax);
vn = zeros(imax, jmax);
pn = ones(imax, jmax);
Tn = zeros(imax, jmax);

u_star = NaN(imax, jmax);
v_star = NaN(imax, jmax);

p_prime = NaN(imax, jmax);

unp1 = NaN(imax, jmax);
vnp1 = NaN(imax, jmax);
pnp1 = NaN(imax, jmax);
Tnp1 = NaN(imax, jmax);

F1n = zeros(imax, jmax);
F2n = zeros(imax, jmax);
F3n = zeros(imax, jmax);
F1nm1 = zeros(imax, jmax);
F2nm1 = zeros(imax, jmax);
F3nm1 = zeros(imax, jmax);

index = @(i,j) (i-1) * jmax + j;

deltat = zeros(imax, jmax);

ifix = 4;
jfix = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficent Matrix for Poisson Equation
A = eye(imax*jmax, imax*jmax);
for i = 2 : imax-1
    for j = 2: jmax-1
        if (i~= ifix || j ~= jfix)
            idx = index(i,j);
            k1 = deltax/deltay;
            k2 = deltay/deltax;
            A(idx, idx) = -2 *(k1 + k2);
            A(idx, idx+1) = k1;
            A(idx, idx-1) = k1;
            A(idx, idx+jmax) = k2;
            A(idx, idx-jmax) = k2;
        end
    end
end

% Pressure Boundaries
for i = 2:imax-1
    for j = 1:jmax
        idx = index(i,j);
        if j == 1
            A(idx, idx+1) = -1;
        elseif j == jmax
            A(idx, idx-1) = -1;
        end
    end
end

for i = 1:imax
    for j = 2:jmax-1
        idx = index(i,j);
        if i == 1
            A(idx, idx+jmax) = -1;
            
        elseif i == imax
            A(idx, idx-jmax) = -1;
        end
    end
end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot A Matrix
figure(4);
heatmap(A, "Colormap", jet)
title("Koeffizentenmatrix");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
Re = 200;
Pr = 1;
Cfl = 0.35;

alpha_relax = 0.9;
beta_relax = 0.5;

itr_max = 80000;

contiplot = zeros(1, itr_max);
p_primeplot = zeros(1, itr_max);
conti_conv = 3*1e-10;
p_conv = 3*1e-10;

tol = 4e-4;
gs_itr_max = 1e4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
for itr = 1:itr_max
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation of lokal Timestep
    for i = 2 : imax-1
        for j = 2 :jmax-1
            % Timestep calculation
            temp = (abs(un(i,j)) / deltax + ...
                abs(vn(i,j)) / deltay + 2 / Re / deltax^2 + 2 / Re / deltay^2);
            deltat(i,j) = Cfl * temp^-1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Projection Step
    % calculation of u_star
    for i = 2 : imax-2
        for j = 2: jmax-1            
            % RHS Term for X Impulse
            F1n_1 = (((un(i+1, j) + un(i,j))/2)^2 -  ...
                ((un(i, j) + un(i-1,j))/2)^2)* deltay;
                        
            F1n_2 = (((un(i, j) + un(i,j+1))/2) * ((vn(i, j) + vn(i+1,j))/2) - ...
                ((un(i, j-1) + un(i,j))/2) * ((vn(i, j-1) + vn(i+1,j-1))/2)) * deltax;

            F1n_3 = (pn(i+1,j)-pn(i,j))*deltay;

            F1n_4 = (((un(i+1, j) - un(i,j))/deltax) -  ...
                ((un(i, j) - un(i-1,j))/deltax))* deltay;

            F1n_5 = (((un(i, j+1) - un(i,j))/deltay) -  ...
                ((un(i, j) - un(i,j-1))/deltay))* deltax;

            F1n(i,j) = - F1n_1 - F1n_2 - F1n_3 + (F1n_4 + F1n_5) / Re;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prediction Step Values
            u_star(i,j) = (3/2 * F1n(i,j) - 1/2 * F1nm1(i,j)) * deltat(i,j) / (deltax * deltay) + un(i,j);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Projection Step
    % calculation of v_star
    for i = 2 : imax-1
        for j = 2: jmax-2
            % RHS Term for Y Impulse
            F2n_1 = (((un(i, j) + un(i,j+1))/2) * ((vn(i, j) + vn(i+1,j))/2) - ...
                ((un(i-1, j+1) + un(i-1,j))/2) * ((vn(i-1, j) + vn(i,j))/2)) * deltay;

            F2n_2 = (((vn(i, j+1) + vn(i,j))/2)^2 -  ...
                ((vn(i, j) + vn(i,j-1))/2)^2)* deltax;

            F2n_3 = (pn(i,j+1)-pn(i,j))*deltax;

            F2n_4 = (((vn(i+1, j) - vn(i,j))/deltax) -  ...
                ((vn(i, j) - vn(i-1,j))/deltax))* deltay;

            F2n_5 = (((vn(i, j+1) - vn(i,j))/deltay) -  ...
                ((vn(i, j) - vn(i,j-1))/deltay))* deltax;

            F2n_6 = (Tn(i, j+1) + Tn(i, j)) / 2 * deltax * deltay;

            F2n(i,j) = - F2n_1 - F2n_2 - F2n_3 + (F2n_4 + F2n_5) / Re + F2n_6;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prediction Step Values
            v_star(i,j) = (3/2 * F2n(i,j) - 1/2 * F2nm1(i,j)) * deltat(i,j) / (deltax * deltay) + vn(i,j);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set Star Velocity Boundaries
    % Left Wall
    u_star(1,2:jmax-1) = 0;
    % Right Wall
    u_star(imax-1,2:jmax-1) = 0;
    % Bottom Wall
    v_star(2:imax-1,1) = 0;
    % Top Wall
    v_star(2:imax-1,jmax-1) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Poisson Solver
    % 
    b = zeros(imax*jmax, 1);
    x0 = zeros(imax*jmax, 1);
    for i = 2 : imax-1
        for j = 2: jmax-1
            if (i~= ifix || j ~= jfix)
                idx = index(i,j);
                b(idx) = ((u_star(i, j) - u_star(i-1, j)) * deltay + ...
                           (v_star(i, j) - v_star(i, j-1)) * deltax) /  deltat(i,j);
            end
        end
    end
    

    x = GaussSeidel(A, b, x0, tol, gs_itr_max);
    %x = A\b;
    for i = 2 : imax-1
        for j = 2: jmax-1
            idx = index(i,j);
            p_prime(i,j) = x(idx);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correction Step
    % Calculation of Timestep n + 1 
    for i = 2 : imax-1
        for j = 2: jmax-1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % RHS Term for Energy Impulse
            F3n_1 = (un(i,j) * (0.5 * (Tn(i+1,j) + Tn(i,j))) - ...
                un(i-1,j) * (0.5 * (Tn(i,j) + Tn(i-1,j)))) * deltay;

            F3n_2 = (vn(i,j) * (0.5 * (Tn(i,j+1) + Tn(i,j))) - ...
                vn(i,j-1) * (0.5 * (Tn(i,j) + Tn(i,j-1)))) * deltax;

            F3n_3 = ((Tn(i+1,j) - Tn(i,j))/deltax - ...
                (Tn(i,j) - Tn(i-1,j))/deltax)*deltay;

            F3n_4 = ((Tn(i,j+1) - Tn(i,j))/deltay - ...
                (Tn(i,j) - Tn(i,j-1))/deltay)*deltax;

            F3n(i,j) = - F3n_1 - F3n_2 + (F3n_3 + F3n_4) /Re / Pr;

            %%%%%%%%%%%%%%%%%%%%%
            % Calculation of n+1 Values of T
            Tnp1(i,j) = (3/2 * F3n(i,j) - 1/2 * F3nm1(i,j)) * deltat(i,j) / (deltax * deltay) + Tn(i,j);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation of n+1 Values of u
    for i = 2 : imax-2
        for j = 2: jmax-1
            unp1(i,j) = u_star(i,j) - deltat(i,j) * (p_prime(i+1,j) - p_prime(i,j)) / deltax * alpha_relax;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation of n+1 Values of v
    for i = 2 : imax-1
        for j = 2: jmax-2
            vnp1(i,j) = v_star(i,j) - deltat(i,j) * (p_prime(i,j+1) - p_prime(i,j)) / deltay * alpha_relax;
        end
    end
    pnp1 = pn + p_prime * beta_relax;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Change Timestep
    unm1 = un;
    vnm1 = vn;
    pnm1 = pn;
    Tnm1 = Tn;

    un = unp1;
    vn = vnp1;
    pn = pnp1;
    Tn = Tnp1;

    F1nm1 = F1n;
    F2nm1 = F2n;
    F3nm1 = F3n;

    disp(itr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Wall Boundaries
    % Left Wall (Wall with Constant Temperature T_H)
    un(1,2:jmax-1) = 0;
    vn(1,2:jmax-2) = - vn(2,2:jmax-2);
    Tn(1,2:jmax-1) = 2 * (T_H - T_M) / (T_H - T_C) - Tn(2,2:jmax-1);
    pn(1,2:jmax-1) = pn(2,2:jmax-1);
    
    % Right Wall (Wall with Constant Temperature T_C)
    un(imax-1,2:jmax-1) = 0;
    vn(imax,2:jmax-2) = - vn(imax-1,2:jmax-2);
    Tn(imax,2:jmax-1) = 2 * (T_C - T_M) / (T_H - T_C) - Tn(imax-1,2:jmax-1);
    pn(imax,2:jmax-1) = pn(imax-1,2:jmax-1);
    
    % Bottom Wall   (adiabat Wall)
    un(2:imax-2,1) = - un(2:imax-2,2);
    vn(2:imax-1,1) = 0;    
    Tn(2:imax-1,1) = Tn(2:imax-1,2);
    pn(2:imax-1,1) = pn(2:imax-1,2);
    
    % Top Wall      (adiabat Wall)
    un(2:imax-2,jmax) = - un(2:imax-2,jmax-1);
    vn(2:imax-1,jmax-1) = 0;
    Tn(2:imax-1,jmax) = Tn(2:imax-1,jmax-1);
    pn(2:imax-1,jmax) = pn(2:imax-1,jmax-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check Continuity
    conti = NaN(imax,jmax);
    for i = 2 : imax-1
        for j = 2: jmax-1
            conti(i,j) = abs((un(i, j) - un(i-1, j)) * deltay + ...
                       (vn(i, j) - vn(i, j-1)) * deltax);
        end
    end
    contiplot(itr) = max(abs(conti), [], "all");
    p_primeplot(itr) = max(abs(p_prime), [], "all");

    if mod(itr, 200) == 0
        % Plot Temperature
        figure(1);
        Tdisp = Tn(1:imax, 1:jmax)';
        Tdisp = flipud(Tdisp);
        heatmap(Tdisp, "Colormap", jet)
        title("Temperature");
        
        % Plot Velocities
        figure(2);
        x = linspace(1,imax, imax);
        y = linspace(jmax,1, jmax);
        xu = x + 0.5;
        yv = y - 0.5;
        uy = zeros(imax, jmax);
        vx = zeros(imax, jmax);
        quiver(xu, y, flipud(un'), uy);
        grid on;
        
        hold on;
        quiver(x, yv, vx, flipud(vn'));
        hold off;
        title("Velocity");
        % Plot Pressure
        figure(3);
        pdisp = pn(1:imax, 1:jmax);
        pdisp = flipud(pdisp');
        heatmap(pdisp, "Colormap", jet)
        title("Pressure");
        % Plot Continuity Residual
        figure(5);
        iterations = linspace(1,itr,itr);
        semilogy(iterations,contiplot(1:itr));
        title("Continuity residual");
        % Plot Correction Pressure Residual
        figure(6);
        semilogy(iterations,p_primeplot(1:itr));
        title("Correction pressure residual");
        drawnow;
    end
    if contiplot(itr) < conti_conv && itr > 5e2 && p_primeplot(itr) < p_conv
        break;
    end

end
disp('finished calculation');
disp('Continuity residual: ');
disp(contiplot(itr));
disp('Correction Pressure residual: ');
disp(p_primeplot(itr));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Temperature
figure(1);
Tdisp = Tn(1:imax, 1:jmax);
Tdisp = flipud(Tdisp');
heatmap(Tdisp, "Colormap", jet)
title("Temperature");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Velocity components
figure(2);
x = linspace(1,imax, imax);
y = linspace(jmax,1, jmax);
xu = x + 0.5;
yv = y + 0.5;
uy = zeros(imax, jmax);
vx = zeros(imax, jmax);
hold off;
quiver(xu, y, flipud(un'), uy);
hold on;
quiver(x, yv, vx, flipud(vn'));
hold off;
title("Velocity");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Pressure
figure(3);
pdisp = pn(1:imax, 1:jmax);
pdisp = flipud(pdisp');
heatmap(pdisp, "Colormap", jet)
title("Pressure");

% Plot Continuity Residual
figure(5);
iterations = linspace(1,itr,itr);
semilogy(iterations,contiplot(1:itr));
title("Continuity residual");
% Plot Correction Pressure Residual
figure(6);
semilogy(iterations,p_primeplot(1:itr));
title("Correction pressure residual");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot cell velocity
figure(7);
x = linspace(0,1, imax-2);
y = linspace(2,0, jmax-2);
udisp = flipud(un');
vdisp = flipud(vn');
quiver(x, y, udisp(2:imax-1, 2:jmax-1), vdisp(2:imax-1, 2:jmax-1));
title("Velocity");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot streamline
figure(8);
[X,Y] = meshgrid(2:imax-1, 2:jmax-1); 
udisp = un';
vdisp = vn';
streamslice(X,Y,udisp(2:imax-1, 2:jmax-1), vdisp(2:imax-1, 2:jmax-1));
title("Velocity");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculated Nusselt Number
Nu = 0;
for j = 2 : jmax-1
    Nu = Nu + (Tn(1,j) - Tn(2,j))/deltax*deltay;
end
disp(Nu);


