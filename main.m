clear vars; clc; close all;
% 
L = 100;
H = 2 * L;
rho_bar = 1;
T_H = 40;
T_C = 10;
T_M = (T_H + T_C) / 2;
g = 9.81;
beta = 0.5;
V_C = sqrt(g*beta*(T_H-T_C)*L);
p_0 = rho_bar * V_C^2;

Re = 100;
Pr = 1;

% Initialziation of Grid
ni = 10;           % Number of Cells in X Direction
nj = 10;           % Number of Cells in Y Direction
imax = ni + 2;      % Number of Array Elements in X Direction
jmax = nj + 2;      % Number of Array Elements in Y Direction

deltax = L / ni;
deltay = H / nj;

% Initialization of Fluid parameter arrays
unm1 = zeros(imax, jmax);
vnm1 = zeros(imax, jmax);
pnm1 = zeros(imax, jmax);
Tnm1 = zeros(imax, jmax);

un = ones(imax, jmax);
vn = ones(imax, jmax);
pn = ones(imax, jmax);
Tn = zeros(imax, jmax);

u_star = zeros(imax, jmax);
v_star = zeros(imax, jmax);

p_prime = zeros(imax, jmax);

unp1 = zeros(imax, jmax);
vnp1 = zeros(imax, jmax);
pnp1 = zeros(imax, jmax);
Tnp1 = zeros(imax, jmax);

F1n = zeros(imax, jmax);
F2n = zeros(imax, jmax);
F3n = zeros(imax, jmax);
F1nm1 = zeros(imax, jmax);
F2nm1 = zeros(imax, jmax);
F3nm1 = zeros(imax, jmax);


% Time Iteration
t_sim = 0;
itr_max = 10;
deltat = 1e-5;
for itr = 0:itr_max
    % Left Wall (Wall with Constant Temperature T_H)
    un(1,:) = 0;
    vn(1,:) = - vn(2,:);
    Tn(1,:) = 2 * (T_H - T_M) / (T_H - T_C) - Tn(2,:);
    p_prime(1,:) = p_prime(2,:);

    
    % Right Wall (Wall with Constant Temperature T_C)
    un(imax,:) = 0;
    vn(imax,:) = - vn(imax-1,:);
    Tn(imax,:) = 2 * (T_C - T_M) / (T_H - T_C) - Tn(imax-1,:);
    p_prime(imax,:) = p_prime(imax-1,:);

    % Bottom Wall   (adiabat Wall)
    un(:,1) = - un(:,2);
    vn(:,1) = 0;
    Tn(:,1) = Tn(:,2);
    p_prime(:,1) = p_prime(:,2);

    % Top Wall      (adiabat Wall)
    un(:,jmax) = - un(:,jmax-1);
    vn(:,jmax) = 0;
    Tn(:,jmax) = Tn(:,jmax-1);
    p_prime(:,jmax) = p_prime(:,jmax-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Projection Step
    % calculation of u_star and v_star
    for i = 2 : imax-1
        for j = 2: imax-1
            % RHS Term for X Impulse
            F1n_1 = (((un(i+1, j) + un(i,j))/2)^2 -  ...
                ((un(i, j) + un(i-1,j))/2)^2)* deltay;
                        
            F1n_2 = (((un(i, j) + un(i,j+1))/2) * ((vn(i, j) + vn(i+1,j))/2) - ...
                ((un(i, j-1) + un(i,j))/2) * ((vn(i, j-1) + vn(i+1,j-1))/2)) * deltax;

            F1n_3 = (pn(i+1,j)-pn(i,j))*deltay;

            F1n_4 = (((un(i+1, j) - un(i,j))/deltax) +  ...
                ((un(i, j) - un(i-1,j))/deltax))* deltay;

            F1n_5 = (((un(i, j+1) - un(i,j))/deltay) +  ...
                ((un(i, j) - un(i,j-1))/deltay))* deltax;

            F1n(i,j) = F1n_1 + F1n_2 + F1n_3 + F1n_4 + F1n_5;

            % RHS Term for Y Impulse
            F2n_1 = (((un(i, j) + un(i,j+1))/2) * ((vn(i, j) + vn(i+1,j))/2) - ...
                ((un(i-1, j+1) + un(i-1,j))/2) * ((vn(i-1, j) + vn(i,j))/2)) * deltax;

            F2n_2 = (((vn(i, j+1) - vn(i,j))/2)^2 -  ...
                ((vn(i, j) + vn(i,j-1))/2)^2)* deltay;

            F2n_3 = (pn(i,j+1)-pn(i,j))*deltax;

            F2n_4 = (((vn(i+1, j) - vn(i,j))/deltax) +  ...
                ((vn(i, j) - vn(i-1,j))/deltax))* deltay;

            F2n_5 = (((vn(i, j+1) - vn(i,j))/deltay) +  ...
                ((vn(i, j) - vn(i,j-1))/deltay))* deltax;

            F2n_6 = (Tn(i, j+1) - Tn(i, j)) * deltax * deltay;

            F2n(i,j) = F2n_1 + F2n_2 + F2n_3 + F2n_4 + F2n_5 + F2n_6;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prediction Step Values

            u_star(i,j) = (3/2 * F1n(i,j) - 1/2 * F1nm1(i,j)) * deltat / (deltax * deltay) + un(i,j);
            v_star(i,j) = (3/2 * F2n(i,j) - 1/2 * F2nm1(i,j)) * deltat / (deltax * deltay) + vn(i,j);
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Poisson Solver
    A = zeros(imax*jmax, imax*jmax);
    b = zeros(imax*jmax, 1);
    x = zeros(imax*jmax, 1);
    for i = 2 : imax-1
        for j = 2: imax-1
            index = (i-1) * jmax + j;
            b(index) = (u_star(i, j) - u_star(i-1, j)) * deltay + ...
                       (v_star(i, j) - v_star(i, j-1)) * deltax;
            A(index, index) = -2 *(deltax/deltay + deltay/deltax);
            A(index, index+1) = deltax/deltay;
            A(index, index-1) = deltax/deltay;
            A(index, index+jmax) = deltay/deltax;
            A(index, index-jmax) = deltay/deltax;


        end
    end
    



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correction Step
    % Calculation of Timestep n + 1 
    for i = 2 : imax-1
        for j = 2: imax-1
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

            F3n(i,j) = F3n_1 + F3n_2 +F3n_3 + F3n_4;

            %%%%%%%%%%%%%%%%%%%%%
            % Calculation of n+1 Values
            Tnp1(i,j) = (3/2 * F3n(i,j) - 1/2 * F3nm1(i,j)) * deltat / (deltax * deltay) + Tn(i,j);





        end
    end





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    t_sim = t_sim + deltat;

end






