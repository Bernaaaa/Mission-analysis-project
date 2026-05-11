clear
clc
close all

% : Characterize the transfer orbit with the given data
% goal: Characterize the transfer orbit with the given data

% INPUT:
% % Initial orbital parameters [𝑎𝑖, 𝑒𝑖, 𝑖𝑖, Ω𝑖, 𝜔𝑖 ]
% % Final orbit parameters [𝑎𝑓, 𝑒𝑓, 𝑖𝑓, Ω𝑓, 𝜔𝑓 ]

% OUTPUT:
 % [𝑎𝑇, 𝑒𝑇, 𝑖𝑇, Ω𝑇, 𝜔𝑇] : keplerian paramters of the transfert orbit (shape, orientation, plane of the orbit)
 % 𝜃1𝑇, 𝜃2𝑇 : initial and final points on the transfer orbit
 % 𝜃1𝑖, 𝜃2𝑓 : initial and final(optimal) points on the initial and final orbits

%----------------------------------------------------------------------------------------------------------

%Degree of freedom of the problem: 3 (om_tn, th_i, th_f)
% thi, thf, om_tn (pericenteer argument of the transfer orbit)

% Constants:
mu = 1.32712440041279419e+11; % m^3/s^2 from https://ssd.jpl.nasa.gov/astro_par.html (NASA JPL)

% 1) Definition of the state vector (position and velocity) at points 1 and 2:

%26 3908 Nyx (1980 PA)

%Intial Data: 
r_xi = -7090.590200;
r_yi = -5612.557300;
r_zi = 3948.902900;
v_xi = 5.698000;
v_yi = -5.995000;
v_zi = 1.710000;

r1i = [r_xi; r_yi; r_zi];
v1i = [v_xi; v_yi; v_zi];

[a_i, e_i, i_i, OM_i, om_i] = car2par(r1i, v1i, mu);

%Final Data: 
a_f = 1.927925;
e_f = 0.459072;
i_f = 2.19 ;
OM_f = 261.19 ;
om_f = 126.66 ;
mass = 1.0472E+12 ;
D = 1.00;

%1) scelta om finale = 0 per avere pericentro in direzione del sole, 180 per avere apocentro in direzione del sole (?)


%2) Initialize the vectors for the variation of the parameters of the transfer orbit:
om_tn = linspace(0, 2*pi, 100);  % pericenter argument of the transfer orbit varying from 0 to 360 degrees
th_i_val = linspace(0, 2*pi, 100); % true anomaly of the initial orbit varying from 0 to 360 degrees
th_f_val = linspace(0, 2*pi, 100); % true anomaly of the final orbit varying from 0 to 360 degrees

% Create 3D grid for the parameters
[om_tn_grid, th_f_val_grid, th_i_val_grid] = meshgrid(om_tn, th_f_val, th_i_val);

%Create a 3D matrix to store the total Δv for each combination of parameters
DV = zeros(length(om_tn), length(th_f_val), length(th_i_val));

%Loop for the variation of the parameters of the transfer orbit:
for j=1:length(om_tn)
    om_t = om_tn(j);

    % Loop for the variation of the true anomaly of the initial orbit:
    for k =1:length(th_i_val)
        th_i = th_i_val(k);
        
        % Loop for the variation of the true anomaly of the final orbit:
        for m=1:length(th_f_val)
            th_f = th_f_val(m);
            
            % Definition of the state vector (position and velocity) at points 1 and 2 for each combination of parameters:
            [r2f, v2f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);

            r2=r2f;
            r1=r1i;

            % 3) Definition of the transfer orbital plane, identified by the normal unit vector:
            h_t=cross(r1, r2)/norm(cross(r1, r2));

            % 4) Determination of the inclination of the transfer orbit:
            i_t=acos(dot(h_t, [0 0 1]));

            % 5) Determination of the line of nodes:
            N_t=cross([0 0 1], h_t)/norm(cross([0 0 1], h_t));

            % 6) Determination of the right ascension of the ascending node (RAAN):
            if N_t(2)>=0
                OM_t=acos(dot(N_t, [1 0 0]));
            else
                OM_t=2*pi-acos(dot(N_t, [1 0 0]));
            end

            % ROTATION MATRICES:
            R3_OM = [cos(OM_t), sin(OM_t), 0;
                    -sin(OM_t), cos(OM_t), 0; 
                    0, 0, 1];

            R1_i = [1, 0, 0; 
                    0, cos(i_t), sin(i_t); 
                    0, -sin(i_t), cos(i_t)];

            R3_om = [cos(om_t), sin(om_t), 0; 
                    -sin(om_t), cos(om_t), 0; 
                    0, 0, 1];

            % Transformation matrix from the inertial frame to the perifocal frame of the transfer orbit:
            T = R3_om * R1_i * R3_OM;

            % Calculation of the position at point 1 and point 2 in the perifocal system of the transfer orbit:
            r1_t=(T * r1);
            r2_t=(T * r2);

            % Determination of the true anomaly at point 1 and point 2 in the transfer orbit:    
            co_1tn = r1_t(1)/norm(r1_t);
            si_1tn = r1_t(2)/norm(r1_t);
            theta_1t = atan2(si_1tn, co_1tn);

            co_2tn = r2_t(1)/norm(r2_t);
            si_2tn = r2_t(2)/norm(r2_t);
            theta_2t = atan2(si_2tn, co_2tn);   

            % d) Determination of the orbital parameters (semi-major axis and eccentricity) using the conic equation:
            r1 = norm(r1);
            r2 = norm(r2);
            e_t = (r1-r2)/(r1*cos(theta_1t) - r2*cos(theta_2t));
            a_t = r1 *(1+e_t*cos(theta_1t))/(1-e_t^2);

            % e) Determination of the velocity at the two maneuver points for each 𝜔𝑇,𝑛 :
            [~,v1t] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i, mu);
            [~,v2t] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);

            % f) Calculation of the Δ𝑣 for each 𝜔𝑇,𝑛
            delta_v1 = norm(v1t - v1i);
            delta_v2 = norm(v2f - v2t);

            % storage of the total Δv for each combination of parameters in the 3D matrix:
            DV(k, j, m) = delta_v1 + delta_v2;  % N.B. the order of the indices in DV is (th_i, om_tn, th_f) to be consistent with the order of the loops, but the order of the parameters in the grid is (om_tn, th_f, th_i) because of the way meshgrid works. This will be important for the visualization of the results.

        end
    end
end

%% VISUALIZATION OF THE RESULTS:

figure()
slice(om_tn_grid, th_f_val_grid, th_i_val_grid, DV, ...
      [], [], th_i_val(20)) % slice at th_i = 20 (arbitrary choice, can be changed to visualize the results for different values of th_i)
colorbar
shading interp
xlabel('\omega_t')
ylabel('\theta_f')
zlabel('\theta_i')

figure()
slice(om_tn_grid, th_f_val_grid, th_i_val_grid, DV, ...
      [], th_f_val(20), []); % slice at th_f = 20 (arbitrary choice, can be changed to visualize the results for different values of th_f)
colorbar
shading interp
xlabel('\omega_t')
ylabel('\theta_f')
zlabel('\theta_i')

figure()
slice(om_tn_grid, th_f_val_grid, th_i_val_grid, DV, ...
      om_tn(20), [], []) ; % slice at om_tn = 20 (arbitrary choice, can be changed to visualize the results for different values of om_tn)
colorbar
shading interp
xlabel('\omega_t')
ylabel('\theta_f')
zlabel('\theta_i')

figure()
contourf(om_tn_grid(:,:,1), th_f_val_grid(:,:,1), DV(:,:,1), 20) % contour plot at th_i = 0 (arbitrary choice, can be changed to visualize the results for different values of th_i)
colorbar
xlabel('\omega_t')
ylabel('\theta_f')  

%8) Identificazione del valore di 𝜔𝑇,𝑛 che minimizza il Δ𝑣 totale (Δ𝑣1 + Δ𝑣2):
%% OPTIMIZATION OF THE TRANSFER ORBIT PARAMETERS:

