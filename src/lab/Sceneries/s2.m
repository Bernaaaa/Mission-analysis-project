clear
clc

% Obiettivo: Caratterizzare l’orbita di trasferimento con i dati forniti

% INPUT:

% % Parametri orbitali iniziali [𝑎𝑖, 𝑒𝑖, 𝑖𝑖, Ω𝑖, 𝜔𝑖 ]
% % Parametri orbitali finali [𝑎𝑓, 𝑒𝑓, 𝑖𝑓, Ω𝑓, 𝜔𝑓 ] assegnati


% OUTPUT:

 % [𝑎𝑇, 𝑒𝑇, 𝑖𝑇, Ω𝑇, 𝜔𝑇] : parametri kepleriani dell’orbita di trasferimento (forma, orientazione, piano dell’orbita)
 % 𝜃1𝑇, 𝜃2𝑇 : punto iniziale e finale sull’orbita di trasferimento
 % 𝜃1𝑖, 𝜃2𝑓 : punto iniziale sull’orbita di partenza e punto finale sull’orbita di arrivo

%----------------------------------------------------------------------------------------------------------

% 1) Definizione del vettore di stato (posizione e velocità) nei punti 1 e 2:

%26 3908 Nyx (1980 PA)
%Dati iniziali
r_xi = -7090.590200;
r_yi = -5612.557300;
r_zi = 3948.902900;
v_xi = 5.698000;
v_yi = -5.995000;
v_zi = 1.710000;

ri = [r_xi; r_yi; r_zi];
vi = [v_xi; v_yi; v_zi];

%Dati finali 
a_f = 1.927925;
e_f = 0.459072;
i_f = 2.19 ;
OM_f = 261.19 ;
om_f = 126.66 ;
mass = 1.0472E+12 ;
D = 1.00;




[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rf, vf, mu);


[r1i, v1i] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i, mu)

[r1f, v1f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu)

r1=r1i;
r2=r2i;

% 2) Definizione del piano orbitale di trasferimento, identificato dal versore normale:

h_t=cross(r1, r2)/norm(cross(r1, r2));

% 3) Calcolo dell’inclinazione dell’orbita di trasferimento:

i_t=acos(dot(h_t, [0 0 1]));

% 4) Determinazione della linea dei nodi:

N_t=cross([0 0 1], h_t)/norm(cross([0 0 1], h_t));

% 5) Calcolo dell’ ascensione retta del nodo ascendente (RAAN):

if N_t(2)>=0
    OM_t=acos(dot(N_t, [1 0 0]));
else
    OM_t=2*pi-acos(dot(N_t, [1 0 0]));
end

% 6) Ciclo for 

om_tn = linspace(0, 2*pi, 100); % Variazione dell'argomento del pericentro da 0 a 360 gradi


for j=1:length(om_tn)
    om_t=om_tn(j);

% a) Definizione della matrice da sistema di riferimento inerziale a perifocale:
    R3_OM = [cos(OM_t), sin(OM_t), 0;
            -sin(OM_t), cos(OM_t), 0; 
            0, 0, 1];
    R1_i = [1, 0, 0; 
            0, cos(i_t), sin(i_t); 
            0, -sin(i_t), cos(i_t)];

    R3_om = [cos(om_t), sin(om_t), 0; 
            -sin(om_t), cos(om_t), 0; 
            0, 0, 1];

    T = R3_om * R1_i * R3_OM;

%  b) Calcolo di posizione nel punto 1 e punto 2 nel sistema perifocale dell’orbita di trasferimento:
    r1_pf=cross(T,r1);
    r2_pf=cross(T,r2);

%  c) Determinazione delle anomalie vere sull’orbita di trasferimento:    

co_1tn = r1_pf(1)/norm(r1_pf);
si_1tn = r1_pf(2)/norm(r1_pf);
theta_1t = atan2(si_1tn, co_1tn);

co_2tn = r2_pf(1)/norm(r2_pf);
si_2tn = r2_pf(2)/norm(r2_pf);
theta_2t = atan2(si_2tn, co_2tn);  



end


