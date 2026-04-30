clear
clc
% : Characterize the transfer orbit with the given data
% obiettivo: Caratterizzare l’orbita di trasferimento con i dati forniti

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
%mu sole
%Dati iniziali
r_xi = -7090.590200;
r_yi = -5612.557300;
r_zi = 3948.902900;
v_xi = 5.698000;
v_yi = -5.995000;
v_zi = 1.710000;

r1i = [r_xi; r_yi; r_zi];
v1i = [v_xi; v_yi; v_zi];

[a_i, e_i, i_i, OM_i, om_i] = car2par(r1i, v1i, mu)

%Dati finali 
a_f = 1.927925;
e_f = 0.459072;
i_f = 2.19 ;
OM_f = 261.19 ;
om_f = 126.66 ;
mass = 1.0472E+12 ;
D = 1.00;
%1) scelta om finale = 0 per avere pericentro in direzione del sole, 180 per avere apocentro in direzione del sole

[r2f, v2f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu)
%2)

r1=r1i;
r2=r2f;

% 3) Definizione del piano orbitale di trasferimento, identificato dal versore normale:

h_t=cross(r1, r2)/norm(cross(r1, r2));

% 4) Calcolo dell’inclinazione dell’orbita di trasferimento:

i_t=acos(dot(h_t, [0 0 1]));

% 5) Determinazione della linea dei nodi:

N_t=cross([0 0 1], h_t)/norm(cross([0 0 1], h_t));

% 6) Calcolo dell’ ascensione retta del nodo ascendente (RAAN):

if N_t(2)>=0
    OM_t=acos(dot(N_t, [1 0 0]));
else
    OM_t=2*pi-acos(dot(N_t, [1 0 0]));
end

% 6) Ciclo for 

om_tn = linspace(0, 2*pi, 100); % Variazione dell'argomento del pericentro da 0 a 360 gradi

    R3_OM = [cos(OM_t), sin(OM_t), 0;
            -sin(OM_t), cos(OM_t), 0; 
            0, 0, 1];
    R1_i = [1, 0, 0; 
            0, cos(i_t), sin(i_t); 
            0, -sin(i_t), cos(i_t)];

% inizializzo vettore per storare valori delta velocità per ogni om_tn
delta_velocita = zeros(length(om_tn), 1);


for j=1:length(om_tn)
    om_t=om_tn(j);

% a) Definizione della matrice da sistema di riferimento inerziale a perifocale:
    R3_om = [cos(om_t), sin(om_t), 0; 
            -sin(om_t), cos(om_t), 0; 
            0, 0, 1];

    T = R3_om * R1_i * R3_OM;

%  b) Calcolo di posizione nel punto 1 e punto 2 nel sistema perifocale dell’orbita di trasferimento:
    r1_t=cross(T,r1);
    r2_t=cross(T,r2);

%  c) Determinazione delle anomalie vere sull’orbita di trasferimento:    

co_1tn = r1_t(1)/norm(r1_t);
si_1tn = r1_t(2)/norm(r1_t);
theta_1t = atan2(si_1tn, co_1tn);

co_2tn = r2_t(1)/norm(r2_t);
si_2tn = r2_t(2)/norm(r2_t);
theta_2t = atan2(si_2tn, co_2tn);   

% d) Determinazione dei parametri di forma (semiasse maggiore e eccentricità) tramite l’equazione della conica:

e_t = (r1-r2)/(r1*cos(theta_1t) - r2*cos(theta_2t));
a_t = r1*(1+e_t*cos(theta_1t))/(1-e_t^2);

% e) Determinazione della velocità nei due punti di manovra per ogni 𝜔𝑇,𝑛 :

[~,v1t] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i, mu)

[~,v2t] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu)


%f)Calcolo del Δ𝑣 per ogni 𝜔𝑇,𝑛
delta_v1 = norm(v1t - v1i);
delta_v2 = norm(v2f - v2t);

delta_velocita(j) = delta_v1 + delta_v2; % Salvo il Δ𝑣 totale per ogni 𝜔𝑇,𝑛

end

%8) Identificazione del valore di 𝜔𝑇,𝑛 che minimizza il Δ𝑣 totale (Δ𝑣1 + Δ𝑣2):

[dv_best, idx_best] = min(delta_velocita);
om_t_best = om_tn(idx_best);
