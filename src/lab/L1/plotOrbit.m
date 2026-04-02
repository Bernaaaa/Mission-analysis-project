function plotOrbit(a, e, i, Om, om, th0, thf, dth, mu)
    %La funzione plotta l'orbita a partire dai parametri orbitali
    %
    % plotOrbit(a, e, i, Om, om, th0, thf, dth, mu)
    %
    % -----------------------------------------------------------------
    % INPUT:
    % a - semiasse maggiore in [km]
    % e - eccentricità
    % i - inclinazione in [rad]
    % Om - RAAN in [rad]
    % om - argomento del pericentro in [rad]
    % th0 - anomalia vera iniziale in [rad]
    % thf - anomalia vera finale in [rad]
    % dth - passo di integrazione in [rad]
    % mu - parametro gravitazionale standard in [km^3/s^2]

% Define the range of true anomaly
th_value = th0: dth: thf;
r_value = [];
v_value = [];
for th = th_value
    [rr, vv] = par2car(a, e, i, Om, om, th, mu);
    r_value = [r_value, rr];
    v_value = [v_value, vv];
end
plot3(r_value(1,:), r_value(2,:), r_value(3,:))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Orbit Plot')
grid on
axis equal
hold on
plot3(0, 0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
legend('Orbit', 'Central Body')
hold off

end 


    