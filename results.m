function results(PercentCoverage,preposition_DV1,preposition_DV2,DV1,DV2,mass_array1,R2,m_break,mass_payload, power_payload, prop_scheme, R1)




%finds indecies of the best percent coverage
[vals,i1] = max(PercentCoverage);
[maxval,i2] = max(vals);
ii = i1(i2);
jj = i2;

% Find index ii and jj of best system
R2_ii = R2(ii);
dv1 = DV1(ii,jj);
dv2 = DV2(ii,jj);


m3 = mass_array1(1);

[~,power_mass4,power_area4] = panel_power(R2_ii, [], power_payload); 

[mass_array2,power_area2,~] = prop_sizing1(m3, power_area4, R1, dv1, prop_scheme(2,:));
m2 = mass_array2(1);

[mass_array3,power_area3,~] = prop_sizing1(m2, power_area4, R2_ii, dv2, prop_scheme(3,:));
m4 = mass_array3(1);

power_area1 = power_area2;


disp(' ');
disp('Optimal Solution:');
fprintf( 'ISO Coverage: %.2f percent\n', PercentCoverage(ii,jj));
fprintf( 'Maximum Intercept Range: %.2f AU\n', R2_ii);
disp('~~~~~');

disp(' ');
disp('Propulsive Capabilities:');
fprintf('Earth Hyberbolic Excess: %.2f m/s\n', preposition_DV1);
fprintf('Preposition DV: %.2f m/s\n', preposition_DV2);
fprintf('Departure DV: %.2f m/s\n', dv1);
fprintf('Arrival DV: %.2f m/s\n', dv2);
disp('~~~~~');


disp(' ');
disp('Stage 1 - Preposition');
disp(' ');
fprintf('Payload mass: %.2f kg\n', mass_array1(1));
fprintf('Propellant mass: %.2f kg\n', mass_array1(2));
fprintf('Propellant structure mass: %.2f kg\n', mass_array1(3));
fprintf('Power mass: %.2f kg\n', mass_array1(4));
fprintf('Dry mass: %.2f kg\n', mass_array1(5));
fprintf('Total stage mass: %.2f kg\n', mass_array1(6));
fprintf('Total solar array area required: %.2f m2\n', power_area1);
disp('~~~~~');

disp(' ');
disp('Stage 2 - Departure');
disp(' ');
fprintf('Payload mass: %.2f kg\n', mass_array2(1));
fprintf('Propellant mass: %.2f kg\n', mass_array2(2));
fprintf('Propellant structure mass: %.2f kg\n', mass_array2(3));
fprintf('Power mass: %.2f kg\n', mass_array2(4));
fprintf('Dry mass: %.2f kg\n', mass_array2(5));
fprintf('Total stage mass: %.2f kg\n', mass_array2(6));
fprintf('Total solar array area required: %.2f m2\n', power_area2);
disp('~~~~~');

disp(' ');
disp('Stage 3 - Arrival');
disp(' ');
fprintf('Payload mass: %.2f kg\n', mass_array3(1));
fprintf('Propellant mass: %.2f kg\n', mass_array3(2));
fprintf('Propellant structure mass: %.2f kg\n', mass_array3(3));
fprintf('Power mass: %.2f kg\n', mass_array3(4));
fprintf('Dry mass: %.2f kg\n', mass_array3(5));
fprintf('Total stage mass: %.2f kg\n', mass_array3(6));
fprintf('Total solar array area required: %.2f m2\n', power_area3);
disp('~~~~~');

disp(' ');
disp('Payload');
disp(' ');
fprintf('Payload mass: %.2f kg\n', mass_payload);
fprintf('Power mass: %.2f kg\n', power_mass4);
fprintf('Total stage mass: %.2f kg\n', m4);
fprintf('Total solar array area required: %.2f m2\n', power_area4);
disp('~~~~~');


figure
surf(m_break,R2,PercentCoverage);
title('For optimization use, not for external')
xlabel('Mass Breakdown')
ylabel('Design R2 (AU)')
zlabel('ISO Coverage (%)')
colorbar
caxis([0,16]) % fixes colorbar colors to cap at 15, gives consistent coloring regardless of plot
end

