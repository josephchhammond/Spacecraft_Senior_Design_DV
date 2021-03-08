function results2(PercentCoverage,preposition_DV1,preposition_DV2,V,DV1,DV2,m_array,R2,m_break)



%finds indecies of the best percent coverage
[vals,i1] = max(PercentCoverage);
[maxval,i2] = max(vals);
ii = i1(i2);
jj = i2;

if size(PercentCoverage,1) == 1
   jj = ii;
   ii = 1;
end

% Find index ii and jj of best system
R2_ii = R2(ii);
dv1 = DV1(ii,jj);
dv2 = DV2(ii,jj);
v_total = V(ii,jj);
mass_sys = [m_array(:,ii,jj)];



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
fprintf('Total Tank Volume: %.2f m^3\n', v_total);
disp('~~~~~');


disp(' ');
disp('Mass Breakdown:');
fprintf('Payload mass: %.2f kg\n', mass_sys(1));
disp(' ');
fprintf('Preposition Chem Mass: %.2f kg\n', mass_sys(9));
fprintf('Departure Eprop Mass: %.2f kg\n', mass_sys(2));
fprintf('Arrival Chem Mass: %.2f kg\n', mass_sys(3) - mass_sys(9));
fprintf('Mass Staged After Departure Burn: %.2f kg\n', mass_sys(10));
disp(' ');
fprintf('EProp Tank Mass: %.2f kg\n', mass_sys(4));
fprintf('Chem Prop Tank Mass: %.2f kg\n', mass_sys(5));
fprintf('Additional Power Mass (external to payload): %.2f kg\n', mass_sys(6));
fprintf('Thrust Plate Mass : %.2f kg\n', mass_sys(7) - mass_sys(6) - mass_sys(5) - mass_sys(4));
fprintf('Total Prop Dry Mass: %.2f kg\n', mass_sys(7));
disp(' ');
fprintf('Total Mass: %.2f kg\n', mass_sys(8));
disp('~~~~~');



if size(PercentCoverage,1) == 1
    figure
    plot(m_break,PercentCoverage);
    title('For optimization use, not for external')
    xlabel('Mass Breakdown')
    ylabel('ISO Coverage (%)')
elseif size(PercentCoverage,2) == 1
    figure
    plot(R2,PercentCoverage);
    title('For optimization use, not for external')
    xlabel('Design R2 (AU)')
    ylabel('ISO Coverage (%)')
else
    figure
    surf(m_break,R2,PercentCoverage);
    title('For optimization use, not for external')
    xlabel('Mass Breakdown')
    ylabel('Design R2 (AU)')
    zlabel('ISO Coverage (%)')
    colorbar
    caxis([0,16]) % fixes colorbar colors to cap at 15, gives consistent coloring regardless of plot
end

