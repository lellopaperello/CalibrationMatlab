function vt = reference(D, model)
% Implementation of the terminal velocity - diameter relations from the
% literature
switch model
    case 'Brandes - 2008, T = -1°C'
        % Diameter in mm
        D = D*1e3;
        vt = 0.87*D.^0.23;
    case 'Brandes - 2008, T = -5°C'
        % Diameter in mm
        D = D*1e3;
        vt = 0.67*D.^0.25;
    case 'Brandes - 2008, T = -10°C'
        % Diameter in mm
        D = D*1e3;
        vt = 0.55*D.^0.23;
    case 'Kuhn-Martin - 2019, graupel'
        % Diameter in micron
        D = D*1e6;
        vt = 0.0013*D.^0.98;
        vt(D > 1200) = NaN;
    case 'Kuhn-Martin - 2019, rimed needles'
        % Diameter in micron
        D = D*1e6;
        vt = 0.02*D.^0.41;
        vt(D > 2000) = NaN;
    case 'Justo-Bosworth - 1971, dendrites'
        % Radius in cm
        r = D*50;
        vt = 1.23*r.^0.2;
    case 'Justo-Bosworth - 1971, plates & columns'
        % Radius in cm
        r = D*50;
        vt = 1.78*r.^0.2;
    case 'Justo-Bosworth - 1971, average'
        % Radius in cm
        r = D*50;
        vt = 1.50*r.^0.2;
end