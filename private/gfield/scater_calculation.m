%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of photon absorption rate

% Start with Suva's measurement
Power_measurement = 0.5e-6; % W
spot_area = 4e-6; % m^2


% Set physical constants
SpeedOfLight = 3e8; % m/s
PlanksConstant = 6.626e-34; % Js

% set wavelength range (note range is constrated by Baylor nomogram, below)
WaveLengthRange = [778:-2:380] * 1e-9; % in m 

% estimate the laser spectra as a delta function
laser_lambda = 488e-9; % m
laser_spectra = zeros(length(WaveLengthRange),1);
laser_spectra(WaveLengthRange == laser_lambda) = 1;

% M-opsin spectral sensitivty from with Baylor nomogram
WL = 778:-2:380;
WL = WL.*0.001;
a0 = -5.2734;
a1 = -87.403;
a2 = 1228.4;
a3 = -3346.3;
a4 = -5070.3;
a5 = 30881;
a6 = -31607;
lambdaMax = 500;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
RodPhotonSensitivity = 10.^LogPhotonSensitivity;

% plot laser spectra (blue) and cone spectra (black)
plot(WaveLengthRange, RodPhotonSensitivity, 'k', WaveLengthRange, laser_spectra, 'b')
xlabel('wavelength [nm]')
ylabel('relative power (sensitivity)')
legend('M cone sensitivity', 'laser spectra')

%%%%%%%%
% actualy calculation begins here:

% scale laser spectra by the measured power;
scaled_laser_spectra = laser_spectra .* Power_measurement;

Intensity = (scaled_laser_spectra .* WaveLengthRange') ./ (PlanksConstant * SpeedOfLight);
PhotonFlux = Intensity ./ spot_area;

RodCollectingArea = 1e-12; %(Baylor 1984, Table 1) m^2 (photoreceptor collecting areay is ~1µm^2)
EffectivePhotonFlux = dot(PhotonFlux,RodPhotonSensitivity([200:-1:1])');
PhotonCatchRate = EffectivePhotonFlux * RodCollectingArea

