function viralAssemblyModel ()

	%% Initial Conditions
	% We're setting everything to 0 except the first variable
	% The first variable is the copy number of gapped DNA in the nucleus
	% This simulates the initial infection 
	initial = zeros(18,1);
	initial(1) = 0.0001;
	vol = 75; % Volume of cell

	%% Run the simulation
	% Run for 100 time units (minutes?)
	Tend = 60*24*7;
	% I use tic and toc to time how long the simulation takes
	tic
	[t, out] = ode15s(@mathva, [0,Tend], initial);
	toc

	%% Plot!
	% This plots concentration of virions over time
	figure(1);
	plot(t, out(:,1), t, out(:,2), t, out(:,17),'LineWidth',2);
	title('Concentration of Virions over Time')
	xlabel('Time')
	ylabel('Virions')

end

function dS = mathva(~,S)

	%% Provide useful names for things
	D_gap  = S(1);  % Gapped DNA in nucleus
	D_ccc  = S(2);  % cccDNA in nucleus
	D_gmod = S(3);  % Modified gapped DNA in nucleus, unable to produce 19S
	D_cmod = S(4);  % Modified cccDNA in nucleus, unable to produce 19S
	R_19s  = S(5);  % 19S RNA
	R_35   = S(6);  % Unspliced, unmodified 35S RNA
	R_35m  = S(7);  % Modified 35S RNA
	P_1    = S(8);  % P1
	P_2    = S(9);  % P2
	P_3    = S(10); % P3
	P_4    = S(11); % P4
	P_4s   = S(12); % P4 subunits
	P_5    = S(13); % P5
	P_6    = S(14); % P6
	V_i    = S(15); % Intermediate pure virions
	V_mi   = S(16); % Intermediate impure virions
	V      = S(17); % Pure virions
	V_m    = S(18); % Impure virions

	frac_u = 0.3;          % Fraction of unspliced 35S RNA
	frac_s = 1-frac_u;     % Fraction of spliced 35S RNA
	R_35u  = frac_u*R_35;  % Unspliced 35S
	R_35s  = frac_s*R_35;  % Spliced 35S
	R_35mu = frac_u*R_35m; % Modified unspliced 35S
	R_35ms = frac_s*R_35m; % Modified spliced 35S

	%% Constants
	% ALL TIME (SHOULD BE) IN MINUTES
	out_inf_pure  = 0;          % No outside infection (for now)
	out_inf_mod  = 0;          % No outside infection (for now)
	k_v      = 0.01;       % Rate of virions reentering host nucleus
	alpha_c  = 0.1;        % Rate of repair of gapped DNA, value from Nakabayashi
	DNA_degrade  = 0.0001;    % Degradation rate of CCC DNA

	pure2modg = 0.01 %rate in which pure 35s gapped is changed to modified 35s gapped
	pure2modc = 0.01 %rate in which pure 35s ccc is changed to modified 35s ccc

	alpha_19 = 0.05;       % Transcription rate of 19S, value from Nakabayashi (tunable)
	gamma_19 = log(2)/600; % Degradation rate = ln(2)/(half-life)
	alpha_35 = 0.0653;     % Transcription rate of 35S, value from Martiene
	gamma_35 = log(2)/600; % Degradation rate = ln(2)/(half-life)

	beta_1  = 0.1; % From Nakabayashi
	beta_2  = 0.1; % From Nakabayashi
	beta_3  = 0.1; % From Nakabayashi
	beta_4  = 0.1; % From Nakabayashi
	beta_5  = 0.1; % From Nakabayashi
	beta_6  = 0.1; % From Nakabayashi

	delta_1 = 0.001;   % From Nakabayashi
	delta_2 = 0.001;   % From Nakabayashi
	delta_3 = 0.001;   % From Nakabayashi
	delta_4 = 0.001;   % From Nakabayashi
	delta_5 = 0.001;   % From Nakabayashi
	delta_6 = 0.001;   % From Nakabayashi

	delta_v = 0.001;    % Rate of degradation of virions
	v_exit  = log(2)/100;    % Rate that virions leave the cell

	k_p   = 0.1;    % Packaging rate
	k_l   = 1;      % Rate of P2 leaving for elIB
	k_p5s = 1;      % Splicing of P4 (tunable)
	k_anchor = 1;   % Rate of P3 binding to virions

	D_max = 0.5;   % Maximum number of gapped DNA in nucleus

	L = 1; % made up parameters
	k_value = 0.01; % made up parameters
	p6_activation = (9.3/1.6)/10;

	RNAiFactor = L / (1 + exp(-k_value * (P_6 - 0.5)));

	%% Equations

	% Gapped DNA in nucleus
	eq1 = out_inf_pure + k_v*V*(D_max-D_gap-D_gmod-D_cmod-D_ccc) - alpha_c*D_gap - pure2modg*D_gap - DNA_degrade * D_gap ;

	% Covalently closed circular DNA in nucleus
	eq2 = alpha_c*D_gap - DNA_degrade*D_ccc - pure2modc*D_ccc;

	% Modified gapped DNA in nucleus
	eq3 = out_inf_mod + k_v*V_m*(D_max-D_gap-D_gmod-D_cmod-D_ccc) - alpha_c*D_gmod + pure2modg*D_gap - DNA_degrade*D_gmod; % Brandon to update

	% Modified cccDNA in nucleus
	eq4 = alpha_c*D_gmod - DNA_degrade*D_cmod + pure2modc*D_ccc; % Brandon to update

	% 19S RNA
	eq5 = alpha_19*D_ccc - gamma_19*R_19s - RNAiFactor;

	%when multiplied by gamma 19 gets a weird result in DE

	% Total pure 35S RNA (spliced or unspliced)
	eq6 = alpha_35*D_ccc - gamma_35*R_35 - k_p*P_4s*P_5*R_35u - RNAiFactor;

	% Total impure 35S RNA (spliced or unspliced)
	eq7 = alpha_35*D_cmod - gamma_35*R_35m - k_p*P_4s*P_5*R_35mu - RNAiFactor; % Brandon to update

	% Protein 1
	eq8 = beta_1*(1 + P_6/ (P_6 + p6_activation)) * (R_35u + R_35mu) - delta_1 * P_1;

	% Protein 2
	eq9 = beta_2 * (1 + P_6/ (P_6 + p6_activation)) * (R_35u + R_35mu) - delta_2*P_2 - k_l*P_2;

	% Protein 3
	eq10 = beta_3 * (1 + P_6/ (P_6 + p6_activation)) * (R_35 + R_35m) - delta_3*P_3 - k_anchor*P_3*V_i - k_anchor*P_3*V_mi;

	% Protein 4
	eq11 = beta_4 * (1 + P_6/ (P_6 + p6_activation)) * (R_35 + R_35m) - delta_4*P_4 - k_p5s*P_4;

	% Protein 4, subunits
	eq12 = k_p5s * P_4 - delta_4 * P_4s - k_p*P_4s*P_5*R_35u - k_p*P_4s*P_5*R_35mu;

	% Protein 5
	eq13 = beta_5 * (1 + P_6/ (P_6 + p6_activation)) * R_35 - delta_5*P_5 - k_p*P_4s*P_5*R_35u - k_p*P_4s*P_5*R_35mu;

	% Protein 6
	eq14 = beta_6*R_19s - delta_6*P_6;

	% Intermediate pure virions
	eq15 = k_p*P_4s*P_5*R_35u - k_anchor*P_3*V_i;

	% Intermediate impure virions
	eq16 = k_p*P_4s*P_5*R_35mu - k_anchor*P_3*V_mi; % Brandon to update

	% Complete pure virions
	eq17 = k_anchor*P_3*V_i - k_v*V*(D_max-D_gap-D_gmod-D_cmod-D_ccc) - delta_v*V - v_exit*V;

	% Complete impure virions
	eq18 = k_anchor*P_3*V_mi - k_v*V_m*(D_max-D_gap-D_gmod-D_cmod-D_ccc) - delta_v*V_m - v_exit*V_m; % Brandon to update

	% Altogether now!
	dS = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12 eq13 eq14 eq15 eq16 eq17 eq18]';    

end