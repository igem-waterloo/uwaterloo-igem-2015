function dS = viralAssemblyODE(~,S,params)

	%% STATE VARIABLES
    
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
	P_5    = S(12); % P5
	P_6    = S(13); % P6
	V_i    = S(14); % Intermediate pure virions
	V_mi   = S(15); % Intermediate impure virions
	V      = S(16); % Pure virions
	V_m    = S(17); % Impure virions
    
    %% PASSING PARAMETERS
    
    % DNA parameters
    k_v         = params(1); % Rate of virions reentering host nucleus
    alpha_c     = params(2); % Rate of repair of gapped DNA
    DNA_degrade = params(3); % DNA degradation rate
    pure2mod    = params(4); % Rate of Cas9 modification of DNA
    D_max       = params(5); % Maximum genomes allowed in nucleus
    
    % RNA parameters
    alpha_19 = params(6);  % Transcription rate of 19S
    alpha_35 = params(7);  % Transcription rate of 35S
    gamma_19 = params(8);  % Degradation rate of 19S
    gamma_35 = params(9);  % Degradation rate of 35S
    frac_u   = params(10); % Fraction of unspliced 35S RNA
    
    % Protein parameters
    beta          = params(11); % Translation rate for P1-P5
    beta_6        = params(12); % Translation rate for P6
    delta         = params(13); % Degradation rate for P1-P5
    delta_6       = params(14); % Degradation rate for P6
    p6_activation = params(15); % MM coefficient for P6 transactivation
    
    % Virion parameters
    delta_v  = params(16); % Degradation rate of virions
    v_exit   = params(17); % Rate at which virions exit the cell
    k_p      = params(18); % Packaging rate
    k_anchor = params(19); % Rate of P3 binding to virions
    
    % RNAi parameters
    L       = params(20); % 0
    k_value = params(21); % 0
    x0      = params(22); % 0

    
    %% DERIVED PARAMETERS
	frac_s = 1-frac_u;     % Fraction of spliced 35S RNA
	R_35u  = frac_u*R_35;  % Unspliced 35S
	R_35s  = frac_s*R_35;  % Spliced 35S
	R_35mu = frac_u*R_35m; % Modified unspliced 35S
	R_35ms = frac_s*R_35m; % Modified spliced 35S
    
    pure2modg = pure2mod;
    pure2modc = pure2mod;
    
	beta_1 = beta;
	beta_2 = beta;
	beta_3 = beta;
	beta_4 = beta;
	beta_5 = beta;
    
	delta_1 = delta;
	delta_2 = delta;
	delta_3 = delta;
	delta_4 = delta;
	delta_5 = delta;
    
	RNAi_factor = L / (1 + exp(-k_value * (P_6 - 0.5)));
    D_total     = D_gap + D_gmod + D_ccc + D_cmod;
    
	%% DIFFERENTIAL EQUATIONS

	% Gapped DNA in nucleus
	eq1 = k_v*V*(D_max-D_total) - DNA_degrade*D_gap - alpha_c*D_gap - pure2modg*D_gap;

	% Covalently closed circular DNA in nucleus
	eq2 = alpha_c*D_gap - DNA_degrade*D_ccc - pure2modc*D_ccc;

	% Modified gapped DNA in nucleus
	eq3 = k_v*V_m*(D_max-D_total) - DNA_degrade*D_gmod - alpha_c*D_gmod + pure2modg*D_gap;

	% Modified cccDNA in nucleus
	eq4 = alpha_c*D_gmod - DNA_degrade*D_cmod + pure2modc*D_ccc;

	% 19S RNA
	eq5 = alpha_19*D_ccc - (gamma_19+RNAi_factor)*R_19s;

	% Total pure 35S RNA (spliced or unspliced)
	eq6 = alpha_35*D_ccc - (gamma_35+RNAi_factor)*R_35 - k_p*P_4*P_5*R_35u;

	% Total impure 35S RNA (spliced or unspliced)
	eq7 = alpha_35*D_cmod - (gamma_35+RNAi_factor)*R_35m - k_p*P_4*P_5*R_35mu;

	% Protein 1
	eq8 = beta_1 * P_6/(P_6 + p6_activation) * (R_35u + R_35mu) - delta_1*P_1;

	% Protein 2
	eq9 = beta_2 * P_6/(P_6 + p6_activation) * (R_35u + R_35mu) - delta_2*P_2;

	% Protein 3
	eq10 = beta_3 * P_6/(P_6 + p6_activation) * (R_35 + R_35m) - delta_3*P_3 - k_anchor*P_3*(V_i+V_mi);

	% Protein 4
	eq11 = beta_4 * P_6/(P_6 + p6_activation) * (R_35 + R_35m) - delta_4*P_4 - k_p*P_4*P_5*(R_35u+R_35mu);
    
	% Protein 5
	eq12 = beta_5 * P_6/(P_6 + p6_activation) * (R_35 + R_35m) - delta_5*P_5 - k_p*P_4*P_5*(R_35u+R_35mu);

	% Protein 6
	eq13 = beta_6*R_19s - delta_6*P_6;

	% Intermediate pure virions
	eq14 = k_p*P_4*P_5*R_35u - k_anchor*P_3*V_i;

	% Intermediate impure virions
	eq15 = k_p*P_4*P_5*R_35mu - k_anchor*P_3*V_mi;

	% Complete pure virions
	eq16 = k_anchor*P_3*V_i - k_v*V*(D_max-D_total) - delta_v*V - v_exit*V;

	% Complete impure virions
	eq17 = k_anchor*P_3*V_mi - k_v*V_m*(D_max-D_total) - delta_v*V_m - v_exit*V_m;

	% Altogether now!
	dS = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12 eq13 eq14 eq15 eq16 eq17]';    

end
