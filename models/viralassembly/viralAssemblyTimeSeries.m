%% INITIALIZATION

% Initial concentrations
init = zeros(17,1);

num_pure_founders   = 7;
num_impure_founders = 0;

vol = 1;

init(1) = num_pure_founders/vol;
init(3) = num_impure_founders/vol;

% Time to end simulation
Tend = 60*10;

% ODE options
vaODE = @viralAssemblyODE;
options = odeset('Refine', 6);

for ii = 0:1
    % DNA parameters
    k_v         = 0.1;     % Rate of virions reentering host nucleus
    alpha_c     = 0.1;     % Rate of repair of gapped DNA (Nakabayashi)
    DNA_degrade = 0.001;   % DNA degradation rate (Nakabayashi)
    pure2mod    = 0.01*ii; % Rate of Cas9 modification of DNA
    D_max       = 100/vol; % Maximum genomes allowed in nucleus

    % RNA parameters
    alpha_19 = 0.01;  % Transcription rate of 19S
    alpha_35 = 0.09;  % Transcription rate of 35S
    gamma_19 = 0.001; % Degradation rate of 19S
    gamma_35 = 0.001; % Degradation rate of 35S
    frac_u   = 0.3;   % Fraction of unspliced 35S RNA

    % Protein parameters
    beta          = 0.1;         % Translation rate for P1-P5
    beta_6        = 0.1;         % Translation rate for P6
    delta         = 0.001;       % Degradation rate for P1-P5
    delta_6       = 0.001;       % Degradation rate for P6
    p6_activation = (9.3-1.6)/2; % Half-sat. constant for P6 transactivation

    % Virion parameters
    delta_v  = 0.001; % Degradation rate of virions
    v_exit   = 0.1;   % Rate at which virions exit the cell
    k_p      = 0.1;   % Packaging rate
    k_anchor = 0.1;   % Rate of P3 binding to virions

    % RNAi parameters
    L       = 0; % 0
    k_value = 0; % 0
    x0      = 0; % 0

    parameters = [k_v alpha_c DNA_degrade pure2mod D_max ...
                  alpha_19 alpha_35 gamma_19 gamma_35 frac_u ...
                  beta beta_6 delta delta_6 p6_activation ...
                  delta_v v_exit k_p k_anchor ...
                  L k_value x0];

    % Simulation
    tic
    [t, out] = ode23s(vaODE, [0 Tend], init, options, parameters);
    toc

    % Plot!
    figure(1);
    subplot(1,2,ii+1)
    plot(t, out(:,16), t, out(:,17), t, out(:,16)+out(:,17),'LineWidth',2);
    title('Concentration of Virions over Time')
    xlabel('Time (minutes)')
    ylabel('Virions')
    legend('Pure virions','Impure virions','Total Virions')
    axis([0 600 -0.1 1.2])
end