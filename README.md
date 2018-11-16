# MCML_MODEL

## This code is used for MCML simulation
### The code is written in MATLAB OOP format, the whole result could be calculated and showed by methods function
#### Basic Parameters
* iterations = 10000; % iteration times
* ni=1;               % reflection index of input
* nt=1.37;            % reflection index of output
* slice = 200;        % the size of total dividence
* threshold_dead = 0.00001; % the dead threshold for photon package
* prob_survive = 0.01; % the possibility of survive
* x=0;                % position of photon package along x
* y=0;                % position of photon package along y
* z=0;                % position of photon package along z(depth direction)
* vx=0;               % direction of photon package along x
* vy=0;               % direction of photon package along y
* vz=1;               % direction of photon package along z(depth direction)
* weight=1;           % weight of each photon package
* out_bool=0;         % whether the photon package is outside
* dead_bool=0;        % the survive condition of photon package

#### Basic Methods
* MCML_model      % Construction function
* initial_package % Package initialization
* scatter         % photon scatter
* move            % movie the photon
* is_hit_boundary % check hit or not
* hit_boundary    % process hit boundary
* absorption_record % record the absorption
* check_dead      % check dead or not
* simulate_once   % simulate once
* simulate        % simulate iterations time
* imshow_absorption % show the absorption distribution in x,y,z and calculate the decay rate

### free copyright by Zheyuan Zhang
