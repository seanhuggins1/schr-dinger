% function seq_crank(Nspace,Ntime,tau,varargin)
% Monolithic function which integrates the time dependent Schrodinger Equation using a
% Crank-Nicholson or FTCS scheme, and produces an animated plot of the time
% evolution.
%
% Required Input:
% ===============
%
% Nspace (:) The number of spatial grid points
% Ntime (:) The number of time steps to perform the integration for
% tau (:) The time step for each forward step of the wave function.
%
% Output:
% =======
%
% No output
%
%
% Varargin:
% =========
%
% 'Method'	['Implicit'] set to 'Implicit' for the implicit Crank-Nicholson
%                      scheme. The Crank-Nicholson scheme is unconditionally stable for all tau
%                      set to 'explicit' for the explicit
%                      forward-time-centered-space scheme. 
%                      The Crank-Nicholson scheme is unconditionally unstable for all tau
%
% 'Length'	[200] changes the size of the spatial domain given by -L/2 < x < L/2
%
% 'PlotType'    ['psi'] set to 'psi' to plot an animation of the time
%                      evolution of the wave function psi(x,t) given by the time dependent
%                      Schrodinger Equation
%                      set to 'prob' to plot an animation of the time evolotion of the probability
%                      density given by psi(x,t)*conj(psi(x,t))
%                      setting to 'prob' will also produce a plot of total
%                      probability conservation
%                      set to 'party' to simply have a good time
%
% 'Potential'	['[]'] pass in spatial indicies j, where the user wishes to
%                      set V(x) = 1. The spatial indicies j can be between 1 < j < Nspace
%
% 'WParam'	[10 0 0.5] parameters can be passed in to change the initial
%                      condition of the wave function. WParam is a three element array corresponding to
%                      [sigma0 x0 k0].
%
% Requires: no external m-files
% =========
%
% Example Use: seq_crank(300, 250, 0.1, 'Method', 'Implicit', 'Length', 300, 'PlotType', 'psi', 'Potential', [50 70:200 240], 'Wparam', [15 -10 0.4]);
% ============
%
% Author:
% =======
%
% SHuggins 16 Dec. 2018
%
function seq_crank(Nspace,Ntime,tau,varargin)

% If called with no arguments, echo a useage line. 
if nargin == 0
    disp(' ')
    disp('ERROR(seq_crank): usage is seq_crank(Nspace,Ntime,tau,varargin)')
    disp(' ')
    return
end

%Check validity of Nspace input
if ~isnumeric(Nspace) || ~isreal(Nspace) || isinf(Nspace) || isnan(Nspace) || Nspace <= 0 || (floor(Nspace)~= Nspace)
    disp(' ')
    disp('ERROR(seq_crank): The number of spatial grid points must be a finite, positive integer greater than zero')
    disp('Choosing default number of spatial grid points: 300')
    disp(' ')
    Nspace = 300;
end

%Check validity of Ntime input
if ~isnumeric(Ntime) || ~isreal(Ntime) || isinf(Ntime) || isnan(Ntime) || Ntime <= 0 || (floor(Ntime)~= Ntime)
    disp(' ')
    disp('ERROR(seq_crank): The number of iterations must be a finite, positive integer greater than zero')
    disp('Choosing default number of iterations: 500')
    disp(' ')
    Ntime = 500;
end

%Check validity of tau input
if ~isnumeric(tau) || ~isreal(tau) || isinf(tau) || isnan(tau) || tau <= 0
    disp(' ')
    disp('ERROR(seq_crank): The time step must be a finite, positive real number greater than zero')
    disp('Choosing default tau: 0.5')
    disp(' ')
    tau = 0.5;
end

% Check that all varargin come in pairs.
if mod(length(varargin),2) ~= 0
  disp(' ')
  disp('ERROR(seq_crank): mis-match (odd number) of vargargin inputs')
  disp(' ')
  return
end

% Set defaults
method='implicit';
L=200;
plotType='psi';
V=zeros(Nspace, 1);
potential_spikes = [];
WParam=[10 0 0.5];

% Parse the varargin arguments
for j=1:2:length(varargin)
    switch lower(varargin{j})
    case 'potential'
       potential_spikes = round(varargin{j+1});
       
       %Check validity of potential spikes input
       if ~isnumeric(potential_spikes) || ~isreal(potential_spikes) ...
               || any(isinf(potential_spikes)) || any(isnan(potential_spikes)) ...
               || any(potential_spikes <= 0)  || any(potential_spikes > Nspace) || any(floor(potential_spikes)~= potential_spikes)
           disp(' ')
           disp('ERROR(seq_crank): indicies for potential spikes must be finite, positive integers within the spatial grid.');  
           disp('Choosing default potential spike indicies: []'); 
           disp(' ')
           potential_spikes = [];
       end
    case 'length'
       L=varargin{j+1};
       
       %Check validity of length input
       if ~isnumeric(L) || ~isreal(L) || isinf(L) || isnan(L) || L <= 0
           disp(' ')
           disp('ERROR(seq_crank): length must be a finite, positive, real number.'); 
           disp('Choosing default length: 200'); 
           disp(' ')
           L=200;
       end  
    case 'method'
       method=lower(varargin{j+1});
           
       %Allowed strings for method type
       allowedMethods = {'implicit','explicit'};
       
       %Check validity of method type input
       if ~any(strcmp(method,allowedMethods)) 
           disp(' ')
           disp('ERROR(seq_crank): method type not recognized. Please choose "implicit" for Crank-Nicholson or "explicit" for FTCS.');      
           disp('Choosing default method type: "implicit"');
           disp(' ')
           method='implicit';
       end        
    case 'plottype'
       plotType=lower(varargin{j+1});
       
       %Allowed strings for plot type
       allowedPlots = {'psi','prob','party'};
       
       %Check validity of plot type input
       if ~any(strcmp(plotType,allowedPlots))   
           disp(' ')
           disp('ERROR(seq_crank): plot type not recognized. Please choose "psi" to plot the wave function or "prob" to plot probability density and total probability.');      
           disp('the user may also choose "party" to simply have a good time.');
           disp('Choosing default plot type: "psi"');
           disp(' ')
           plotType='psi';
       end       
       
    case 'wparam'
       WParam = varargin{j+1}; 
      
       %Check validity of parameters input
       if length(WParam) ~= 3 
           disp(' ')
           disp('ERROR(seq_crank): WParam must be a vector with 3 parameters [sigma0 x0 k0].');
           disp('Choosing default WParam: [10 0 0.5]');
           disp(' ')
           WParam=[10 0 0.5];
       end
       if ~isnumeric(WParam) || ~isreal(WParam) ...
               || any(isinf(WParam)) || any(isnan(WParam))
           disp(' ')
           disp('ERROR(seq_crank): parameters in WParam [sigma0 x0 k0] must be finite, real numbers.');
           disp('Choosing default WParam: [10 0 0.5]');
           disp(' ')
           WParam=[10 0 0.5];
       end
       %Make sure that x0 is within the domain
       if (WParam(2) >= L/2) || (WParam(2) <= -L/2)
           disp(' ')
           disp('ERROR(seq_crank): x0 must be within the domain -L/2 < x < L/2');
           disp('Choosing default WParam: [10 0 0.5]');
           disp(' ')
           WParam=[10 0 0.5];
       end
       
    otherwise
       disp(' ')
       disp(sprintf('WARNING: unknown varargin <%s> ignored',varargin{j}))
       disp(' ')
    end
end

%Place spikes of potential in potential array
V(potential_spikes,1) = 1;

%Grab WParam constants from inputs
sig0 = WParam(1);
x0 = WParam(2);
k0 = WParam(3);

%Set important constants
hbar = 1;
m = 1/2;

%Set the spatial grid spacing 
h = L/(Nspace-1);



%Compute the Hamiltonian Matrix
I = eye(Nspace);
for j = 1:Nspace
    %Logical indexing to handle periodic boundary conditions
    jm = j-1;
    jp = j+1;
    %If our left end is out of range, set it to the right end
    if jm == 0
        jm = Nspace;
    end
    %If our right end is out of range, set it to the left end
    if jp == Nspace+1
        jp = 1;
    end
    
    %Build the tridiagonal Hamiltonian with periodic boundary conditions
    for k = 1:Nspace
        H(j,k) = ((-hbar^2)/m)*(I(jp,k) + I(jm,k) - 2*I(j,k))/(h^2) + V(j,1)*I(j,k);
    end
end


%Crank-Nicholson Scheme
if strcmp(method,'implicit')
    %Calculate the matrix needed for forward time integration
    M=(inv((I + ((1i*tau)/(2*hbar)).*H)))*(I - ((1i*tau)/(2*hbar)).*H);
    % NO STABILITY CHECK REQUIRED. CRANK NICHOLSON IS UNCONDITIONALLY
    % STABLE
    
%FTCS Scheme
elseif strcmp(method,'explicit')
    %Calculate the matrix needed for forward time integration
    M = (I - ((1i*tau)/(hbar)).*H);
    %Grab spectral radius
    r = max(abs(eig(M)));
    %Check for stability
    if r > 1.000001
        disp(' ')
        disp(sprintf('ERROR(seq_crank): Spectral radius: %.7f is greater than unity, solution will be unstable.', r)); 
        disp('To run FTCS, please choose a time step so as to make the spectral radius <= 1.000001.')
        disp('Choosing default method type: "implicit"');
        disp(' ')
        M=(inv((I + ((1i*tau)/(2*hbar)).*H)))*(I - ((1i*tau)/(2*hbar)).*H);
        method = 'implicit';
    else
        disp(' ')
        disp(sprintf('WARNING(seq_crank): Spectral radius: %.7f is greater than unity, FTCS solution will be unstable at large t. Running anyways.', r));
        disp(' ')
    end
    
        
end

%Initialize x spatial array
x = linspace(-L/2,L/2, Nspace);
xcentered = x - x0;


%Implement periodicity in the initial condition
%I still get some wobblies when I choose x0 near the edge but this is
%better than not having any sort of periodicity
j = find(xcentered < -L/2);
xcentered(j) = xcentered(j) + L;
j = find(xcentered > L/2);
xcentered(j) = xcentered(j) - L;


%Initialize the Wave Function 
psi = zeros(length(x), Ntime);
psi(:,1) = (1/sqrt(sig0*(pi^(1/2)))).*exp((1i*k0).*x).*exp((-(xcentered).^2)./(2*(sig0^2)));

%Clear whatever figure is open 
clf
hold off

figure(1);

%If the user chooses to plot psi
if strcmp(plotType, 'psi')  
    
    %Grab the axisheight as the maximum value of psi if we have a large wave function
    axisheight = max(real(psi(:,1)));
    if axisheight < 0.5
        axisheight = .5;
    end
    
    %Integrate from 1 to the total time
    for t = 1:Ntime
        
        
        %Calculate forward time psi using matrix multiplication
        psi(:,t+1) = M*psi(:,t);
        
        %Plot psi at the current time
        plot(x, real(psi(:,t)), 'b');
        
        %Plot potential spikes as vertical lines
        for j = 1:length(potential_spikes)
            line([x(potential_spikes(j)) x(potential_spikes(j))],[-axisheight axisheight],'Color',[1 0 0], 'LineWidth', 0.5)
        end
        %Set the limits of the axes
        axis([-L/2 L/2 -axisheight axisheight]);
        
        
        
        %Label the axes
        xlabel('x');
        ylabel('psi(x,t)');
        title(sprintf('%s time evolution of the Schrodinger Wave Equation. Psi(x,t) vs. x',upper(method)));  
        %Draw the current plot
        drawnow;
    end
        

     
%If the user chooses to plot probability
elseif strcmp(plotType, 'prob')
    
    %For some reason whenever I called with Ntime == Nspace, my probability
    %conservation plot was off. I couldn't quite figure out why, so I
    %figured I might as well leave a warning if nothing else
    if Ntime == Nspace
        disp(' ')
        disp('WARNING(seq_crank): probability conservation plot does not work for Ntime == Nspace.')
        disp(' ')
    end
    
    %Initialize the probabliity distribution
    prob = zeros(length(x), Ntime);
    prob(:,1) = psi(:,1).*conj(psi(:,1));
    
    %Grab the axisheight as the maximum value of probability if we have a large wave function
    axisheight = max(prob(:,1));
    if axisheight < 0.1
        axisheight = 0.1;
    end
    
    %Initialize the total probability
    totalProb = zeros(length(x), Ntime);
    totalProb(:,1) = trapz(x, prob(:,1));
    
    %Integrate from 1 to the total time
    for t = 1:Ntime
        
        %Calculate forward time psi using matrix multiplication
        psi(:,t+1) = M*psi(:,t);
        
        %Calculate probablility density
        prob(:,t+1) = psi(:,t+1).*conj(psi(:,t+1));
        
        %Calculate the total probability using the inner product 
        totalProb(:,t+1) = trapz(x, prob(:,t+1));
        
        %Plot the probability density at the current time
        plot(x, prob(:,t), 'r');
        
        %Plot potential spikes as vertical lines
        for j = 1:length(potential_spikes)
            line([x(potential_spikes(j)) x(potential_spikes(j))],[-axisheight axisheight],'Color',[0 0 1], 'LineWidth', 0.5)
        end
        
        %Set the axes
        axis([-L/2 L/2 0 axisheight]);  
        
        %Label the axes
        xlabel('x');
        ylabel('P(x,t)');
        title(sprintf('%s time evolution of the Schrodinger Wave Equation. Prob(x,t) vs. x', upper(method)));
        %Draw the current plot
        drawnow;
    end
    %Open up a new figure
    figure(2);
    
    %Make the time array again for plotting total probability
    t = 1:Ntime;

    %Plot the base 10 logarithm of the absolute value of 1 - total
    %probability as a function of time, to show the difference between the
    %expected value of probability vs. the observed value.
    plot(t, log10(abs(1-totalProb(:,t))), 'k');
    title('log10(|1-totalprobability(t)|) vs. t');
    %Label the axes
    xlabel('t');
    ylabel('log10 |1 - |P(x,t)|^2|');
    
%If the user chooses to party
elseif strcmp(plotType, 'party')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'Color','k')
    [y, Fs] = audioread('boogiewonderland.mp3');
    sound(y, Fs, 16);
    figure('units','normalized','outerposition',[0 0 1 1])
    disp(' ');
    disp('To stop party use ctrl+c and type: clear sound');
    disp('THANKS FOR A GREAT SEMESTER DR. SIGUT');
    disp(' ');
    
    %Grab the axisheight as the maximum value of psi if we have a large wave function
    axisheight = max(real(psi(:,1)));
    if axisheight < 0.5
        axisheight = .5;
    end
    
    %Integrate from 1 to the total time    
    for t = 1:Ntime
        %Calculate forward time psi using matrix multiplication
        psi(:,t+1) = M*psi(:,t);
        
        %Plot psi at the current time
        plot(x, real(psi(:,t)), 'LineWidth', 2 + 5*abs(cos(0.1*t)), 'color', [ abs(cos(0.1*t)) (1 - abs(sin(0.075*t))) abs(sin(0.2*t))]);
        
        %Plot potential spikes as vertical lines
        for j = 1:length(potential_spikes)
            line([x(potential_spikes(j)) x(potential_spikes(j))],[-axisheight axisheight],'Color',[round(abs(cos(0.1*t))) round(abs(sin(0.3*t))) round(abs(cos(0.2*t)))], 'LineWidth', 2 + 10*abs(cos(0.1*t)))
        end
        
        %Party Mode
        set(gca,'Color','k')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        %plot(x, V(x))
        txt = 'PARTY ON!!!';
        txt = text((sin(0.05*t)*L)/4,0.3*cos(0.05*t),txt, 'color', [round(abs(cos(0.2*t))) round(abs(sin(0.1*t))) round(abs(cos(0.3*t)))]);
        txt.FontSize = 40 + 10*sin(0.05*t);
        txt.Rotation = t;
        txt.HorizontalAlignment = 'center';
        title('Party vs. Party')
        %Set the limits of the axes
        axis([-L/2 L/2 -axisheight axisheight]);
        title('Party vs. Party')
        
        %Draw the current plot
        drawnow;
    end
    %End Party
    clear sound;
 end




