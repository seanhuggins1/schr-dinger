%% seq_crank
% Integrates the time dependent Schrodinger Equation using a
% Crank-Nicholson or FTCS scheme, and produces an animated plot of the time
% evolution.
%
%% Syntax
%
%   seq_crank(Nspace, Ntime, tau)
%   seq_crank(Nspace, Ntime, tau, varargin)
%
%% Description
% |seq_crank(Nspace, Ntime, tau)|, where Nspace is the number of spatial
% grid points, and Ntime is the number of tau steps, integrates the time
% dependent Schrodinger Equation. Using either the Crank-Nicholson scheme
% or a forward-time-centered-space scheme, seq_crank evolves the wave
% function psi(x,t) forward for a total time of Ntime*tau.
% seq_crank will produce an animation of phi(x,t) vs. x.
%%
% |seq_crank(Nspace, Ntime, tau, varargin)| will produce an animation of
% the time evolution of the Schrodinger Equation with options given by
% varargin. Using the 'Method' varargin, the user may choose between using an implicit Crank-Nicholson
% scheme which is unconditionally stable for all tau, or they may choose to
% use the explicit FTCS scheme given a small enough time step. The length of
% the spatial domain can be changed using the 'Length' varargin. A numeric
% value can be chosen to alter the domain corresponding to 
% |-L/2 < x < L/2|. The user may also choose to produce different plot
% animations via. the 'PlotType' varargin. Animations can be made of the wave
% function psi(x,t), the probability density prob(x,t) and the conservation
% of total probability. The user can also plot 'party' for a great time.
% The user can set spatial indicies where the potential V(x) shall be set
% to 1, using the 'Potential' varargin. Lastly, the user may alter the
% initial condition of the wave function by changing the parameters
% |[sigma0 x0 k0]| using the 'WParam' varargin. See *Examples of Use* for
% more information on how to call |seq_crank| using various varargin.
%
%% Examples of Use
% *Running with default parameters*
seq_crank(300, 250, 0.1); 
%%
% Here we have integrated the wave function using the default parameters,
% as well as Nspace = 300, and Ntime = 250.
% Default parameters for varargin are as follows:
%
% |'Method' = 'Implicit'|
%
% |'Length' = 200|
%
% |'PlotType' = 'psi'|
%
% |'Potential' = []|
%
% |'WParam' = [10 0 0.5]| corresponding to |[sigma0 x0 k0]|
%
%
%%
% *Probability plots*
seq_crank(300, 250, 0.1, 'PlotType', 'prob');
%%
% The first plot shows probability density vs. x at the final time
% step, while the second plot shows how the total probability was conserved
% throughout the integration by plotting the base 10 logarithim of the
% absolute value of 1 - total probability, highlighting the difference
% between the expected probability |1| and the observed probability
% calculated by the inner product of psi and its complex conjugate.
% I highly reccomend the user choose 'PlotType', 'party', 
% with Ntime ~ 3000, as they will not regret their decision.
%
%
%%
% *Plotting the wave function psi(x,t) with potential spikes Vj*
seq_crank(300, 250, 0.1, 'Potential', [25 50]);
%%
% Here seq_crank sets the potential V(x) = 1 at every j index given by
% the user. In this example the chosen j indicies are 25 and 50. Take
% careful note that these are spatial indicies, and they do not
% corresponds to positions along the domain |-L/2 < x < L/2|, rather they
% correspond to an integer index j, where |1 < j < Nspace|.
% Potential spikes are shown as vertical lines.
%
%
%%
% *Plotting the wave function psi(x,t) using the explicit FTCS method*
seq_crank(300, 250, 0.0000001, 'Method', 'explicit');
%%
% In this example we have chosen to plot using the explicit FTCS scheme. As
% explained below in *Varargin 'Method'*, the time integration of the
% schrodinger equation is unconditionally unstable when we use the explicit
% FTCS scheme, which is why we had to use such a small time step here.
% seq_crank will only plot the explicit solution if the time step is small
% enough as to make the spectral radius <= 1.000001. Even then, the program
% has warned us that the solution will be unstable at large Ntime. It is
% reccomended that the user chooses *'Method', 'implicit'* as the implicit
% scheme is unconditionally stable for all tau.
%
%
%%
% *Plotting the wave function psi(x,t) using different parameters*
seq_crank(300, 250, 0.1, 'Wparam', [15 25 0.3]);
%%
% The user may choose to alter the initial condition of the wave function
% by changing the parameters using the *'Wparam'* varargin. The user must
% pass in a 3 element array corresponding to parameters |[sigma0 x0 k0]|. See
% *Varargin 'WParam'* for more details on the allowable values for 
% |[sigma0 x0 k0]|.
%
%
%%
% *Plotting the wave function psi(x,t) and changing the length*
seq_crank(300, 250, 0.1, 'Length', 400);
%%
% Here we have changed the length of our domain corresponding to 
% |-L/2 < x < L/2.|
%%
% *Using a little bit of everything*
seq_crank(300, 250, 0.1, 'Method', 'Implicit', 'Length', 300, 'PlotType', 'psi', 'Potential', [50 70:75 240], 'Wparam', [15 -10 0.4]);
%%
% Above is an example where we have used every possible varargin. We have
% chosen the implicit Crank-Nicholson scheme, a domain length of 300, we
% are plotting psi, we have potential spikes at 50, 70-75, and 0.4, and we
% have changed the default parameters to |[sigma0 x0 k0] = [15 -10 0.4]|.
%
%
%% Input Arguments
% *Nspace* - Number of spatial grid points
% The user may choose the number of spatial grid points to divide the
% domain into. The spatial step size h is calculated using Nspace and
% Length through |h = L/(Nspace-1)|.
% Nspace must be a finite postive integer greater than zero. 
% Choose Nspace to be ~= L for a nice smooth solution.
%
%%
% *Ntime* - Number of time steps
% The user can choose the number of time steps to forwardly integrate the
% wave equation. The wave function will advance from its initial state,
% through a total time of Ntime*tau, in steps of tau.
% Ntime must be a finite postive integer greater than zero. 
% The user may choose to run the simulation for however long they like.
%
%%
% *tau* - Length of each time step
% tau is the time step. A smaller time step will advance the wave function
% slower, and a larger time step will advance the wave function faster.
% For the Crank-Nicholson "implicit" scheme, the user may choose any tau they wish, as
% the Crank-Nicholson scheme is unconditionally stable for all tau.
% As for the FTCS "explicit" scheme, unfortunately the time evolution of the wave
% function is unconditionally UNSTABLE for all tau. The user may choose to
% run the FTCS scheme with a spectral radius <= 1.000001, but this will still be
% unstable and instabilities will become noticable after a long enough
% Ntime. The user can reduce the spectral radius by choosing a smaller tau. 
% A warning will be given if the user chooses such a spectral radius with
% the FTCS scheme. It is reccomended that the user run seq_crank with the
% "implicit" Crank-Nicholson sceme, as they may choose any tau they wish
% and achieve stability.
%
%% Varargin
% *'Method'* - Integration method to use in the time evolution of the wave function
% options for method are *'implicit'* for Crank-Nicholson, and *'explicit'* for forward-time-centered-space integration.
% As stated above in the description for tau, it is reccomended that the
% user choose *'implicit'* as for any choice of tau, the *'explicit'* FTCS
% scheme will be unstable. The user may choose to run *'explicit'* if they
% wish, however a spectral radius <= 1.000001 is required. The user can
% reduce the spectral radius by choosing a smaller tau.
% Default method = |'implicit'|
%
%%
% *'Length'* - Length of the spatial domain
%%
% The user may choose a length L for the spatial domain which is given by 
% |-L/2 < x < L/2|, and divided into Nspace steps.
% The length must be a finite positive real number greater than zero.
% A length >= 200 is reccomended to be able to observe the full behavior of
% the wave function.
% Default length = |200|
%
%%
% *'PlotType'* - Type of plot to show
%%
% Three plotting options are available to the user.
% Choosing *'psi'* brings up an animation of the time evolution of the wave
% function psi(x,t), which runs for Ntime. 
% Choosing *'prob'* brings up an animation of the time evolution of the
% probability density psi(x,t)*conj(psi(x,t)). The total probability,
% being the inner product of psi(x,t) and conj(psi(x,t)), is also shown and
% plotted as the base 10 logarithm of the absolute value of 1 - total probability
% as a function of time to clearly show how well probability is conserved
% throughout the integration.
% If it has not already been made clear, the user may choose to plot
% *'party'* to have nothing but a great time. It is reccomended that
% when calling with the plot type *'party'*, the user chooses a long Ntime,
% say ~ 3000. If the user is not having as much fun as they originally
% hoped for, they may hit |ctrl+c|, and type |clear sound| in the console to
% leave the party. The user should also make sure that they have the file
% 'boogiewonderland.mp3' in the same directory as 'seq_crank.m' before they
% plot *'party'*.
% Default plot type = |'psi'|
%
%%
% *'Potential'* - Spatial indicies j of the potential spikes
%%
% The user may pass in spatial indicies for where they wish to set the potential V(x) = 1
% The indicies should be positive integers between |1 <= j <= Nspace|. Any
% choice of indicies is possible, for example when Nspace is 500, passing
% in 'Potential', |[40:100 250 477:487]| will set V(x) = 1 from 40-100, at
% 250, and from 477-487. Plots will show vertical lines where there are
% these potential spikes.
% Default array for potential spikes = |[]|
%
%%
% *'Wparam'*
%%
% The user may pass in an array of parameters |[sigma0 x0 k0]| to alter the
% initial condition of the wave function. sigma0 x0 and k0 can be finite
% real numbers while x0 is constrained to the domain |-L/2 < x < L/2|.
% Default WParam = |[10 0 0.5]|;
%
%% Output Arguments
% seq_crank produces no output other than animations of various plots (see *Varargin 'PlotType'*)
%