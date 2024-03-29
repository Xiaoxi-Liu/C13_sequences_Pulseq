function [grad,npts] = gradient_lobe(ga_des,dt,gmax,smax,g_off,rnd2even,...
    rev,verb)
% GRADIENT_LOBE Calculate gradient lobe with a initial gradient value 
% for desired area. Results in triangluar or trapezoidal gradient shape.
%
%   grad = gradient_lobe(ga_des,dt,gmax,smax,g_off,rnd2even,rev,verb)
%                                                [unit]    (default)
%  ga_des  Desired area of gradient lobe         [s*T/m]
%      dt  Sampling dwell time                   [s]
%    gmax  Maximum gradient strength             [T/m]
%    smax  Maximum slew rate                     [T/m/s]
%   g_off  Offset gradient value                 [T/m]     (0)
%rnd2even  Round up to even #waveform pts                  (false)
%     rev  Reverse grad (start with 0 and end with g_off)  (false)
%    verb  Verbose mode                                    (false)
%
%    grad  Gradient waveform [1,#pts]            [T/m]
%    npts  #grad waveform points
%
% 8/2019  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('dt','var'),       dt = []; end
if isempty(dt),              dt = 4d-6; end
if ~exist('gmax','var'),     gmax = []; end
if isempty(gmax),            gmax = 33d-3; end
if ~exist('smax','var'),     smax = []; end
if isempty(smax),            smax = 120; end
if ~exist('g_off','var'),    g_off = []; end
if isempty(g_off),           g_off = 0; end
if ~exist('rnd2even','var'), rnd2even = []; end
if isempty(rnd2even),        rnd2even = false; end
if ~exist('rev','var'),      rev = []; end
if isempty(rev),             rev = false; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = false; end

if ~isreal(ga_des)
    warning('ga_des is complex number; taking abs');
    ga_des = abs(ga_des);
end
if gmax<0, error('gmax(=%g)<0',gmax); end
if smax<0, error('smax(=%g)<0',smax); end
if abs(g_off)>gmax
    warning('abs(g_off)(=%g)>gmax(=%g)',abs(g_off),gmax);
end
if abs(ga_des/gmax/dt)>1d10
    error('abs(ga_des/gmax/dt(=%g))>1',ga_des/mgax/dt);
end
if abs(ga_des/gmax/dt)>1d6
    warning('abs(ga_des(=%g))>1d-4',ga_des/gmax/dt);
end


%% calculate positive lobes; invert when necessary
sgn    = sign(ga_des);
if abs(sgn)<1d-10, sgn = 1; end
ga_des = abs(ga_des);
g_off  = sgn*g_off;

g_nom = sqrt(smax*abs(ga_des)+0.5*g_off^2);    % nominal gradient strength
if g_nom<g_off        % if offset larger, invert again for lobe pointing up
     sgn = -sgn;
     ga_des = -ga_des;
     g_off = -g_off;
     g_nom = sqrt(-smax*abs(ga_des)+0.5*g_off^2);
end
if verb, fprintf('g_nom=%g; g_off=%g\n',g_nom,g_off); end


%% differentiate between trapezoids and triangles
if g_nom>gmax % trapezoid: add gradient plateau in middle
    if verb, fprintf('Constructing trapezoid\n'); end
    
    % points for ramp up and down
    n_up   = ceil((gmax-g_off)/smax/dt)+1;
    n_down = ceil(gmax/smax/dt)+1;
    
    % area for ramp up and down
    ga_up   = 0.5*(gmax^2-g_off^2)/smax;
    ga_down = 0.5*gmax^2/smax;
    
    n_plt = ceil((ga_des-(ga_up+ga_down))/gmax/dt);
    % gmax_act = gmax;
else          % triangular gradient
    if verb, fprintf('Constructing triangle\n'); end
    n_up   = ceil((g_nom-g_off)/smax/dt);
    if n_up==1
        warning('n_up==1; setting to 0');
        n_up = 0;
    else
        n_up = n_up+1;
    end
    n_down = ceil(g_nom/smax/dt)+1;
    % g_nom_act = g_nom;
    
    if n_up<0, warning('n_up<0'); end
    
    n_plt = 0;
end
nn = n_up+n_down+n_plt;
% if rnd2even && (abs(nn/2-floor(nn/2))>eps)
%         n_plt = n_plt+1;
% end


%% calculate exact g_nomax_act to prevent rounding errors
g_nom_act = (2*ga_des/dt-n_up*g_off)/(n_up+n_down+2*n_plt);


%% construct actual gradient waveform
grad = sgn*[linspace(g_off,g_nom_act,n_up),...
    g_nom_act*ones(1,n_plt),...
    linspace(g_nom_act,0,n_down)];

if rnd2even && (abs(nn/2-floor(nn/2))>eps)
    grad = [grad 0];
end


%% reverse gradient waveform
if rev, grad = grad(1,end:-1:1); end


%% check waveform
ga_act = sum(grad)*dt;
% ga_err = (sgn*ga_des-ga_act)/ga_des;
ga_err = sgn*ga_des-ga_act;
smax_act = max(abs(diff(grad)))/dt;
gmax_act = max(grad);
if smax_act>smax, warning('smax_act(=%g)>smax(=%g)',smax_act,smax); end
if gmax_act>gmax, warning('gmax_act(=%g)>gmax(=%g)',gmax_act,gmax); end
if abs(ga_err)>1d-6
    warning('abs(ga_err(=%g))>1d-6; ga_des=%g; ga_act=%g',...
        ga_err,ga_des,ga_act); 
end
% npts = n_up+n_plt+n_down;
npts = size(grad,2);

%% print info
if verb
    fprintf('gmax=%g; gmax_act=%g [T/m]\n',gmax,gmax_act);
    fprintf('smax=%g; smax_act=%g [T/m]\n',smax,smax_act);    
    fprintf('n_up=%g; n_plt=%g; n_down=%g\n',n_up,n_plt,n_down);
    fprintf('ga_des=%g; ga_act=%g; err=%g[%%]\n',...
        sgn*ga_des,ga_act,ga_err*100);
    fprintf('duration=%g [ms]\n',npts*dt*1d3);
end
if rnd2even && (abs(size(grad,2)/2-floor(size(grad,2)/2))>eps)
    warning('waveform not even');
end


end   % end main function gradient_lobe.m
