function sol = multigaussfit(x, y, varargin)
% multigaussfit is a gaussian fit function. User set constrained or not
% sol = multigaussfit(x, y, varargin)

%% Parse inputs
p = inputParser;

% Common parameters
addOptional(p, 'ngauss', []); % Number of gaussians
addOptional(p, 'initial', []); % Initial parameters with each column defining a gaussaian function
addOptional(p, 'unconstrained', false); % Unconstrained fit
addOptional(p, 'logy', false); % Fit with log scale yaxis

% Auto parameters
addOptional(p, 'auto', false); % Automode that uses local peaks to determine the number of gaussians and iniatial points.
addOptional(p, 'signamethod', 'scale'); % Method to determine sigma: 'scale' uses a scale of mean, 'sqrt' uses sqrt(mean)
addOptional(p, 'meanscale', 0.5); % Scale of mean for sigma

% User input parameters
addOptional(p, 'userinput', true); % Use user inputs for initial points. Only applies if no initial points were given.

% Plot
addOptional(p, 'plotresults', true); % Plot results
addOptional(p, 'plotcomponents', true); % Plot components
addOptional(p, 'displayfunctions', true); % Display function texts

% Smoothing
addOptional(p, 'smoothingwin', 0); % Window size for smoothing (0 to disable)

% fminsearch parameters (unconstrained fitting)
addOptional(p, 'ToxX', 1e-4);
addOptional(p, 'MaxFunEvals', 10^12);
addOptional(p, 'TolFun', 1e-4);
addOptional(p, 'MaxIter', 100000);

% fmincon parameters (constrained)
addOptional(p, 'A', []); % Inequality A * x <= b
addOptional(p, 'b', []); % Inequality A * x <= b
addOptional(p, 'Aeq', []); % Inequality Aeq * x = beq
addOptional(p, 'beq', []); % Inequality Aeq * x = beq
addOptional(p, 'lb', []); % Lower bound (all 0s if left empty)
addOptional(p, 'ub', []); % Upper bound (limited by xmax and ymax if left empty)
addOptional(p, 'nonlcon', []); % Nonlinear constraint
addOptional(p, 'StepTolerance', 1e-9); % Step tolerance for constrained fitting
% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;


%% Fix inputs
% Fix inputs if needed
x = reshape(x, 1, []);
y = reshape(y, 1, []);

if ~p.auto && ~isempty(p.initial)
    l = length(p.initial(:));
    
    if mod(l, 3) ~= 0
        error('Initial condition has a wrong number of components')
    elseif l/3 ~= p.ngauss
        error('Number of gaussians and initial conditions do not match')
    end
end


%% Preprocessing
% Smoothing
if p.smoothingwin > 0
    y = smooth(y, p.smoothingwin);
end

% Logy
if p.logy
    y = log(y);
end

%% Initialize functions
% Single gauss
singlegauss = @(x, p) p(1)*exp(-((x-p(2))/p(3)).^2);

% Multi gauss function
function yf = multigauss(x, p)
    % Parameters
    ngauss = size(p, 2);
    yf = zeros(ngauss, length(x));
    
    % Calculate each gaussian
    for ii = 1 : ngauss
        yf(ii,:) = singlegauss(x, p(:,ii));
    end
    
    % Combine
    yf = sum(yf, 1);
    
end


% Error function
if p.logy
    gausserror = @(p) norm(log(multigauss(x, p)) - y);
else
    gausserror = @(p) norm(multigauss(x, p) - y);
end

%% Initial fit parameters
if p.auto
    % Fit peaks (Min prominence is 5% of max y)
    maxy = max(y);
    [pks, locs, ~, ~] = findpeaks(y, x, 'MinPeakProminence', maxy/20);
    
    % Find the number of gauss
    p.ngauss = length(pks);
    
    % Intial condition
    params = zeros(3, p.ngauss);
    params(1,:) = pks;
    params(2,:) = locs;
    switch p.signamethod
        case 'scale'
            params(3,:) = locs * p.meanscale;
        case 'sqrt'
            params(3,:) = sqrt(locs);
    end
    
elseif ~isempty(p.initial)
    % Take input
    params = reshape(p.initial, 3, []);
    
elseif p.userinput && ~isempty(p.ngauss)
    % Plot
    figure
    plot(x, y);
    
    params = zeros(3, p.ngauss);
    
    % Get points
    for i = 1 : p.ngauss
        % Title
        title(sprintf('Draw boxes to choose Peak #%i out of %i', i, p.ngauss));
        
        % Draw box
        hbox = imrect();
        coord = wait(hbox);
        
        % Box parameters
        s = coord(3)/2;
        gx = coord(1) + s;
        gy = coord(4);
        
        % Label
        text(gx, gy, num2str(i));
        
        % Log
        params(1,i) = gy;
        params(2,i) = gx;
        params(3,i) = s;
        
        % Delete
        delete(hbox)
    end
    
    close gcf
else
    error('No method for determining initial conditions') 
end

%% Unconstrained fitting
% Unconstrained
if p.unconstrained
    % First, set some options for fminsearch().
    options = optimset('TolX', p.ToxX, 'MaxFunEvals', p.MaxFunEvals, 'TolFun', ...
        p.TolFun, 'MaxIter', p.MaxIter);  % Determines how close the model must fit the data

    % Fit
    sol = fminsearch(gausserror, params, options);
end

%% Constrained fitting
if ~p.unconstrained
    % Lower bounds
    if isempty(p.lb)
        p.lb = zeros(size(params));
    end
    
    % Upper bounds
    if isempty(p.ub)
        p.ub = zeros(size(params));
        p.ub(1,:) = inf;
        p.ub(2:3,:) = max(x);
    end
    
    options.StepTolerance = p.StepTolerance;
    
    % Fit
    sol = fmincon(gausserror, params, p.A, p.b, p.Aeq, p.beq, p.lb, p.ub, p.nonlcon, options);
end

%% Plotting
if p.plotresults
    % Plot
    plot(x, [y', multigauss(x, sol)']);
    
    title('Fit result')
    legend({'Data', 'Fit'})
    
    if p.plotcomponents
        hold on
        for i = 1 : p.ngauss
            plot(x, singlegauss(x, sol(:,i)));
        end
        hold off
    end
    
    if p.displayfunctions
        yaxis = ylim();
        xaxis = xlim();
        
        for i = 1 : p.ngauss
            % p(1)*exp(-((x-p(2))/p(3)).^2);
            fdisplay = sprintf('y%i = %i*exp(-((x-%i)/%i)^2', i, ...
                round(sol(1,i)), round(sol(2,i)), round(sol(3,i)));
            
            text(xaxis(2)*0.5, yaxis(2)*0.5 - (i - 1) * yaxis(2)/20, fdisplay);
        end
    end
    
end

end