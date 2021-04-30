# multigaussfit
 Multi gaussian fit
 
 Automode: flag auto to be true
 e.g., sol = multigaussfit(x, y, 'auto', true);
 
 Numerical initial condition: specify the exact number of gaussians and initial parameters (p) in the form of (in rows as a vector): Amplitude, mean, spread (scaled sigma).
 e.g., sol = multigaussfit(x, y, 'auto', false, 'ngauss', 2, 'initial', [1 2; 2 10; 2 2]); # Initial amplitudes = 1 and 2; means = 2 and 10; spreads = 2 and 2
 
 Hand-drawn initial condition: specify the exact number of gaussians but leave initial parameters empty. Draw rectangles to specify initial conditions when prompt to.
 e.g., sol = multigaussfit(x, y, 'auto', false, 'ngauss', 2, 'initial', []); # Hand draw initial conditoins
 
 
