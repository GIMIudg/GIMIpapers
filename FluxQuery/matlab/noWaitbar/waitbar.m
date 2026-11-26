function h = waitbar(varargin)
% Silent replacement of MATLAB's waitbar, used only by samplesGenFluxQuery.m.
% The CHRR sampler of the COBRA Toolbox (genSamples) opens a waitbar and updates it
% after every sample; on macOS this graphics call can hang MATLAB, and in a terminal
% run (matlab -batch) there is no display for it. It does not use random numbers,
% so the samples are identical with or without it. Progress is printed by
% samplesGenFluxQuery.m after each set.
h = [];
end
