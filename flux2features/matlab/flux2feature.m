function features = flux2feature(model, fluxData, flux2Compare)
% flux2feature Extract interpretable features from metabolic flux data
%
% This function computes a set of mathematical and biologically meaningful 
% features from a single flux vector or a collection of flux samples 
% (e.g., obtained via CHRR/ACHR sampling) for a given COBRA model.
%
% USAGE:
%    features = flux2feature(model, fluxData)
%    features = flux2feature(model, fluxData, flux2Compare)
%
% INPUTS:
%    model:        COBRA model structure with fields:
%                     * .S    - Stoichiometric matrix (m x n)
%                     * .b    - Right-hand side vector (m x 1)
%                     * .lb   - Lower bounds (n x 1)
%                     * .ub   - Upper bounds (n x 1)
%                     * .mets - Metabolite identifiers
%
%    fluxData:     Flux data (n x 1) or (n x k)
%                     * Single flux vector, or
%                     * Matrix of k flux samples
%
%                     The function automatically detects the input type.
%
%    flux2Compare: (optional) Reference flux vector (n x 1)
%                     Used to compute comparative features 
%                     (e.g., distances, similarity metrics)
%
% OUTPUT:
%    features:     Table containing extracted features
%                     * Rows correspond to samples (k)
%                     * Columns correspond to features
%
% NOTES:
%    - Designed for flux samples generated
%    - Includes intrinsic, comparative, and biochemical features
%    - If flux2Compare is not provided, comparative features are omitted
%    - Intended as a core component for metabolic state analysis pipelines
%
% .. Author: German Preciat (2026)