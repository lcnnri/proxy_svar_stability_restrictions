**********************************************************************
Replication package for

Angelini, G., Fanelli, L., and Neri, L. (2025)
"Invalid Proxies and Volatility Changes"
Journal of Business & Economic Statistics
**********************************************************************

----------------------------------------------------------------------
I. OVERVIEW
----------------------------------------------------------------------

This package contains MATLAB code and supporting material to replicate
the empirical results reported in the paper above.

The replication reproduces the main figures and tables of the paper,
including the proxy-SVAR estimation with stability restrictions.

Link to the repository:
https://github.com/lcnnri/proxy_svar_stability_restrictions

----------------------------------------------------------------------
II. HOW TO RUN
----------------------------------------------------------------------

1. Open MATLAB (tested with R2025b Update 1 or later).
2. Set the working directory to the folder "matlab".
3. Run the main replication script:

       main_fiscal.m

The script automatically:
   - loads the data and configuration settings
   - estimates the proxy-SVARs (full-sample and regime-dependent)
   - produces impulse responses, tables, and plots
   - saves the results to the folder "results"

Results are written to "matlab/results/", including:
   - figures/ : impulse responses and multipliers
   - tables/  : estimated parameters and confidence intervals
   - logs/    : console output and run diagnostics

----------------------------------------------------------------------
III. FILE ORGANIZATION
----------------------------------------------------------------------

.
├── matlab/
│   ├── main_fiscal.m       Main replication script (run this)
│   ├── scripts/            High-level routines
│   ├── functions/          Auxiliary functions
│   ├── data/               Placeholder for data files
│   └── results/            Logs, tables, figures, workspaces
├── Supplement-AFN.pdf      Supplementary material to the paper
├── LICENSE                 License file
└── README.txt              This file

----------------------------------------------------------------------
IV. SOFTWARE INFORMATION
----------------------------------------------------------------------

Tested environment:
   MATLAB Version: R2025b Update 1 (25.2.0.3042426)
   Operating System: Windows 11 Home, Build 26100
   Java Version: 1.8.0_202-b08, Oracle HotSpot 64-Bit

Required MATLAB toolboxes:
   - Optimization Toolbox
   - Statistics and Machine Learning Toolbox
   - Global Optimization Toolbox
   - Econometrics Toolbox

**********************************************************************
End of README.txt
**********************************************************************
