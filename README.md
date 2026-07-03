# Proxy-SVAR with Stability Restrictions

This repository contains MATLAB code to replicate the empirical results in:

> **G. Angelini, L. Fanelli, and L. Neri (2026),** *Invalid Proxies and Volatility Changes*,
> *Journal of Business & Economic Statistics*, 1–27.
> DOI: [10.1080/07350015.2025.2602857](https://doi.org/10.1080/07350015.2025.2602857)

Replication guide:
https://lcnnri.github.io/proxy_svar_stability_restrictions/

## Overview

This package reproduces the empirical results reported in the main body of the paper, including the proxy-SVAR estimation with stability restrictions. It generates Figure 1 and Table 3.

Repository:
https://github.com/lcnnri/proxy_svar_stability_restrictions

## How to Run

1. Open MATLAB. The code was tested with MATLAB R2025b.
2. Set the working directory to the folder `matlab/`.
3. Run the main replication script:

```matlab
main_fiscal;
```

The script automatically loads the data and configuration settings, estimates the proxy-SVARs, including the full-sample and regime-dependent specifications, produces impulse responses, tables, and plots, and saves all results under `matlab/results/`.

Output folders under `matlab/results/`:

```text
figures/   Figure 1
tables/    Table 3
logs/      Console output and run diagnostics
```

## Repository Layout

```text
.
├── matlab/
│   ├── main_fiscal.m      # main replication script
│   ├── scripts/           # high-level routines
│   ├── functions/         # auxiliary functions
│   ├── data/              # data placeholder
│   └── results/           # logs, tables, workspaces, figures
├── docs/
│   └── index.html         # HTML replication guide
├── Supplement-AFN.pdf     # supplementary material to the paper
├── LICENSE                # license file
├── README.txt             # txt version for journal submission
└── README.md              # this file
```

## Software Information

Tested environment:

* MATLAB Version: R2025b Update 1, 25.2.0.3042426
* Operating System: Windows 11 Home, Build 26100
* Java Version: 1.8.0_202-b08, Oracle HotSpot 64-Bit

Required MATLAB toolboxes:

* Optimization Toolbox
* Statistics and Machine Learning Toolbox
* Global Optimization Toolbox
* Econometrics Toolbox

All scripts are self-contained and require only the standard MATLAB toolboxes listed above. No third-party MATLAB packages are used.

## Citation

If you use this replication package or parts of the code, please cite:

```bibtex
@article{AngeliniFanelliNeri2026,
  author    = {Giovanni Angelini and Luca Fanelli and Luca Neri},
  title     = {Invalid Proxies and Volatility Changes},
  journal   = {Journal of Business \& Economic Statistics},
  year      = {2026},
  pages     = {1--27},
  publisher = {Taylor \& Francis},
  doi       = {10.1080/07350015.2025.2602857},
  url       = {https://doi.org/10.1080/07350015.2025.2602857}
}
```

An earlier working-paper version is available at:

https://arxiv.org/abs/2403.08753

## License

This work is licensed under the [MIT License](LICENSE).
