# app-template-python

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Removes the background field from a QSM total field map with V-SHARP (QSM.jl in Julia), producing a tissue (local) field map and the eroded V-SHARP mask.**

Installs Julia 1.9.3 and the QSM.jl, NIfTI.jl and MriResearchTools packages at run time, masks the field map, converts it to Hz using the field map units (from input_units or the fieldmap JSON) and the field strength, and runs V-SHARP background field removal. Use it between phase unwrapping/field estimation and dipole inversion in a QSM pipeline. Needs internet access on the compute node; the README is an unfilled template.

| | |
|---|---|
| 📥 **Inputs** | `neuro/fmap` · `neuro/mask` |
| 📤 **Outputs** | `neuro/fmap` · `neuro/mask` |
| ⚙️ **Resources** | 16 GB memory · 1 h walltime |
| 🔗 **Run it** | [doi.org/10.25663/brainlife.app.835](https://doi.org/10.25663/brainlife.app.835) |

#### 📥 Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `fieldmap` | `neuro/fmap` (single) | yes | Total B0 field map with magnitude and JSON sidecars. Field strength must be in the input metadata (MagneticFieldStrength) or B0; units must be in the fieldmap JSON (Units) or input_units. |
| `mask` | `neuro/mask` | yes | Brain mask restricting the field map. |

#### 📤 Outputs

| Output | Datatype | Description |
|---|---|---|
| `fmap` | `neuro/fmap` | Magnetic field variations caused by tissues (background field removed). |
| `mask` | `neuro/mask` | V-SHARP mask |

<!-- brainlife-app:end -->

[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.444-blue.svg)](https://doi.org/10.25663/bl.app.444)

## Documentation

This is a template for a python-based brainlife.io/app. Please add the description of the app here.

You can include: 
1. What the App does, and how it does it at the basic level. 
2. Briefly description of input / output files.

## Authors
- [Author1](email), Affiliation

## Contributors
- [Author2](email), Affiliation
- [Author3](email), Affiliation

## Citations
We kindly ask that you cite the following articles when publishing papers and code using this app. 

1. **brainlife. io: A decentralized and open source cloud platform to support neuroscience research**. Hayashi, S., Caron, B. A., et al. & Pestilli, F. (2023). ArXiv. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10274934/

2. **Add any other citations based on the app** 

## Funding Acknowledgement
brainlife.io is publicly funded. For the sustainability of the project, we kindly ask that you acknowledge the following funding sources:

[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-ACI-1916518](https://img.shields.io/badge/NSF_ACI-1916518-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1916518)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NIH-NIBIB-R01EB030896](https://img.shields.io/badge/NIH_NIBIB-R01EB030896-green.svg)](https://grantome.com/grant/NIH/R01-EB030896-01)


#### MIT Copyright (c) 2021 brainlife.io The University of Texas at Austin and Indiana University
