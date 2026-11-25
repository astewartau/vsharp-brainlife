<h1 align="center">V-SHARP Background Field Removal</h1>

<p align="center">
  Removes the background field from a QSM total field map with V-SHARP (QSM.jl in Julia), producing a tissue (local) field map and the eroded V-SHARP mask.
</p>

<p align="center">
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg" alt="ABCD Spec v1.1"></a>
  <a href="https://doi.org/10.25663/brainlife.app.835"><img src="https://img.shields.io/badge/Brainlife-brainlife.app.835-blue.svg" alt="Run on brainlife.io"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-see_LICENSE-lightgrey.svg" alt="License"></a>
</p>

## At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Removes the background field from a QSM total field map with V-SHARP (QSM.jl in Julia), producing a tissue (local) field map and the eroded V-SHARP mask.**

Installs Julia 1.9.3 and the QSM.jl, NIfTI.jl and MriResearchTools packages at run time, masks the field map, converts it to Hz using the field map units (from input_units or the fieldmap JSON) and the field strength, and runs V-SHARP background field removal. Use it between phase unwrapping/field estimation and dipole inversion in a QSM pipeline. Needs internet access on the compute node; the README is an unfilled template.

| | |
|---|---|
| **Inputs** | `neuro/fmap` · `neuro/mask` |
| **Outputs** | `neuro/fmap` · `neuro/mask` |
| **Resources** | 16 GB memory · 1 h walltime |
| **Run it** | [doi.org/10.25663/brainlife.app.835](https://doi.org/10.25663/brainlife.app.835) |

#### Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `fieldmap` | `neuro/fmap` (single) | yes | Total B0 field map with magnitude and JSON sidecars. Field strength must be in the input metadata (MagneticFieldStrength) or B0; units must be in the fieldmap JSON (Units) or input_units. |
| `mask` | `neuro/mask` | yes | Brain mask restricting the field map. |

#### Outputs

| Output | Datatype | Description |
|---|---|---|
| `fmap` | `neuro/fmap` | Magnetic field variations caused by tissues (background field removed). |
| `mask` | `neuro/mask` | V-SHARP mask |

<!-- brainlife-app:end -->

## Running the App

### On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.835](https://doi.org/10.25663/brainlife.app.835), choose `neuro/fmap`, `neuro/mask` from your project and press **Execute**. The results are saved to your project as `neuro/fmap`, `neuro/mask`.

### Locally

You need [`jq`](https://jqlang.github.io/jq/).

1. Clone this repository.
2. Create `config.json` with the paths to your inputs and the parameters:

   ```json
   {
       "fieldmap": "/path/to/fieldmap/fieldmap",
       "fieldmap_json": "/path/to/fieldmap/fieldmap_json",
       "magnitude": "/path/to/fieldmap/magnitude",
       "magnitude_json": "/path/to/fieldmap/magnitude_json",
       "mask": "/path/to/mask/mask"
   }
   ```

3. Run `./main`.

## Troubleshooting

- **Values for either 'input_units' or 'fieldmap_json' are required!**: Provide a fieldmap JSON with a Units field.
- **MagneticFieldStrength not found in _inputs for fieldmap**: Set MagneticFieldStrength in the field map metadata.

## Citing

If you use this App, please cite it, the software it builds on and brainlife.io:

- V-SHARP Background Field Removal. brainlife.io App. [doi:10.25663/brainlife.app.835](https://doi.org/10.25663/brainlife.app.835)
- **Add any other citations based on the app**
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## Authors

- [Author1](email), Affiliation

**Contributors**

- [Author2](email), Affiliation
- [Author3](email), Affiliation

## Funding

[![Wellcome-226486/Z/22/Z](https://img.shields.io/badge/Wellcome-226486%2FZ%2F22%2FZ-purple.svg)](https://wellcome.org/)
[![NIH-NINDS-U24NS140384](https://img.shields.io/badge/NIH_NINDS-U24NS140384-green.svg)](https://reporter.nih.gov/project-details/11261037)
[![NIH-NINDS-UM1NS132207](https://img.shields.io/badge/NIH_NINDS-UM1NS132207-green.svg)](https://reporter.nih.gov/project-details/11411328)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://reporter.nih.gov/project-details/10262925)
[![NIH-NIBIB-R01EB030896](https://img.shields.io/badge/NIH_NIBIB-R01EB030896-green.svg)](https://reporter.nih.gov/project-details/10428625)
[![NIH-NIMH-R01MH126699](https://img.shields.io/badge/NIH_NIMH-R01MH126699-green.svg)](https://reporter.nih.gov/project-details/10460628)
[![NIH-NIDCD-R01DC020088](https://img.shields.io/badge/NIH_NIDCD-R01DC020088-green.svg)](https://reporter.nih.gov/project-details/11140284)
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-2203524](https://img.shields.io/badge/NSF_BCS-2203524-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2203524)
[![NSF-IIS-1636893](https://img.shields.io/badge/NSF_IIS-1636893-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-OAC-2148729](https://img.shields.io/badge/NSF_OAC-2148729-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148729)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NSF-IIS-2148700](https://img.shields.io/badge/NSF_IIS-2148700-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148700)
[![NSF-OAC-1916518](https://img.shields.io/badge/NSF_OAC-1916518-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1916518)

<sub>brainlife.io, which runs this App, is publicly funded by Wellcome, the U.S. National Institutes of Health and the U.S. National Science Foundation; acknowledging that support helps keep it running.</sub>

## License

See [LICENSE](LICENSE).
