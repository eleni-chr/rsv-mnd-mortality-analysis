# RSV and Motor Neuron Disease Mortality Analysis

This repository contains seven R scripts for an ecological analysis of US monthly deaths with motor neuron disease (MND) as the underlying cause and national RSV and influenza laboratory surveillance. The analysis covers January 2011 through December 2019, with additional methodological checks and comparator outcomes.

The scripts are complementary analysis code for the associated report:

- **Report:** DOI forthcoming
- **Code repository:** DOI forthcoming
- **Analysis dataset:** DOI forthcoming

## What the scripts do

| Script | Purpose | Main output folder |
| --- | --- | --- |
| 01_ALS_virus_mortality_pilot.R | Feasibility analysis of monthly national MND mortality and influenza/RSV surveillance. | pilot_results/ |
| 02_ALS_virus_sensitivity_analysis.R | Sensitivity analyses using alternative RSV definitions, lags, and model specifications. | sensitivity_results/ |
| 03_RSV_technical_validation.R | Technical validation using Poisson models with HAC uncertainty estimates and model diagnostics. | technical_check_results/ |
| 04_RSV_MND_robustness_analysis.R | National comparator analyses and regional robustness analyses. | analysis_results/ |
| 05_RSV_MND_regional_validation.R | Regional validation, coverage checks, and paired regional resampling. | regional_validation_results/ |
| 06_RSV_MND_final_harmonised_validation.R | Harmonised national and regional comparisons on matched observations. | final_method_validation_results/ |
| 07_CDC_WONDER_RSV_MND_complementary_figures.R | Complementary figures and figure-data provenance files based on results from the preceding scripts. | complementary_figures/ |

Script 07 writes a figure-data manifest. If a final results table is unavailable or cannot be matched unambiguously, some Figure 2 estimates may come from rounded reference values embedded in the script; the manifest and plotting-data files record the source. Review these files against the analysis outputs before reusing the figures or values.

## Requirements and use

Install R and the packages used by the scripts. From an R session, the following installs the external packages used across the workflow:

    install.packages(c(
      "dplyr", "tidyr", "readr", "lubridate", "ggplot2", "MASS",
      "sandwich", "fixest", "patchwork", "ragg", "remotes"
    ))

The splines package is included with R. Script 01 can install some dependencies and attempt automated data retrieval when local files are absent; the automated requests may fail if a source service or its interface changes. The remotes package is used only if Script 01 needs to install the WONDER API helper.

Obtain the source data described below and place it in a raw_data/ folder at the repository root. Run scripts from that root directory, in numerical order:

    source("01_ALS_virus_mortality_pilot.R")
    source("02_ALS_virus_sensitivity_analysis.R")
    source("03_RSV_technical_validation.R")
    source("04_RSV_MND_robustness_analysis.R")
    source("05_RSV_MND_regional_validation.R")
    source("06_RSV_MND_final_harmonised_validation.R")
    source("07_CDC_WONDER_RSV_MND_complementary_figures.R")

Each later script uses output files from earlier scripts. Scripts 05 and 06 also require the HHS-region WONDER exports described below. Script 07 expects the generated results folders or their relevant CSV files to be available under the project root.

## Obtaining the source data

The repository does not include the source data. Create raw_data/ and retrieve the data directly from the public sources. Keep the original downloaded files and metadata unchanged so the scripts can check their contents and provenance.

### CDC WONDER mortality data

Use the [CDC WONDER Underlying Cause of Death query](https://wonder.cdc.gov/ucd-icd10.html) and its [API documentation](https://wonder.cdc.gov/wonder/help/wonder-api.html). The analysis uses underlying-cause-of-death ICD-10 selections, not multiple-cause selections:

- G12.2 — motor neuron disease (not ALS alone)
- G20 — Parkinson disease comparator
- V01–V99 — transport-injury comparator
- J09–J18 — influenza and pneumonia comparator

For each required outcome, select 2011–2019, all US residents, all ages, both sexes, and group by Year/Month. Export the complete tab-delimited WONDER result with its query metadata and footnotes. The primary national MND export should be saved as raw_data/wonder_G12_2_monthly.txt. National comparator exports used by Script 04 should be saved as:

- raw_data/Parkinson_G20_national_monthly.txt
- raw_data/Transport_national_monthly.txt
- raw_data/Respiratory_national_monthly.txt

For regional analyses, repeat the queries grouped by HHS Region and Year/Month and save the four exports with the following names:

- raw_data/MND_G12_2_HHS_monthly_2011_2019.txt
- raw_data/Parkinson_G20_HHS_monthly_2011_2019.txt
- raw_data/Transport_V01_V99_HHS_monthly_2011_2019.txt
- raw_data/Respiratory_J09_J18_HHS_monthly_2011_2019.txt

Script 01 can also use a two-column CSV named raw_data/wonder_G12_2_monthly.csv for the national MND series. The full tab-delimited export is preferred for the later scripts because they validate the query metadata. Script 04 can use a complete set of ten HHS-region monthly rows to derive a national total only when every region cell is observed. Script 05 requires exports grouped by HHS Region and Month; some filenames accepted by that script contain “national” for historical compatibility, but their contents must still be HHS-region tables.

CDC WONDER suppresses small counts. Preserve suppressed or absent cells as missing. Do not replace them with zero, infer them by subtracting regional totals from national totals, or redistribute them. Follow the [CDC WONDER data-use restrictions](https://wonder.cdc.gov/ucd-icd10.html) when using and redistributing these statistics.

### CDC RSV laboratory surveillance

Download the complete CSV from CDC's [Respiratory Syncytial Virus Laboratory Data (NREVSS)](https://data.cdc.gov/d/52kb-ccu2) dataset and save it as:

    raw_data/cdc_nrevss_rsv_52kb-ccu2.csv

The scripts use the weekly RSV tests and detections by HHS region and diagnostic method, including PCR. NREVSS is laboratory surveillance data and does not identify individual patients or deaths. Script 01 attempts to download this CSV automatically if it is not present.

### CDC influenza laboratory surveillance

Script 01 can retrieve national WHO/NREVSS influenza laboratory data through the cdcfluview R package. If retrieving it automatically is unavailable, download the data from CDC's [FluView Interactive surveillance portal](https://gis.cdc.gov/grasp/fluview/FluPortalDashboard.html) or the [FluView Interactive information page](https://www.cdc.gov/fluview/overview/fluview-interactive.html).

Use the ICL/NREVSS laboratory stream for the National area, not ILINet. The source is combined WHO/NREVSS laboratory reporting through the 2014–2015 season and clinical-laboratory reporting beginning in 2015–2016. Download both ranges, including 2010–2011 for the lagged exposure and 2019–2020 to cover the end of the study period. Save the original exports as:

- raw_data/WHO_NREVSS_Combined_prior_to_2015_16.csv
- raw_data/WHO_NREVSS_Clinical_Labs.csv

Alternatively, Script 01 accepts a harmonised weekly file named raw_data/cdc_fluview_weekly_harmonized.csv with the columns week_start, tests, positive, and source.

## Outputs and data release

The scripts create their results in the folders listed above, including tables, summaries, diagnostics, and figures. These generated files are not included in the code repository by default. The processed analysis data will be deposited separately:

- **Zenodo analysis dataset:** DOI forthcoming

The Apache 2.0 license in this repository applies to the R code only. It does not apply to CDC or other third-party source data, the separately deposited analysis dataset, or the report. Those materials remain subject to their own licenses and terms. Confirm that any data deposited or redistributed comply with the source terms and disclosure restrictions, including CDC WONDER suppression rules.

## Interpretation and disclaimer

This is exploratory ecological research. Group-level associations do not establish that individuals who died were infected with RSV or influenza, and they do not establish causation or individual risk. ICD-10 G12.2 covers motor neuron disease and is not specific to ALS. Surveillance methods, testing volume, reporting, confounding, and model assumptions may affect the estimates. The scripts and outputs are not medical advice or a clinical decision tool.

The software is provided “as is”, without warranties or conditions of any kind, express or implied, subject to the terms of the Apache License 2.0. No guarantee is made that the scripts will run unchanged with future versions of R, packages, APIs, or source data. See LICENSE for the complete warranty disclaimer and limitation of liability.

## Authors, copyright, and license

- **Author:** Dr Eleni Christoforidou
- **Copyright:** © 2026 Dr Eleni Christoforidou
- **Code license:** [Apache License 2.0](LICENSE) (Apache-2.0)

Third-party R packages are distributed under their respective licenses.

## Citation

If you use the code in this repository, please cite it using (DOI forthcoming).
