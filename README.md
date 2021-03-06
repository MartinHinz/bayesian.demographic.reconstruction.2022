
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesian.demographic.reconstruction.2022

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MartinHinz/bayesian.demographic.reconstruction.2022/master?urlpath=rstudio)
[![DOI](https://zenodo.org/badge/477622472.svg)](https://zenodo.org/badge/latestdoi/477622472)

This repository contains the data and code for our paper:

> Hinz, M., Roe, J., Laabs, J., Heitz, C., Kolar, J. (2022). *Bayesian
> inference of prehistoric population dynamics from multiple proxies: a
> case study from the North of the Swiss Alps*. Name of journal/book
> <https://doi.org/xxx/xxx>

Our pre-print is online here:

> Hinz, M., Roe, J., Laabs, J., Heitz, C., Kolar, J. (2022). *Bayesian
> inference of prehistoric population dynamics from multiple proxies: a
> case study from the North of the Swiss Alps*. Name of journal/book,
> Accessed 30 May 2022. Online at <https://doi.org/xxx/xxx>

### How to cite

Please cite this compendium as:

> Hinz, M., Roe, J., Laabs, J., Heitz, C., Kolar, J., (2022).
> *Compendium of R code and data for Bayesian inference of prehistoric
> population dynamics from multiple proxies: a case study from the North
> of the Swiss Alps*. Accessed 30 May 2022. Online at
> <https://doi.org/10.5281/zenodo.6594498>

## Authors

-   Martin Hinz (<martin.hinz@iaw.unibe.ch>), Institute of
    Archaeological Sciences & Oeschger Centre for Climate Change
    Research, University of Bern
    [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--9904--6548-green.svg)](http://orcid.org/0000-0002-9904-6548)
-   Joe Roe (<joe@joeroe.io>), Institute of Archaeological Sciences,
    University of Bern
    [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--1011--1244-green.svg)](http://orcid.org/0000-0002-1011-1244)
-   Julian Laabs (<julian.laabs@ufg.uni-kiel.de>), CRC 1266 - Scales of
    Transformation, University of Kiel
    [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0003--1378--3908-green.svg)](http://orcid.org/0000-0003-1378-3908)
-   Caroline Heitz (<caroline.heitz@iaw.unibe.ch>), CRC 1266 - Scales of
    Transformation, University of Kiel
    [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7188--6775-green.svg)](http://orcid.org/0000-0001-7188-6775)
-   Jan Kol???? (<jan.kolar@ibot.cas.cz>), Department of Vegetation
    Ecology, Institute of Botany of the Czech Academy of Sciences &
    Institute of Archaeology and Museology, Faculty of Arts, Masaryk
    University
    [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--8013--6992-green.svg)](http://orcid.org/0000-0001-8013-6992)

## Abstract

Robust estimates of population are essential to the study of
human???environment relations and socio-ecological dynamics in the past.
Population size and density can directly inform reconstructions of
prehistoric group size, social organisation, economic constraints,
exchange, and political and social institutions. In this pilot study, we
present an approach that we believe can be usefully transferred to other
regions, as well as refined and extended to greatly advance our
understanding of prehistoric demography.

Here, we present a Bayesian hierarchical model that uses Poisson
regression and state-space representation to produce absolute estimates
of past population size and density. Using the area North of the main
ridge of the Swiss Alps in prehistoric times (6000???1000 BCE) as a case
study, we show that combining multiple proxies (site counts, radiocarbon
dates, dendrochronological dates, and landscape openness) produces a
more robust reconstruction of population dynamics than any single proxy
alone. The model???s estimates of the credibility of its prediction, and
the relative weight it affords to individual proxies through time, give
further insights into the relative reliability of the evidence currently
available for paleodemographic research. Our prediction of population
development of the case study area accords well with the current
understanding in the wider literature, but provides a more precise and
higher-resolution estimate that is less sensitive to spurious
fluctuations in the proxy data than existing approaches, especially the
popular summed probability distribution of radiocarbon dates.

The archaeological record provides several potential proxies of human
population dynamics, but individually they are inaccurate, biased, and
sparse in their spatial and temporal coverage. Similarly, current
methods for estimating past population dynamics are often simplistic:
they work on limited spatial scales, tend to rely ona single proxy, and
are rarely able to infer population size or density in absolute terms.
In contemporary demography, it is becoming increasingly common to use
Bayesian statistics to estimate population trends and project them into
the future. The Bayesian approach is popular because offers the
possibility of combining heterogenous data, and at the same time
qualifying the uncertainty and credibility attached to forecasts. These
same characteristics make it well-suited to applications to
archaeological data in paleodemographic studies.

## Highlights

-   Bayesian modelling can integrate multiple, heterogeneous population
    proxies from the archaeological record
-   Our initial model produces more robust, high-resolution estimates of
    past population dynamics than previous, single-proxy approaches
-   We provide absolute estimates of population size and density on the
    area north of the Swiss Alpes in prehistoric times (6000???1000 BCE)

## Keywords

Prehistoric demography; Bayesian modelling; Multi-proxy; Settlement
dynamics

## Contents

The **repository** consists of:

-   [:file_folder: manuscript](/manuscript): R Markdown source document
    for manuscript.
-   [:file_folder: paper](/paper): A rendered version, of the submitted
    manuscript as
    [`paper_submission.docx`](/paper/paper_submission.docx) and
    [`paper_submission.pdf`](/paper/paper_submission.pdf), suitable for
    reading (the code is replaced by figures and tables in this file).
    You also find (as an editors cut) those versions before we trimmed
    them down to the word count of the journal
    ([`paper_extended_version.docx`](/paper/paper_extended_version.docx)
    and
    [`paper_extended_version.pdf`](/paper/paper_extended_version.pdf)
    respectively).
-   [:file_folder: analysis](/analysis): R Markdown source document for
    the actual analysis. Includes code to reproduce the figures and
    tables generated by the analysis. You can look at the
    [Github-Markdown Version](/analysis/analysis.md) for a rendered
    version of the analytical run incl.??figures.
-   [:file_folder: code](/code): Functions used in the analysis.
-   [:file_folder: data](/data): Data used in the analysis.
-   [:file_folder: figures](/figures): Plots and other illustrations

## How to run in your browser or download and run locally

This research compendium has been developed using the statistical
programming language R. To work with the compendium, you will need
installed on your computer the [R
software](https://cloud.r-project.org/) itself and optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

To perform the actual analysis, we recommend a powerful computer with a
multi-core processor, ideally with the Linux operating system, which has
at least 64GB RAM memory. In addition, a hard disk space of at least 5GB
should be reserved for the process.

### Licenses

**Text, code and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.
