[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# CompEvo User Guide
## Comparative and Evolutionary Genomics Workflow

Snakemake workflow for comparative and evolutionary genomics. It is specially adapted to analyze the predictions generated [EXOGAP](https://github.com/dorinemerlat/EXOGAP).

## Authors
Dorine Merlat, iCube, University of Strasbourg, France, dorine.merlat@etu.unistra.fr

## What is CompEvo ?


## Getting Started
### Prerequisites

Les protéomes doivent avoir les en-têtes au même format que Uniprot, chaque protéome doit être placé dans un fichier. Pas d'isoformes, uniquement les formes canoniques.

### Installation
1. Install [SnakeMake](https://snakemake.readthedocs.io/en/stable/) following the procedure indicated [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

1. Clone the repository
    
    ```bash
    clone git https://github.com/dorinemerlat/CompEvo.git
    ```
    
## Usage

```bash
snakemake --use-conda --cores 64 [endfile]
```

## License

Distributed under the MIT License. See `LICENSE.txt` for more information.