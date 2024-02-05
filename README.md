# BaroclinicAdjustment

[![Build Status](https://github.com/simone-silvestri/BaroclinicAdjustment.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simone-silvestri/BaroclinicAdjustment.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://mit-license.org)

[![Documentation](https://img.shields.io/badge/documentation-stable%20release-red?style=flat-square)](https://simone-silvestri.github.io/BaroclinicAdjustment.jl/dev)

## Installation

To install BaroclinicAdjustment, follow these steps:

1. Clone the repository:

    ```shell
    git clone https://github.com/simone-silvestri/BaroclinicAdjustment.jl.git
    ```

2. Change to the project directory:

    ```shell
    cd BaroclinicAdjustment.jl
    ```

3. Install the required dependencies:

    ```shell
    julia --project -e 'import Pkg; Pkg.instantiate()'
    ```

## Running the Test Cases

To run the test cases, execute the following command from the project directory:

```shell
julia --project=run_and_visualize run_and_postprocess.jl
```
