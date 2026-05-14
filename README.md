# solarr

## Author

- **Beniamino Sartini**
- **Email**: beniamino.sartini2@unibo.it

## Overview
**solarr** is an R package designed for implementing stochastic models and option pricing using solar radiation data. This package provides a range of tools and functions to analyze solar radiation, forecast its behavior, and evaluate options based on solar data.

## Installation

You can install the **solarr** package from GitHub using the following command:

```R
# Install devtools if you haven't already
install.packages("devtools")

# Install solarr from GitHub
devtools::install_github("beniamino98/solarr")
```

An online documentation is available [here](https://solarr-package.netlify.app/).

## Features

- Implementation of stochastic models for solar radiation.
- Option pricing tools specifically designed for solar data.
- Functions for seasonal adjustments and clear sky model estimations.
- Various statistical functions to analyze solar radiation data.

## Dependencies
The **solarr** package depends on the following R packages:
- `R (>= 4.4.0)`
- `ggplot2`
- `dplyr (>= 1.1.3)`
- `mclust`

Additionally, it imports the following packages:

- `broom`
- `formula.tools`
- `lifecycle`
- `lubridate (>= 1.8.0)`
- `mixtools`
- `numDeriv`
- `purrr (>= 1.0.2)`
- `rugarch (>= 1.4.1)`
- `stringr (>= 1.5.0)`
- `tidyr (>= 1.2.0)`

## License

This package is licensed under the GPL-3 License.

## Contributing

If you would like to contribute to **solarr**, please fork the repository and submit a pull request. 

## Issues

For any issues or feature requests, please open an issue on the GitHub repository.