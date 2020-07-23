# Plate Modeling

The notebooks in this directory enable a precision materials index (PMI) study
of a laminate angle-ply plate. Note that the notebooks must be run in a
particular order for successful operation:

1. `plate_doe.ipynb` generates realizations from the laminate structural response.
2. `plate_modeling.Rmd` uses these realizations to fit a monomial-standard model.
3. `plate_pmi.Rmd` uses the monomial-standard model to compute PMIs for different scenarios.

## Dependencies

- `plate_doe.ipynb` relies on the open-source [py_grama](https://github.com/zdelrosario/py_grama) package, as well as [Jupyter](https://jupyter.org/) to run the notebook.
  - I recommend [Anaconda](https://www.anaconda.com/products/individual) as a convenient way to install most Python dependencies.
- Both `.Rmd` notebooks rely on the [Tidyverse](https://www.tidyverse.org/) for data wrangling.
  - I recommend [Rstudio](https://rstudio.com/products/rstudio/) as a convenient way to work with Rmarkdown files.
