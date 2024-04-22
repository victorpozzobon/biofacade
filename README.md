# Microalgae bio-reactive façade

Repository hosting the Python codes (please be kind with the typos in the comments) and data associated with the articles:
[Microalgae bio-reactive façade: A radiative-convective model powered by hourly illumination computation and historical weather data.](https://doi.org/10.1016/j.jobe.2024.109407) Pozzobon, V. (2024).  _Journal of Building Engineering_, 111352. [(Preprint)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_d.pdf)

[Microalgae bio-reactive façade: Location and weather-based systematic optimization.](https://doi.org/10.1016/j.buildenv.2024.111352) Pozzobon, V. (2024).  _Building and Environment_, 111352. [(Preprint)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_b.pdf)

It has been tested successfully on April 2024.

## How to run

__Step 1: execute the code__

Run _Biofacade_model_explicit.py_ (it will get the data from the _data_ folder) and generate two outputs:
- _results_vector_extended_Marignane_dt_0.5_azi_0.npy_ contains the outputs for all the time steps (several MB)
- _results_vector_Marignane_dt_0.5_azi_0.npy_ contains agglomerated results/indicators

__Step 2: sort the wavelengths by occurence__

You can plot the data with your favorite framework/software

![Image not found](./results/IllumYear.png?raw=true)
![Image not found](./results/TYear.png?raw=true)

## Data generation

Only sample data are distributed here. Indeed, [Météo France database](https://donneespubliques.meteofrance.fr/?fond=produit&id_produit=90&id_rubrique=32) represents more than 1.3 GB. You can find all the necessary tools to extract data from this base and turn it into input files for the code in the _helpFunctions_ folder (the routine is _Biofacade_preprocessing_).

## Contact

Please feel free to contact me. You can find my details on: https://victorpozzobon.github.io/
