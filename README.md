# Microalgae bio-reactive façade

Repository hosting the Python codes (please be kind with the typos in the comments) and data associated with the articles:
- [Microalgae bio-reactive façade: A radiative-convective model powered by hourly illumination computation and historical weather data.](https://doi.org/10.1016/j.jobe.2024.109407) Pozzobon, V. (2024).  _Journal of Building Engineering_, 111352. [(Preprint)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_d.pdf)
- [Microalgae bio-reactive façade: Location and weather-based systematic optimization.](https://doi.org/10.1016/j.buildenv.2024.111352) Pozzobon, V. (2024).  _Building and Environment_, 111352. [(Preprint)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_b.pdf) [(Supplementary materials)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_b_Supplementary_Materials.pdf)
- [Microalgae bio-reactive façade: A model coupling weather, illumination, temperature, and cell growth over the year.](https://doi.org/10.1016/j.renene.2024.121545) **Pozzobon, V.** (2024).  _Renewable Energy_, 237B, 121545. [(Publisher open Access)](https://www.sciencedirect.com/science/article/pii/S0960148124014459/pdfft?md5=635ad712361fd1fd9967b9cbca3410c2&pid=1-s2.0-S0960148124014459-main.pdf) [(PDF file)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_g.pdf) [(Supplementary materials)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_g_Supplementary_Materials.pdf)
- [Microalgae bio-reactive façade: System thermal–biological optimization.](https://doi.org/10.1016/j.renene.2024.121377) **Pozzobon, V.** (2024).  _Renewable Energy_, 235, 121377. [(Publisher open Access)](https://www.sciencedirect.com/science/article/pii/S0960148124014459/pdfft?md5=635ad712361fd1fd9967b9cbca3410c2&pid=1-s2.0-S0960148124014459-main.pdf) [(PDF file)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_f.pdf) [(Supplementary materials)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_f_Supplementary_Materials.pdf)

It has been tested successfully on October 2024.

## Differences between the two codes

The first one is hosted in the _Temperature_and_light_only_ folder. This code predicts the temperature of the biofaçade module reservoir and the volume-averaged illumination within it. The second one is hosted in the repository's root folder (up there ↑). In addition to predicting temperature, it also predicts microalgae production and quality (pigment content). Both are described separately below.


It is supplied

## Temperature and light only

It is supplied with data for a South-oriented module located in Marseille. The input dataset span from 2013 to 2022. _illuminations_vector_Marignane_dt_0.5_azi_0.csv_ contains inicident illumination data. It has been generated using [Recommended Practice for the Calculation of Daylight Availability.](https://doi.org/10.1080/00994480.1984.10748791) (1984). _Journal of the Illuminating Engineering Society_, 13(4), implemented in _helpFunctions/Solar_Time.py_. _meteo_vector_Marignane_dt_0.5_azi_0.csv_ contains the relevant weather data extrcat from [Météo France database](https://donneespubliques.meteofrance.fr/?fond=produit&id_produit=90&id_rubrique=32) (see below).

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


## How to cite

Please cite the article presenting the complete biological model when reusing this work: [Microalgae bio-reactive façade: A model coupling weather, illumination, temperature, and cell growth over the year.](https://doi.org/10.1016/j.renene.2024.121545) **Pozzobon, V.** (2024).  _Renewable Energy_, 237B, 121545. [(Publisher open Access)](https://www.sciencedirect.com/science/article/pii/S0960148124014459/pdfft?md5=635ad712361fd1fd9967b9cbca3410c2&pid=1-s2.0-S0960148124014459-main.pdf) [(PDF file)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_g.pdf) [(Supplementary materials)](https://victorpozzobon.github.io/assets/preprints/Pozzobon_2024_g_Supplementary_Materials.pdf)

## Contact

Please feel free to contact me. You can find my details on: https://victorpozzobon.github.io/
