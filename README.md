[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4279402.svg)](https://doi.org/10.5281/zenodo.4279402)

# Mapping the languages of Twitter in Finland: richness and diversity in space and time 
This repository contains the scripts used in the article published in *Neuphilologische Mitteilungen* in 2020. The scripts are written in Python 3.7 using Anaconda distribution specifically.

### Requirements
* Python version 3.7 or newer (preferrably Anaconda distribution)
* A Twitter dataset containing user-specific timelines. To gather a similar data set, please see [Samuli Massinen's repo](https://github.com/DigitalGeographyLab/cross-border-mobility-twitter)
* Language identification model binary (lid.176.bin) for fastText available from here: https://fasttext.cc/docs/en/language-identification.html
* Area polygons. The administrative areas used in the article can be found from National Land Survey's [download service](https://tiedostopalvelu.maanmittauslaitos.fi/tp/kartta?lang=en) 

### How to use the scripts
The scripts can be run from the command line in a terminal.
1. Create a python environment by typing `conda create --name environment_name` and activate it `conda activate environment_name`
2. Instal packages using the [requirements.txt](requirements.txt) by typing `conda install --file requirements.txt`
3. Run scripts by typing `python scpritname.py` and add flags if required by the script

### Recommended running order of scripts
| Step | Script | Description | Input | Output |
| ---- | :----- | :---------- | :---- | :----- |
| 0 | [requirements.txt](requirements.txt) | Python environment text file | `conda install --file requirements.txt` | Python environment |
| 1 | [twitter_multilangid.py](twitter_multilangid.py) | Detects languages of tweets in database | database table | database table |
| 2 | [get_geotagged_posts.py](get_geotagged_posts.py) | Downloads geotagged tweets to disk | database table | Pickled dataframe |
| 3 | [get_user_langprofiles.py](get_user_lang_profiles.py) | Downloads user language profile to disk | Database table | Pickled dataframe |
| 4 | [post_location_detector.py](post_location_detector.py) | Detects location of geotagged posts | Output from step 1 & area polygons | Pickled dataframe |
| 5 | [detect_homes.py](detect_homes.py) | Detects likely home location of Twitter users using Samuli Massinen's HDA | Output from step 3 | Pickled dataframe and csv files |
| 6 | [user_language_diversity.py](user_language_diversity.py) | Calculates user language diversity metrics and diversity class | Output from step 2 & csv from step 4 | Pickled dataframe and pdf plot |
| 7 | [area_language_diversity.py](area_language_diversity.py) | Calculates language diversity in area polygons | Output from step 1 and area polygons | Polygon geopackage |
| 8 | [div_class_homes.py](div_class_homes.py) | Calculates home location distribution for low and moderate diversity users | Output from step 5 and home location csv from step 4 | Pickled dataframe |
| 9 | [temporal_plot.py](temporal_plot.py) | Plots temporal changes in language diversity across selected regions | Output from step 1 and areas polygon | PDF plot |
