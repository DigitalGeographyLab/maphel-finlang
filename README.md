# Mapping the languages of Twitter in Finland: richness and diversity in space and time 
A repository containing the scripts used in the article published in Neuphilologische Mitteilungen in 2020. The scripts use common python libraries

# Pre-requisites
* A Twitter dataset containing user-specific timelines. To gather this data set, please see Samuli Massinen's repo: https://github.com/DigitalGeographyLab/cross-border-mobility-twitter
* Language identification model binary for fastText available from here: https://fasttext.cc/docs/en/language-identification.html


### Recommended running order of scripts
| Step | Script | Description | Input | Output |
| ---- | :----- | :---------- | :---- | :----- |
| 0 | [twitter_multilangid.py](twitter_multilangid.py) | Detects languages of tweets in database | database table | database table |
| 1 | [get_geotagged_posts.py](get_geotagged_posts.py) | Downloads geotagged tweets to disk | database table | Pickled dataframe |
| 2 | [get_user_langprofiles.py](get_user_lang_profiles.py) | Downloads user language profile to disk | Database table | Pickled dataframe |
| 3 | [post_location_detector.py](post_location_detector.py) | Detects location of geotagged posts | Output from step 1 & area polygons | Pickled dataframe |
| 4 | [detect_homes.py](detect_homes.py) | Detects likely home location of Twitter users using Samuli Massinen's HDA | Output from step 3 | Pickled dataframe and csv files |
| 5 | [user_language_diversity.py](user_language_diversity.py) | Calculates user language diversity metrics and diversity class | Output from step 2 & csv from step 4 | Pickled dataframe and pdf plot |
| 6 | [area_language_diversity.py](area_language_diversity.py) | Calculates language diversity in area polygons | Output from step 1 and area polygons | Polygon geopackage |
| 7 | [div_class_homes.py](div_class_homes.py) | Calculates home location distribution for low and moderate diversity users | Output from step 5 and home location csv from step 4 | Pickled dataframe |
| 8 | [temporal_plot.py](temporal_plot.py) | Plots temporal changes in language diversity across selected regions | Output from step 1 and areas polygon | PDF plot |
