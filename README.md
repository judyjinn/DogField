# DogField

An experiment where search and rescues dogs are outfitted with GPS trackers and meters which track their head up and down position. We also take meteorological measurements over the course of a 0.5-1.0 mile trail. This script analyzes dog behavior and interprets how meteorological factors may affect dog trails deviating from the original trail laid by a human.

`field.py` is a script which manipulates GPS data from the arduino into a longitudinal and latitudinal coordinates to be graphed. It then pulls head up and down behavior and maps them over the trails. Finally, weather patterns are integreated to show temparature, wind, and humidity fluctuations throughout trails. Analyses are then performed to measure distance from dog trails to human trails.

Original data is not included from repository for privacy.
