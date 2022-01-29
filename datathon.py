#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 02:49:53 2022

@author: nikhaz
"""

import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon

geodata = gpd.read_file('/Users/nikhaz/Downloads/tl_2020_48_tabblock20/tl_2020_48_tabblock20.shp')
geo = geodata[geodata['COUNTYFP20'] == '201']
data = pd.read_csv('/Users/nikhaz/Downloads/nonlabeled.csv')
labels = pd.read_csv('/Users/nikhaz/Downloads/Block Labels v1.csv', names = ['Labels'])
data['label'] = labels['Labels']
data['GEOID'] = data['GEOID'].astype(str)
new = pd.merge(geo, data, left_on='GEOID20', right_on='GEOID')
new.explore(column = new['label'], categorical = True).save('blocks_map.html')