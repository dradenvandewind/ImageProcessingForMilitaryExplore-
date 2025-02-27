
"""

# KEY = 'client_secret_705965637696-lpdn9g1tdeil70sq22sm7rjrklm17au7.apps.googleusercontent.com.json'
# PROJECT = 'plenary-agility-451209-b9'
# SERVICE_ACCOUNT='testimagesat@plenary-agility-451209-b9.iam.gserviceaccount.com'
# !gcloud auth login --project {PROJECT}
# !gcloud iam service-accounts keys create {KEY} --iam-account {SERVICE_ACCOUNT}


https://console.cloud.google.com/apis/credentials/wizard?api=earthengine.googleapis.com&previousPage=%2Fapis%2Fapi%2Fearthengine.googleapis.com%2Fmetrics%3Fproject%3Dplenary-agility-451209-b9%26inv%3D1%26invt%3DAbpzVw&inv=1&invt=AbpzVw&project=plenary-agility-451209-b9

705965637696-lpdn9g1tdeil70sq22sm7rjrklm17au7.apps.googleusercontent.com

gcloud auth login --project plenary-agility-451209-b9

gcloud iam service-accounts keys create client_secret_705965637696-lpdn9g1tdeil70sq22sm7rjrklm17au7.apps.googleusercontent.com.json  --iam-account testimagesat@plenary-agility-451209-b9.iam.gserviceaccount.com



"""

# Folium package provides a map widget.
import folium
import webbrowser
import json
import os

import requests
from PIL import Image
from io import BytesIO

from IPython import display
import matplotlib.pyplot as plt

import io
import rasterio
import numpy as np
from scipy.optimize import nnls  # Non-Negative Least Squares




KEY = 'client_secret_705965637696-lpdn9g1tdeil70sq22sm7rjrklm17au7.apps.googleusercontent.com.json'


def setup_folium(location, zoom=None, width=None, height=None):
  # Using folium.Figure instead of just folium.Map allows better control over
  # the map width and height.
  fig = folium.Figure(width=width, height=height)
  map = folium.Map(location=location, zoom_start=zoom)
  map.add_to(fig)
  return fig, map

from google.auth.transport.requests import AuthorizedSession
from google.oauth2 import service_account

credentials = service_account.Credentials.from_service_account_file(KEY)
scoped_credentials = credentials.with_scopes(
    ['https://www.googleapis.com/auth/cloud-platform'])
session = AuthorizedSession(scoped_credentials)
print(session)



EE_API = 'https://earthengine.googleapis.com/v1'
EE_PUBLIC = f'{EE_API}/projects/earthengine-public'


#gaza
"""
AOI_JSON = {
  "type": "Point",
  "coordinates": [31.40279742867142, 34.431847364499305]  # longitude, latitude
}
"""
AOI_JSON = {
  "type": "Point",
  "coordinates": [34.40279742867142 ,31.431847364499305]  # longitude, latitude
}




# GeoJSON coordinates are longitude and latitude, but folium accepts latitude 
# first and longitude second.
MAP_CENTER = list(reversed(AOI_JSON['coordinates']))
print('POI coordinates (lon, lat):', AOI_JSON['coordinates'])
print('Map center (lat, lon):', MAP_CENTER)

# Display the point on the map
folium_fig, folium_map = setup_folium(location=MAP_CENTER, zoom=25, width=1200, height=1200)
folium.Marker(MAP_CENTER, popup='FIIT STU').add_to(folium_map)
folium_fig

map_path_in = "gaza.html"
folium_map.save(map_path_in)

# Ouvrir automatiquement dans le navigateur
#webbrowser.open(map_path_in)

"""
# Search for the available images in catalog

Once we have the geometry (either a point, rectangle, or any kind of polygon) we can search through the catalog to get the available image data. Usually, we are not interested in all the images which intersect with the geometry. The most common image properties to filter are:
1. Dataset
2. Intersection with AOI
3. Time period when the image was captured
4. Cloud coverage (not needed in radar datasets)

"""

# Construct the filter
max_cloud_cover_pct = 30  # use 20-30 to select reasonably low cloud coverage
filters = [
    'startTime > "2024-07-01T00:00:00.000Z"',
    'endTime < "2025-09-01T00:00:00.000Z"',
    f'properties.CLOUDY_PIXEL_PERCENTAGE < {max_cloud_cover_pct}',
    f'intersects({json.dumps(json.dumps(AOI_JSON))})',
]
print('Filters:', *filters, sep='\n - ')

# Send the API request
query = {'filter': ' AND '.join(filters)}
url = f'{EE_PUBLIC}/assets/COPERNICUS/S2:listAssets'
response = session.get(url, params=query)
assets = response.json()
#print("response",response)
#print("assets",assets)


try:
    data = response.json()  # Convertir en dictionnaire
    print(data)

    error_code = data['error']['code']
    error_message = data['error']['message']
    error_status = data['error']['status']

    print(f"Erreur {error_code}: {error_status}")
    print(f"Détails: {error_message}")

except ValueError:
    print("La réponse n'est pas un JSON valide.")
except KeyError:
    print("Les clés d'erreur ne sont pas présentes dans la réponse.")


# Extract error details
if "error" in assets:

    """
    error_code = response["error"].get("code", "Unknown Code")
    error_message = response["error"].get("message", "No message provided")
    error_status = response["error"].get("status", "Unknown Status")

    print(f"Error Code: {error_code}")
    print(f"Message: {error_message}")
    print(f"Status: {error_status}")
    exit()
    """


assert(assets)
    


# List the filtered assets
print('Found image assets:')
cloud_covers = []
for asset in assets['assets']:
    id = asset['id']
    cloud_cover = asset['properties']['CLOUDY_PIXEL_PERCENTAGE']
    cloud_covers.append(cloud_cover)
    start_time = asset['startTime']
    print(f'{id} | {start_time} | {cloud_cover}')


 
ASSET_ID = 'COPERNICUS/S2/20250218T082031_20250218T082936_T36RXV'

# Send the API request
url = f'{EE_PUBLIC}/assets/{ASSET_ID}'
response = session.get(url)
asset = response.json()

print('Bands Count:', len(asset['bands']))
print('Band Names:', ', '.join(band['id'] for band in asset['bands']))
print('Second Band (B2, Blue):', json.dumps(asset['bands'][1], indent=2, sort_keys=True))

print('Geometry:', json.dumps(asset['geometry']))


# Show the asset geometry on the map
folium_fig, folium_map = setup_folium(location=MAP_CENTER, zoom=25, width=800, height=800)
folium.Marker(MAP_CENTER, popup='FIIT STU').add_to(folium_map)
folium.GeoJson(asset['geometry'], name='Asset geometry').add_to(folium_map)
#folium_fig

map_path = "map.html"
folium_map.save(map_path)

# Ouvrir automatiquement dans le navigateur
#webbrowser.open(map_path)


dimensions = asset['bands'][1]['grid']['dimensions']
print('Size:', round(dimensions['height'] * dimensions['width'] / 1e6, 2), 'Mpx')
print('Height:', round(dimensions['height'] * 10 / 1000, 2), 'km')
print('Width:', round(dimensions['width'] * 10 / 1000, 2), 'km')

url = f'{EE_PUBLIC}/assets/{ASSET_ID}:getPixels'
body = {
    'fileFormat': 'PNG',
    'bandIds': ['B4', 'B3', 'B2','B8', 'B11'],
    'region': asset['geometry'],
    'grid': {
        'dimensions': {'width': 1200, 'height': 1200},
    },
    'visualizationOptions': {
        #'ranges': [{'min': 0, 'max': 3000}],
        'ranges': [{'min': 1200, 'max': 3500}],
    },
}

#print('Geometry:', asset['geometry'])


image_response = session.post(url, json=body)
image_content = image_response.content

#display.Image(image_content)

image = Image.open(BytesIO(image_content))
plt.imshow(image)
plt.axis('off')  # Supprime les axes
#plt.show()

image_response = session.post(url, json=body)

print(image_response)

# Vérifier si la requête a réussi
if image_response.status_code == 200:
    image_content = image_response.content  # Obtenir les données de l'image

    # Charger et afficher l'image avec PIL
    image = Image.open(BytesIO(image_content))
    image.save("output_gaza.jpg")
    
    #image.show()  # Ouvre l'image avec la visionneuse d'images par défaut

else:
    print(f"Erreur {image_response.status_code}: {image_response.text}")


import pyproj  # Python interface to PROJ (https://proj.org/)
import shapely.geometry, shapely.ops

# Convert GeoJSON to shapefile geometry.
aoi_shp = shapely.geometry.shape(AOI_JSON)
print('Point in EPSG:4326 (degrees):', aoi_shp.wkt)

# Projection of the "standard" WGS/GPS CRS
WGS_PROJ = pyproj.Proj('EPSG:4326')
# Projection of the CRS of the image
IMG_PROJ = pyproj.Proj('EPSG:32633')

# Define the transformation between two projections:
transformer = pyproj.Transformer.from_proj(WGS_PROJ, IMG_PROJ, always_xy=True)

# Project the point in WGS84 CRS to the CRS of the image.
aoi_imcrs_shp = shapely.ops.transform(transformer.transform, aoi_shp)
print('Point in EPSG:32633 (meters):', aoi_imcrs_shp.wkt)

# Create rectangular polygon (a square) from the AOI in EPSG:32633
tile_size = 256 * 10  # 256 pixels x 10 m resolution
# Top-left corner of tile
xmin = aoi_imcrs_shp.x - tile_size / 2  # same as translate_x
ymin = aoi_imcrs_shp.y + tile_size / 2  # same as translate_y
# Bottom-right corner of the tile
xmax = xmin + tile_size
ymax = ymin - tile_size
# Helper function to generate square from the top-left and bottom-right coords.
tile_imgcrs_shp = shapely.geometry.box(xmin, ymin, xmax, ymax)
print('Tile in EPSG:32633:', tile_imgcrs_shp.wkt)
tile_imgcrs_shp

# Transformation from the image CRS to the WGS84 CRS
transformer = pyproj.Transformer.from_proj(IMG_PROJ, WGS_PROJ, always_xy=True)
wgs_tile_shp = shapely.ops.transform(transformer.transform, tile_imgcrs_shp)

print('Tile in EPSG:32633:', tile_imgcrs_shp.wkt)
print('Tile in EPSG:4326:', wgs_tile_shp.wkt)

# Convert the tile from shapefile to GeoJSON
wgs_tile_json = shapely.geometry.mapping(wgs_tile_shp)
print('Tile in EPSG:4326 (GEOJSON):', wgs_tile_json)

# Preview the tile on the map
folium_fig, folium_map = setup_folium(location=MAP_CENTER, zoom=18, width=1200, height=1200)
folium.Marker(MAP_CENTER, popup='FIIT STU').add_to(folium_map)
folium.GeoJson(wgs_tile_json, name='Tile area').add_to(folium_map)
folium_map.add_child(folium.LayerControl())
folium_fig

map_path_2 = "map_2.html"
folium_map.save(map_path_2)

# Ouvrir automatiquement dans le navigateur
#webbrowser.open(map_path_2)
#extract

# Set the region
#wgs_tile_json = {'type': 'Polygon', 'coordinates': (((-0.6199187254396578, 47.96739427270335), (-0.6130633970836579, 47.94521268231126), (-0.6460747622094035, 47.9406050114084), (-0.652943115266234, 47.96278303786113), (-0.6199187254396578, 47.96739427270335)),)}
#wgs_tile_json = {'type': 'Polygon', 'coordinates': (((34.40218625113778, 31.431946964975594), (34.40233463874827, 31.43625171617842), (34.40794638839432, 31.43633228450345), (34.407905919045646, 31.431728269130005), (34.402105312441364, 31.4319239443836)),)}

#wgs_tile_json = {'type': 'Polygon', 'coordinates': (((34.40218625113778, 31.431946964975594), (34.40233463874827, 31.43625171617842), (34.40794638839432, 31.43633228450345), (34.407905919045646, 31.431728269130005), (34.402105312441364, 31.4319239443836)),)}
wgs_tile_json = {
    "type": "Polygon",
    "coordinates": [
        [
            [34.3842009887357, 31.44476543638217],
            [34.41977052561313, 31.445092321824134],
            [34.41887649775694, 31.415177583130117],
            [34.380497159043216, 31.41566806563202],
            [34.36817234644448, 31.43572115362693],
            [34.384520284398434, 31.444874398322867],
            [34.38541431225579, 31.44492887924619],
            [34.3842009887357, 31.44476543638217]
        ]
    ]
}



# Preview the tile on the map
folium_fig, folium_map = setup_folium(location=MAP_CENTER, zoom=20, width=800, height=800)
folium.Marker(MAP_CENTER, popup='FIIT STU').add_to(folium_map)
folium.GeoJson(wgs_tile_json, name='Tile area').add_to(folium_map)
folium_map.add_child(folium.LayerControl())
folium_fig

url = f'{EE_PUBLIC}/assets/{ASSET_ID}:getPixels'
body = {
    'fileFormat': 'PNG',
    'bandIds': ['B4', 'B3', 'B2'],
    'region': wgs_tile_json,
    'grid': {
        'dimensions': {'width': 1600, 'height': 1600},
    },
    'visualizationOptions': {
        'ranges': [{'min': 200, 'max': 6500}],
    },
}
print('Region:', json.dumps(wgs_tile_json))
response = session.post(url, json=body)



try:
    data = response.json()  # Convertir en dictionnaire
    print(data)

    error_code = data['error']['code']
    error_message = data['error']['message']
    error_status = data['error']['status']

    print(f"Erreur {error_code}: {error_status}")
    print(f"Détails: {error_message}")

except ValueError:
    print("La réponse n'est pas un JSON valide.")
except KeyError:
    print("Les clés d'erreur ne sont pas présentes dans la réponse.")
    
    
print(response)

image_content = response.content
#display.Image(response.content)

image = Image.open(BytesIO(image_content))
#image = image.convert("RGB")
image.save("output_gaza_1.png")

import pyproj  # Python interface to PROJ (https://proj.org/)
import shapely.geometry, shapely.ops

# Convert GeoJSON to shapefile geometry.
aoi_shp = shapely.geometry.shape(AOI_JSON)
print('Point in EPSG:4326 (degrees):', aoi_shp.wkt)

# Projection of the "standard" WGS/GPS CRS
WGS_PROJ = pyproj.Proj('EPSG:4326')
# Projection of the CRS of the image
IMG_PROJ = pyproj.Proj('EPSG:32633')

# Define the transformation between two projections:
transformer = pyproj.Transformer.from_proj(WGS_PROJ, IMG_PROJ, always_xy=True)

# Project the point in WGS84 CRS to the CRS of the image.
aoi_imcrs_shp = shapely.ops.transform(transformer.transform, aoi_shp)
print('Point in EPSG:32633 (meters):', aoi_imcrs_shp.wkt)

# Create rectangular polygon (a square) from the AOI in EPSG:32633
tile_size = 256 * 10  # 256 pixels x 10 m resolution
# Top-left corner of tile
xmin = aoi_imcrs_shp.x - tile_size / 2  # same as translate_x
ymin = aoi_imcrs_shp.y + tile_size / 2  # same as translate_y
# Bottom-right corner of the tile
xmax = xmin + tile_size
ymax = ymin - tile_size
# Helper function to generate square from the top-left and bottom-right coords.
tile_imgcrs_shp = shapely.geometry.box(xmin, ymin, xmax, ymax)
print('Tile in EPSG:32633:', tile_imgcrs_shp.wkt)
tile_imgcrs_shp

# Transformation from the image CRS to the WGS84 CRS
transformer = pyproj.Transformer.from_proj(IMG_PROJ, WGS_PROJ, always_xy=True)
wgs_tile_shp = shapely.ops.transform(transformer.transform, tile_imgcrs_shp)

print('Tile in EPSG:32633:', tile_imgcrs_shp.wkt)
print('Tile in EPSG:4326:', wgs_tile_shp.wkt)

# Convert the tile from shapefile to GeoJSON
wgs_tile_json = shapely.geometry.mapping(wgs_tile_shp)
print('Tile in EPSG:4326 (GEOJSON):', wgs_tile_json)

# Preview the tile on the map
folium_fig, folium_map = setup_folium(location=MAP_CENTER, zoom=18, width=300, height=300)
folium.Marker(MAP_CENTER, popup='FIIT STU').add_to(folium_map)
folium.GeoJson(wgs_tile_json, name='Tile area').add_to(folium_map)
folium_map.add_child(folium.LayerControl())
folium_fig

map_path_2 = "map_2.html"
folium_map.save(map_path_2)

# Ouvrir automatiquement dans le navigateur
#webbrowser.open(map_path_2)


## Load the image data with rasterio


url = f'{EE_PUBLIC}/assets/{ASSET_ID}:getPixels'
body = {
    'fileFormat': 'GEO_TIFF',
    'bandIds': ['B2', 'B3', 'B4', 'B8', 'B11'],
    'region': wgs_tile_json,
    'grid': {
        'dimensions': {'width': 800, 'height': 800},
    },
}
response = session.post(url, json=body)
# Reading in memory data needs bit more massaging:
geotiff = rasterio.MemoryFile(io.BytesIO(response.content)).open()
type(geotiff)


print('Size:', f'{geotiff.width} x {geotiff.height}')
print('Number of bands:', geotiff.count)
print('Band indexes:', geotiff.indexes)
print('Bands:', dict(zip(['B2', 'B3', 'B4', 'B8', 'B11'], geotiff.indexes)))
print('CRS:', geotiff.crs)
print('Affine transformation:', geotiff.transform)

# The image is small, so we can preload all the pixels from the bands:
pixels = {
    'B2': geotiff.read(1),
    'B3': geotiff.read(2),
    'B4': geotiff.read(3),
    'B8': geotiff.read(4),
    'B11': geotiff.read(5)
}

# Investigate data in a single band:
print('Band data type:', type(pixels['B2']))
print('Band data shape:', pixels['B2'].shape)
print('Band data:', pixels['B2'])

import matplotlib
#matplotlib.use('TkAgg')  # Or 'Qt5Agg' if you have PyQt5 installed
matplotlib.use('Qt5Agg') 
import matplotlib.pyplot as pyplot
import rasterio
import rasterio.plot

# Create figure and subplots
fig, axes = pyplot.subplots(2, 3, figsize=(14, 14))
ax_iter = axes.flat

# Plot images
rasterio.plot.show(pixels['B4'], ax=next(ax_iter), cmap='Reds_r', title='B4 (Red)')
rasterio.plot.show(pixels['B8'], ax=next(ax_iter), cmap='Reds_r', title='B8 (Near infra-red)')
rasterio.plot.show(pixels['B3'], ax=next(ax_iter), cmap='Greens_r', title='B3 (Green)')
rasterio.plot.show(pixels['B2'], ax=next(ax_iter), cmap='Blues_r', title='B2 (Blue)')
rasterio.plot.show(pixels['B11'], ax=next(ax_iter), cmap='Oranges', title='B11 (Short-Wave Infrared))')

# Show the figure
pyplot.show()


# Make the band referencing easier
red, green, nir, swir = 'B4', 'B3', 'B8', 'B11'

# NDVI = (NIR-RED) / (NIR+RED)
ndvi = (pixels[nir].astype(float) - pixels[red].astype(float)) / (pixels[nir] + pixels[red])

# NDWI = (GREEN-NIR) / (GREEN+NIR)
ndwi = (pixels[green].astype(float) - pixels[nir].astype(float)) / (pixels[nir] + pixels[green])
#NDBI = SWIR -NIR / SIWR + NIR
swir_band = pixels['B11'].astype(float)  # SWIR (B11)
nir_band = pixels['B8'].astype(float)    # NIR (B8)

ndbi = (swir_band - nir_band) / (swir_band + nir_band)


# Display the results
fig, axes = pyplot.subplots(1, 3, figsize=(12, 6))
ax_iter = axes.flat
rasterio.plot.show(ndvi, ax=next(ax_iter), cmap='summer_r', title='NDVI')
rasterio.plot.show(ndwi, ax=next(ax_iter), cmap='winter_r', title='NDWI')
rasterio.plot.show(ndbi, ax=next(ax_iter), cmap='Oranges', title='NDBI')


#pyplot.show()



plt.figure(figsize=(10, 6))
plt.imshow(ndbi, cmap='Oranges', vmin=-1, vmax=1)
plt.colorbar(label="NDBI Value")
plt.title("NDBI (Normalized Difference Built-up Index)")
#plt.show()


green_band = pixels['B3'].astype(float)

"""
1. MNDWI (Modified Normalized Difference Water Index)

MNDWI is an improved version of NDWI, better suited for detecting water bodies in urban areas. It is often used in military operations to identify potential water hazards or map water bodies in urban settings, where vegetation may obstruct detection.
MNDWI=G−SWIRG+SWIR
MNDWI=G+SWIRG−SWIR​

MNDWI = (G - SWIR)/ (G + SWIR)

"""
print("compute MNDWI")

mndwi = (green_band - swir_band) / (green_band + swir_band)





"""
BAI (Burn Area Index)

BAI is used to detect burned areas, which can help in assessing damage after military strikes or natural wildfires. It is particularly useful for identifying scorched earth tactics or for disaster response.
BAI=1/ ((0.1 - REd)2 + (0.06 - SWIR)2)


    SWIR → B11
    Red → B4
    
print("compute BAI")

"""
red_band = pixels['B3'].astype(float)

bai =  1 / ((0.1 - red_band) ** 2 + (0.06 - swir_band) ** 2)

"""

 NIR / SWIR Ratio

This ratio is used to analyze vegetation and soil conditions. It can also help detect certain types of military vehicles or structures based on thermal and reflective properties. It is often used in combat zone analysis.
NIR/SWIR Ratio=NIRSWIR
NIR/SWIR Ratio=SWIRNIR​

    NIR → B8
    SWIR → B11

Use Case in Military:

    Detection of military vehicles or installations based on different reflectance in NIR and SWIR bands.
    Detecting changes in vegetation or land use after military activities.


"""

print("compute NIRSWIR")


rationir_swir = nir_band / swir_band

print("compute SWIRNIR​")

ratioswir_nir = swir_band / nir_band

"""
TCT is a transformation technique used to analyze vegetation, soil, and wetness in remote sensing data. It helps in terrain classification, and monitoring of military operations.

    Brightness (Overall light reflectance)
    Greenness (Vegetation-related)
    Wetness (Soil and water-related)


"""
print("compute tct")

url = f'{EE_PUBLIC}/assets/{ASSET_ID}:getPixels'
body = {
    'fileFormat': 'GEO_TIFF',
    'bandIds': ['B2', 'B3', 'B4', 'B5', 'B6','B7'],
    'region': wgs_tile_json,
    'grid': {
        'dimensions': {'width': 800, 'height': 800},
    },
}
response = session.post(url, json=body)
# Reading in memory data needs bit more massaging:
geotiff = rasterio.MemoryFile(io.BytesIO(response.content)).open()
type(geotiff)

print(" transformation Tasseled Cap ")
print('Size:', f'{geotiff.width} x {geotiff.height}')
print('Number of bands:', geotiff.count)
print('Band indexes:', geotiff.indexes)
print('Bands:', dict(zip(['B2', 'B3', 'B4', 'B5', 'B6', 'B7'], geotiff.indexes)))
print('CRS:', geotiff.crs)
print('Affine transformation:', geotiff.transform)

# The image is small, so we can preload all the pixels from the bands:
ttc_pixels = {
    'B2': geotiff.read(1),
    'B3': geotiff.read(2),
    'B4': geotiff.read(3),
    'B5': geotiff.read(4),
    'B6': geotiff.read(5),
    'B7': geotiff.read(6)
}


B2_band = ttc_pixels['B2'].astype(float)
B3_band = ttc_pixels['B3'].astype(float)
B4_band = ttc_pixels['B4'].astype(float)
B5_band = ttc_pixels['B5'].astype(float)
B6_band = ttc_pixels['B6'].astype(float)
B7_band = ttc_pixels['B7'].astype(float)


Br=0.3029*B2_band+0.2786*B3_band+0.4733*B4_band+0.5599*B5_band+0.5080*B6_band+0.1872*B7_band

Gr= -0.2941*B2_band - 0.243*B3_band - 0.5424*B4_band + 0.7276*B5_band + 0.071*B6_band -0.1608*B7_band

We = 0.1511*B2_band+0.1973*B3_band+0.3283*B4_band+0.3407*B5_band -0.7117*B6_band - 0.4559*B7_band






fig, axes = pyplot.subplots(2, 4, figsize=(12, 6))
ax_iter = axes.flat
rasterio.plot.show(ndvi, ax=next(ax_iter), cmap='summer_r', title='NDVI')
rasterio.plot.show(ndbi, ax=next(ax_iter), cmap='Oranges', title='NDBI')
rasterio.plot.show(mndwi, ax=next(ax_iter), cmap='summer_r', title='MNDWI')
rasterio.plot.show(bai, ax=next(ax_iter), cmap='RdYlBu', title='BAI')
rasterio.plot.show(rationir_swir, ax=next(ax_iter), cmap='turbo_r', title='rationir_swir')
rasterio.plot.show(ratioswir_nir, ax=next(ax_iter), cmap='terrain_r', title='ratioswir_nir')



pyplot.show()

fig, axes = pyplot.subplots(1, 3, figsize=(12, 6))
ax_iter = axes.flat
rasterio.plot.show(Br, ax=next(ax_iter), cmap='summer_r', title='TCT Brightness')
rasterio.plot.show(Gr, ax=next(ax_iter), cmap='Oranges', title='TCT Greenness')
rasterio.plot.show(We, ax=next(ax_iter), cmap='summer_r', title='TCT Wetness')



pyplot.show()

"""
 Spectral Mixture Analysis (SMA)
 """
 
print("compute Spectral Mixture Analysis")
 
 
url = f'{EE_PUBLIC}/assets/{ASSET_ID}:getPixels'
body = {
    'fileFormat': 'GEO_TIFF',
    'bandIds': [ 'B4', 'B8', 'B11'],
    'region': wgs_tile_json,
    'grid': {
        'dimensions': {'width': 800, 'height': 800},
    },
}
response = session.post(url, json=body)
# Reading in memory data needs bit more massaging:
geotiff = rasterio.MemoryFile(io.BytesIO(response.content)).open()
type(geotiff)

print(" transformation Tasseled Cap ")
print('Size:', f'{geotiff.width} x {geotiff.height}')
print('Number of bands:', geotiff.count)
print('Band indexes:', geotiff.indexes)
print('Bands:', dict(zip(['B4', 'B8', 'B11'], geotiff.indexes)))
print('CRS:', geotiff.crs)
print('Affine transformation:', geotiff.transform)

# The image is small, so we can preload all the pixels from the bands:
sma_pixels = {
    'B4': geotiff.read(1),
    'B8': geotiff.read(2),
    'B11': geotiff.read(3)
}



red = sma_pixels['B4'].astype(np.float32)
nir = sma_pixels['B8'].astype(np.float32)
swir = sma_pixels['B11'].astype(np.float32)

 
 # Définition des spectres des endmembers (sol nu, végétation, eau)
endmembers = np.array([
    [0.3, 0.4, 0.2],  # Sol nu (Red, NIR, SWIR)
    [0.05, 0.6, 0.1],  # Végétation (Red, NIR, SWIR)
    [0.02, 0.03, 0.5]  # Eau (Red, NIR, SWIR)
])

# Empiler les bandes pour chaque pixel
rows, cols = sma_pixels['B4'].shape
pixels = np.stack([red.ravel(), nir.ravel(), swir.ravel()], axis=1)

# Appliquer l'unmixing spectral (NNLS pour imposer f_i ≥ 0)
fractions = np.array([nnls(endmembers.T, pixel)[0] for pixel in pixels])

# Reshape en images 2D
fractions = fractions.reshape(rows, cols, -1)

# Sauvegarde des fractions en TIFF
print(type(sma_pixels['B4']))
print(dir(sma_pixels['B4']))
"""
from collections import namedtuple

RedBand = namedtuple("RedBand", ["data", "profile"])
red = RedBand(data=sma_pixels['B4'].astype(np.float32), profile=sma_pixels['B4'].profile)

# Accès
red_data = red.data
red_profile = red.profile



#profile = sma_pixels['B4'].profile
red_profile.update(count=3, dtype=rasterio.float32)
"""
"""
with rasterio.open("fractions.tif", "w", **red_profile) as dst:
    for i in range(3):
        dst.write(fractions[:, :, i], i + 1)
"""
# Définition des métadonnées (ici, le géoréférencement est fictif)
meta = {
    'driver': 'GTiff',  # Format GeoTIFF
    'count': 3,  # Nombre de bandes
    'dtype': 'float32',  # Type de données
    'crs': 'EPSG:4326',  # Système de référence spatiale (WGS84)
    'transform': rasterio.transform.from_origin(west=0, north=10, xsize=1, ysize=1),  # Géoréférencement fictif
    'width': 800,  # Largeur de l'image (en pixels)
    'height': 800  # Hauteur de l'image (en pixels)
}

# Ouvrir le fichier GeoTIFF d'entrée pour récupérer les métadonnées
with rasterio.open("output_fractions.tif", 'w', **meta) as dst:
        # Écrire chaque bande (sol nu, végétation, eau)
        dst.write(fractions[:, :, 0], 1)  # Bande 1: Sol nu
        dst.write(fractions[:, :, 1], 2)  # Bande 2: Végétation
        dst.write(fractions[:, :, 2], 3)  # Bande 3: Eau
        
#webbrowser.open(output_fractions.tif)

with rasterio.open('output_fractions.tif') as src:
    # Lire plusieurs bandes (par exemple, les bandes 1, 2 et 3)
    band1 = src.read(1)
    band2 = src.read(2)
    band3 = src.read(3)

# Afficher chaque bande séparément
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].imshow(band1, cmap='gray')
axes[0].set_title('Bande 1')
axes[1].imshow(band2, cmap='gray')
axes[1].set_title('Bande 2')
axes[2].imshow(band3, cmap='gray')
axes[2].set_title('Bande 3')

plt.tight_layout()
plt.show()
       
       
with rasterio.open('output_fractions.tif') as src:
    # Lire les bandes rouge (1), verte (2) et bleue (3)
    red = src.read(1)
    green = src.read(2)
    blue = src.read(3)

# Créer une image RGB en empilant les bandes
rgb = np.stack([red, green, blue], axis=-1)

# Afficher l'image RGB
plt.imshow(rgb)
plt.title('Affichage RGB (Bande 1: Rouge, Bande 2: Verte, Bande 3: Bleue)')
plt.show()

"""
Vegetation/Impervious Surface Fraconti (VIF)

This index measures the proportion of vegetation and impervious surfaces (e.g., roads, buildings) in an area. It's important for urban combat and civilian infrastructure monitoring in post-conflict zones.
Use Case in Military:

    Urban warfare operations: Analyzing vegetation and impervious surfaces in cities to assess line-of-sight, troop movement, and cover.
    Post-conflict assessments: Monitoring damage to infrastructure.
"""

print("compute Vegetation/Impervious Surface Fraconti (VIF)")

url = f'{EE_PUBLIC}/assets/{ASSET_ID}:getPixels'
body = {
    'fileFormat': 'GEO_TIFF',
    'bandIds': ['B4', 'B8', 'B11'],
    'region': wgs_tile_json,
    'grid': {
        'dimensions': {'width': 800, 'height': 800},
    },
}
response = session.post(url, json=body)
# Reading in memory data needs bit more massaging:
geotiff = rasterio.MemoryFile(io.BytesIO(response.content)).open()
type(geotiff)


print('Size:', f'{geotiff.width} x {geotiff.height}')
print('Number of bands:', geotiff.count)
print('Band indexes:', geotiff.indexes)
print('Bands:', dict(zip(['B4', 'B8', 'B11'], geotiff.indexes)))
print('CRS:', geotiff.crs)
print('Affine transformation:', geotiff.transform)

# The image is small, so we can preload all the pixels from the bands:
pixels = {
    'B4': geotiff.read(1),
    'B8': geotiff.read(2),
    'B11': geotiff.read(3)
}

# Investigate data in a single band:
print('Band data type:', type(pixels['B4']))
print('Band data shape:', pixels['B4'].shape)
print('Band data:', pixels['B4'])

"""
import rasterio.plot
from matplotlib import pyplot
"""


import rasterio
import rasterio.plot
import matplotlib.pyplot as pyplot



# Make the band referencing easier
red, nir, swir = 'B4', 'B8', 'B11'

# NDVI = (NIR-RED) / (NIR+RED)
ndvi = (pixels[nir].astype(float) - pixels[red].astype(float)) / (pixels[nir] + pixels[red])

#NDBI = SWIR -NIR / SIWR + NIR
swir_band = pixels['B11'].astype(float)  # SWIR (B11)
nir_band = pixels['B8'].astype(float)    # NIR (B8)
red_band = pixels['B4'].astype(float)    # RED (B4)

# Calcul du NDVI

ndbi = (swir_band - nir_band) / (swir_band + nir_band)

# Calcul du NDBI

ndvi = (nir_band - red_band) / (nir_band + red_band)

# Calcul du VIF
vif = ndvi - ndbi

# 'bandIds': ['B4', 'B8', 'B11'],
fig, axes = pyplot.subplots(2, 3, figsize=(14, 14))
ax_iter = axes.flat
rasterio.plot.show(pixels['B4'], ax=next(ax_iter), cmap='gray', title='B4 (Red)')
rasterio.plot.show(pixels['B8'], ax=next(ax_iter), cmap='gray', title='B8 (Near infra-red)')
rasterio.plot.show(pixels['B11'], ax=next(ax_iter), cmap='gray', title='B11')
ndbi = np.clip(ndbi, -1, 1)
ndvi = np.clip(ndvi, -1, 1)
vif = np.clip(vif, -1, 1)

rasterio.plot.show(ndbi, ax=next(ax_iter), cmap='gray', title='NDBI)')
rasterio.plot.show(ndvi, ax=next(ax_iter), cmap='gray', title='NVDI')
rasterio.plot.show(vif, ax=next(ax_iter), cmap='gray', title='VIF Vegetation/Impervious Surface Fraconti')



pyplot.show(block=True) 
