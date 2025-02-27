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
import rasterio.plot
from datetime import datetime

import numpy as np
from scipy.optimize import nnls
import pyproj
import shapely.geometry, shapely.ops
from google.auth.transport.requests import AuthorizedSession
from google.oauth2 import service_account

# Constants
KEY = 'client_secret_705965637696-lpdn9g1tdeil70sq22sm7rjrklm17au7.apps.googleusercontent.com.json'
EE_API = 'https://earthengine.googleapis.com/v1'
EE_PUBLIC = f'{EE_API}/projects/earthengine-public'
"""
# AOI (Area of Interest)
AOI_JSON = {
    "type": "Point",
    "coordinates": [34.40279742867142, 31.431847364499305]  # longitude, latitude
}
# Polygon around AOI
TILE_JSON = {
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
"""

def calculate_centroid(polygon):
    x_coords = [point[0] for point in polygon]
    y_coords = [point[1] for point in polygon]
    centroid_x = sum(x_coords) / len(polygon)
    centroid_y = sum(y_coords) / len(polygon)
    return [centroid_x, centroid_y]

def format_point_coordinates(coords):
    if isinstance(coords[0], list):
        return coords[0]
    return coords


def extract_coordinates(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)
    
    coordinates_polygone = []
    coordinates_aio = []
    polygon_centroid = None
    
    if "features" in data:
        for feature in data["features"]:
            if "geometry" in feature and "coordinates" in feature["geometry"]:
                geom_type = feature["geometry"].get("type", "")
                coords = feature["geometry"]["coordinates"]
                
                # Filter only Point and LineString types
                if geom_type == "Point":
                    coordinates_aio.append(format_point_coordinates(coords))
                elif geom_type == "LineString":
                    # Ensure the LineString is closed
                    if len(coords) > 1 and coords[0] != coords[-1]:
                        print("Warning: LineString coordinates are not closed. Forcing closure.")
                        coords.append(coords[0])
                    coordinates_polygone.append(coords)
                    polygon_centroid = calculate_centroid(coords)  # Calculate centroid for LineString if needed
                elif geom_type == "Polygon":
                    polygon_centroid = calculate_centroid(coords[0])
    
    # If no Point is found, return centroid of the Polygon or LineString
    if not coordinates_aio and polygon_centroid:
        coordinates_aio.append(polygon_centroid)
    
    return coordinates_aio, coordinates_polygone


# Setup Folium Map
def setup_folium(location, zoom=None, width=None, height=None):
    fig = folium.Figure(width=width, height=height)
    map = folium.Map(location=location, zoom_start=zoom)
    map.add_to(fig)
    return fig, map

# Authenticate with Google Earth Engine
def authenticate_gee(key_path):
    credentials = service_account.Credentials.from_service_account_file(key_path)
    scoped_credentials = credentials.with_scopes(['https://www.googleapis.com/auth/cloud-platform'])
    session = AuthorizedSession(scoped_credentials)
    return session

# Fetch Assets from GEE
def fetch_assets(session, aoi_json, max_cloud_cover_pct=30):
    filters = [
        'startTime > "2024-07-01T00:00:00.000Z"',
        'endTime < "2025-09-01T00:00:00.000Z"',
        f'properties.CLOUDY_PIXEL_PERCENTAGE < {max_cloud_cover_pct}',
        f'intersects({json.dumps(json.dumps(aoi_json))})',
    ]
    query = {'filter': ' AND '.join(filters)}
    url = f'{EE_PUBLIC}/assets/COPERNICUS/S2:listAssets'
    response = session.get(url, params=query)
    return response.json()

def filtered_assets(assets_input):
    if "assets" not in assets_input or not isinstance(assets_input["assets"], list):
        print("Format incorrect des assets.")
        return None

    assets = []
    
    for asset in assets_input["assets"]:
        try:
            asset_id = asset["id"]
            start_time = asset["startTime"]
            cloud_cover = float(asset["properties"]["CLOUDY_PIXEL_PERCENTAGE"])

            assets.append({
                "id": asset_id,
                "start_time": datetime.fromisoformat(start_time.replace("Z", "")),  # Convertir en datetime
                "cloud_cover": cloud_cover
            })
        except (KeyError, ValueError) as e:
            print(f"Erreur lors du traitement d'un asset : {e}")
            continue

    if not assets:
        print("Aucun asset valide trouvé.")
        return None

    # Trier par couverture nuageuse croissante, puis par start_time décroissant
    best_asset = min(assets, key=lambda x: (x["cloud_cover"], -x["start_time"].timestamp()))

    print(f"Meilleur asset trouvé : {best_asset['id']} | {best_asset['start_time']} | {best_asset['cloud_cover']}% de couverture nuageuse")
    
    return best_asset["id"]




def create_json_tile(data_input):
    # Polygon around AOI
    TILE_JSON = {
        "type": "Polygon",
        "coordinates": data_input
    }
    return(TILE_JSON)
    
def create_json_aoi(data_input):
    # PoINT AOI
    
    AOI_JSON = {
    "type": "Point",
    "coordinates": data_input[0]  # longitude, latitude
    }
    return(AOI_JSON)
       

def fetch_image_data(session, asset_id, region, bands, width=800, height=800):
    url = f'{EE_PUBLIC}/assets/{asset_id}:getPixels'
    #print('Region:', json.dumps(region))
    
    body = {
        'fileFormat': 'GEO_TIFF',
        'bandIds': ['B3', 'B4', 'B8', 'B11'],
        'region': region,
        'grid': {
            'dimensions': {'width': width, 'height': height},
        },
    }
    response = session.post(url, json=body)
    return response.content


# Process Image Data
def process_image_data(image_content):
    geotiff = rasterio.MemoryFile(io.BytesIO(image_content)).open()
    type(geotiff)
    pixels = {
    'B3': geotiff.read(1),    
    'B4': geotiff.read(2),
    'B8': geotiff.read(3),
    'B11': geotiff.read(4)
    }

    #pixels = {band: geotiff.read(i+1) for i, band in enumerate(geotiff.indexes)}
    return geotiff, pixels

# Calculate Indices
def calculate_indices(pixels):
    red = pixels['B4'].astype(float)
    nir = pixels['B8'].astype(float)
    swir = pixels['B11'].astype(float)
    
    ndvi = (nir - red) / (nir + red)
    ndbi = (swir - nir) / (swir + nir)
    mndwi = (pixels['B3'].astype(float) - swir) / (pixels['B3'].astype(float) + swir)
    
    bai = 1 / ((0.1 - red) ** 2 + (0.06 - swir) ** 2)
    rationir_swir = nir / swir
    ratioswir_nir = swir / nir
    
    return ndvi, ndbi, mndwi, bai, rationir_swir, ratioswir_nir 

# Main Function
def main():
    input_from_geoson = "./input_from_geoson.io.json"
    #input_from_geoson = "./input_from_geoson_io_with_aoi.json"
    
    session = authenticate_gee(KEY)
    
    coordinates_aio, coordinates_polygon = extract_coordinates(input_from_geoson)
    """
    print(coordinates_polygon)
      
    print(coordinates_aio)
    """
   
    polygon_json = create_json_tile(coordinates_polygon)
    AOI_JSON = create_json_aoi(coordinates_aio)
     
    # Fetch Assets
    print("Fetch Assets")
    assets = fetch_assets(session, AOI_JSON)
    if "error" in assets:
        print(f"Error: {assets['error']['message']}")
        return
    
    ASSET_ID = filtered_assets(assets)
    assert(ASSET_ID)
    
    # Fetch Image Data
    print("Fetch Image Data")
    
    image_content = fetch_image_data(session, ASSET_ID, polygon_json, ['B3','B4', 'B8', 'B11'])
   
    

    # Process Image Data
    print("Process Image Data")

    geotiff, pixels = process_image_data(image_content)
    
    # Calculate Indices
    print("Calculate Indices")

    ndvi, ndbi, mndwi, bai, rationir_swir, ratioswir_nir = calculate_indices(pixels)
    #ndvi, ndbi, mndwi = calculate_indices(pixels)
    
    # Plot Results
    print("Plot Results")

    fig, axes = plt.subplots(2, 3, figsize=(14, 14))
    ax_iter = axes.flat
    rasterio.plot.show(ndvi, ax=next(ax_iter), cmap='gray', title='NDVI')
    rasterio.plot.show(ndbi, ax=next(ax_iter), cmap='gray', title='NDBI')
    rasterio.plot.show(mndwi, ax=next(ax_iter), cmap='gray', title='MNDWI')
    
    rasterio.plot.show(bai, ax=next(ax_iter), cmap='gray', title='BAI')
    rasterio.plot.show(rationir_swir, ax=next(ax_iter), cmap='gray', title='NIR/SWIR Ratio')
    rasterio.plot.show(ratioswir_nir, ax=next(ax_iter), cmap='gray', title='SWIR/NIR Ratio')
    
    
    plt.show()

if __name__ == "__main__":
    main()
