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
    #todo : add input argument validators
    x_coords = [point[0] for point in polygon]
    y_coords = [point[1] for point in polygon]
    centroid_x = sum(x_coords) / len(polygon)
    centroid_y = sum(y_coords) / len(polygon)
    return [centroid_x, centroid_y]

def format_point_coordinates(coords):
    #todo : add input argument validators
    if isinstance(coords[0], list):
        return coords[0]
    return coords


def extract_coordinates(json_file):
    #todo : add input argument validators
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
    #todo : add input argument validators
    fig = folium.Figure(width=width, height=height)
    map = folium.Map(location=location, zoom_start=zoom)
    map.add_to(fig)
    return fig, map

# Authenticate with Google Earth Engine
def authenticate_gee(key_path):
    #todo : add input argument validators
    credentials = service_account.Credentials.from_service_account_file(key_path)
    scoped_credentials = credentials.with_scopes(['https://www.googleapis.com/auth/cloud-platform'])
    session = AuthorizedSession(scoped_credentials)
    return session

# Fetch Assets from GEE
def fetch_assets(session, aoi_json, max_cloud_cover_pct=30):
    #todo : add input argument validators
    #rfc3339_date = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:000Z')
    #rfc3339_date = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S.000Z')
    #2025-09-01T00:00:00.000Z"
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
    #todo : add input argument validators
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

    print(f"Best asset found : {best_asset['id']} | {best_asset['start_time']} | {best_asset['cloud_cover']}% of cloud cover")
    
    return best_asset["id"]




def create_json_tile(data_input):
    #todo : add input argument validators
    # Polygon around AOI
    TILE_JSON = {
        "type": "Polygon",
        "coordinates": data_input
    }
    return(TILE_JSON)
    
def create_json_aoi(data_input):
    #todo : add input argument validators
    # PoINT AOI
    
    AOI_JSON = {
    "type": "Point",
    "coordinates": data_input[0]  # longitude, latitude
    }
    return(AOI_JSON)
       

def fetch_image_data(session, asset_id, region, bands, width=800, height=800):
    #todo : add input argument validators
    url = f'{EE_PUBLIC}/assets/{asset_id}:getPixels'
    #print('Region:', json.dumps(region))
    
    
    body = {
        'fileFormat': 'GEO_TIFF',
        'bandIds': bands,
        'region': region,
        'grid': {
            'dimensions': {'width': width, 'height': height},
        },
    }
    response = session.post(url, json=body)
    return response.content


# Process Image Data
def process_image_data(image_content, band):
    #todo : add input argument validators
    geotiff = rasterio.MemoryFile(io.BytesIO(image_content)).open()
    type(geotiff)
    print('Size:', f'{geotiff.width} x {geotiff.height}')
    print('Number of bands:', geotiff.count)
    print('Band indexes:', geotiff.indexes)
    #print('Bands:', dict(zip(['B2', 'B3', 'B4', 'B5', 'B6', 'B7'], geotiff.indexes)))
    print('Bands:', dict(zip(band, geotiff.indexes)))
    print('CRS:', geotiff.crs)
    print('Affine transformation:', geotiff.transform)
    pixels = {
    'B2': geotiff.read(1),
    'B3': geotiff.read(2),
    'B4': geotiff.read(3),
    'B5': geotiff.read(4),
    'B6': geotiff.read(5),
    'B7': geotiff.read(6),
    'B8': geotiff.read(7),
    'B11': geotiff.read(8)
    }

    #pixels = {band: geotiff.read(i+1) for i, band in enumerate(geotiff.indexes)}
    return geotiff, pixels

# Calculate Indices
def calculate_indices(pixels):
    #todo : add input argument validators
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
"""
TCT is a transformation technique used to analyze vegetation, soil, and wetness in remote sensing data. It helps in terrain classification, and monitoring of military operations.

    Brightness (Overall light reflectance)
    Greenness (Vegetation-related)
    Wetness (Soil and water-related)


"""

def compute_TCT(pixels):
    #todo : add input argument validators
    B2_band = pixels['B2'].astype(float)
    B3_band = pixels['B3'].astype(float)
    B4_band = pixels['B4'].astype(float)
    B5_band = pixels['B5'].astype(float)
    B6_band = pixels['B6'].astype(float)
    B7_band = pixels['B7'].astype(float)

    Br=0.3029*B2_band+0.2786*B3_band+0.4733*B4_band+0.5599*B5_band+0.5080*B6_band+0.1872*B7_band
    Gr= -0.2941*B2_band - 0.243*B3_band - 0.5424*B4_band + 0.7276*B5_band + 0.071*B6_band -0.1608*B7_band
    We = 0.1511*B2_band+0.1973*B3_band+0.3283*B4_band+0.3407*B5_band -0.7117*B6_band - 0.4559*B7_band
    return Br, Gr, We

def display_TCT(r_in,g_in,w_in):
    #todo : add input argument validators
    fig, axes = plt.subplots(1, 3, figsize=(12, 6))
    ax_iter = axes.flat
    rasterio.plot.show(r_in, ax=next(ax_iter), cmap='gray', title='TCT Brightness')
    rasterio.plot.show(g_in, ax=next(ax_iter), cmap='gray', title='TCT Greenness')
    rasterio.plot.show(w_in, ax=next(ax_iter), cmap='gray', title='TCT Wetness')
    plt.show()

"""
 Spectral Mixture Analysis (SMA)
"""
def compute_sma(in_pixels):
    #todo : add input argument validators
    red = in_pixels['B4'].astype(np.float32)
    nir = in_pixels['B8'].astype(np.float32)
    swir = in_pixels['B11'].astype(np.float32)
     # Définition des spectres des endmembers (sol nu, végétation, eau)
    endmembers = np.array([
    [0.3, 0.4, 0.2],  # Sol nu (Red, NIR, SWIR)
    [0.05, 0.6, 0.1],  # Végétation (Red, NIR, SWIR)
    [0.02, 0.03, 0.5]  # Eau (Red, NIR, SWIR)
    ])

    # Empiler les bandes pour chaque pixel
    rows, cols = in_pixels['B4'].shape
    new_pixels = np.stack([red.ravel(), nir.ravel(), swir.ravel()], axis=1)

    # Appliquer l'unmixing spectral (NNLS pour imposer f_i ≥ 0)
    fractions = np.array([nnls(endmembers.T, in_pixel)[0] for in_pixel in new_pixels])

    # Reshape en images 2D
    fractions = fractions.reshape(rows, cols, -1)
    return fractions
 
def save_sma(in_fractions, outfile_tif):
    #todo : add input argument validators
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
    with rasterio.open(outfile_tif, 'w', **meta) as dst:
        # Écrire chaque bande (sol nu, végétation, eau)
        dst.write(in_fractions[:, :, 0], 1)  # Bande 1: Sol nu
        dst.write(in_fractions[:, :, 1], 2)  # Bande 2: Végétation
        dst.write(in_fractions[:, :, 2], 3)  # Bande 3: Eau
    return True  
   
def display_tif_picture(input_file):
    #todo : add input argument validators
    
    with rasterio.open(input_file) as src:
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
       

def computeCenterMapCenter(aoi):
    #todo : add input argument validators
    MAP_CENTER = list(reversed(aoi['coordinates']))
    return MAP_CENTER
    

def show_aoi_map(session_id, ASSET_ID, map_center ):
    #todo : add input argument validators
    # Send the API request
    url = f'{EE_PUBLIC}/assets/{ASSET_ID}'
    response = session_id.get(url)
    asset = response.json()
    # Show the asset geometry on the map
    folium_fig, folium_map = setup_folium(location=map_center, zoom=15, width=800, height=800)
    folium.Marker(map_center, popup='FIIT STU').add_to(folium_map)
    folium.GeoJson(asset['geometry'], name='Asset geometry').add_to(folium_map)
    #folium_fig

    map_path = "map.html"
    folium_map.save(map_path)

    # Ouvrir automatiquement dans le navigateur
    webbrowser.open(map_path)


# Main Function
def main():
    #input_from_geoson = "./input_from_geoson.io.json"
    input_from_geoson = "./input_from_geoson_io_with_aoi.json"
    
    session = authenticate_gee(KEY)
    coordinates_aio, coordinates_polygon = extract_coordinates(input_from_geoson)
   
    polygon_json = create_json_tile(coordinates_polygon)
    AOI_JSON = create_json_aoi(coordinates_aio)
    
    MAP_CENTER = computeCenterMapCenter(AOI_JSON)
    print(MAP_CENTER)
     
    # Fetch Assets
    print("Fetch Assets")
    assets = fetch_assets(session, AOI_JSON)
    if "error" in assets:
        print(f"Error: {assets['error']['message']}")
        return
    #print(assets)
    ASSET_ID = filtered_assets(assets)
    assert(ASSET_ID)
    
    # Fetch Image Data
    print("Fetch Image Data")
    
    image_content = fetch_image_data(session, ASSET_ID, polygon_json, ['B2','B3','B4','B5','B6','B7','B8','B11'])
  
    # Process Image Data
    print("Process Image Data")

    geotiff, pixels = process_image_data(image_content,['B2','B3','B4','B5','B6','B7','B8','B11'])
    
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
    
    
    #plt.show()
    """
    show_aoi_map(session,ASSET_ID,MAP_CENTER)
    
    b_red, b_green, b_blue = compute_TCT(pixels)
    display_TCT(b_red, b_green, b_blue)
    
    
    r_fractions = compute_sma(pixels)
    save_sma(r_fractions, "./outfile_sma.tif")
    display_tif_picture("./outfile_sma.tif")
    """
    
    
    
    
    
    

if __name__ == "__main__":
    main()
