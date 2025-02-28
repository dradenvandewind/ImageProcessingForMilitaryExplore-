
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
KEY = 'client_secret_705965637696-lpdn9g1tdeil70sq22sm7rjrklm17au7.apps.googleusercontent.com.json'
EE_API = 'https://earthengine.googleapis.com/v1'
EE_PUBLIC = f'{EE_API}/projects/earthengine-public'

class GeoProcessor:
    def __init__(self, key_path):
        self.key_path = key_path
        self.session = self.authenticate_gee()

    def authenticate_gee(self):
        credentials = service_account.Credentials.from_service_account_file(self.key_path)
        scoped_credentials = credentials.with_scopes(['https://www.googleapis.com/auth/cloud-platform'])
        return AuthorizedSession(scoped_credentials)

    def extract_coordinates(self, json_file):
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
                    
                    if geom_type == "Point":
                        coordinates_aio.append(self.format_point_coordinates(coords))
                    elif geom_type == "LineString":
                        if len(coords) > 1 and coords[0] != coords[-1]:
                            print("Warning: LineString coordinates are not closed. Forcing closure.")
                            coords.append(coords[0])
                        coordinates_polygone.append(coords)
                        polygon_centroid = self.calculate_centroid(coords)
                    elif geom_type == "Polygon":
                        polygon_centroid = self.calculate_centroid(coords[0])
        
        if not coordinates_aio and polygon_centroid:
            coordinates_aio.append(polygon_centroid)
        
        return coordinates_aio, coordinates_polygone

    def calculate_centroid(self, polygon):
        #todo : add input argument validators
        x_coords = [point[0] for point in polygon]
        y_coords = [point[1] for point in polygon]
        centroid_x = sum(x_coords) / len(polygon)
        centroid_y = sum(y_coords) / len(polygon)
        return [centroid_x, centroid_y]

    def format_point_coordinates(self, coords):
        if isinstance(coords[0], list):
            return coords[0]
        return coords

    def create_json_tile(self, data_input):
        #todo : add input argument validators
        return {
            "type": "Polygon",
            "coordinates": data_input
        }
    
    def create_json_aoi(self, data_input):
        #todo : add input argument validators
        return {
            "type": "Point",
            "coordinates": data_input[0]
        }

    def fetch_assets(self, aoi_json, max_cloud_cover_pct=30):
        #todo : add input argument validators
        filters = [
            'startTime > "2024-07-01T00:00:00.000Z"',
            'endTime < "2025-09-01T00:00:00.000Z"',
            f'properties.CLOUDY_PIXEL_PERCENTAGE < {max_cloud_cover_pct}',
            f'intersects({json.dumps(json.dumps(aoi_json))})',
        ]
        query = {'filter': ' AND '.join(filters)}
        url = f'{EE_PUBLIC}/assets/COPERNICUS/S2:listAssets'
        response = self.session.get(url, params=query)
        return response.json()

    def filtered_assets(self, assets_input):
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
                    "start_time": datetime.fromisoformat(start_time.replace("Z", "")),
                    "cloud_cover": cloud_cover
                })
            except (KeyError, ValueError) as e:
                print(f"Erreur lors du traitement d'un asset : {e}")
                continue

        if not assets:
            print("Aucun asset valide trouvé.")
            return None

        best_asset = min(assets, key=lambda x: (x["cloud_cover"], -x["start_time"].timestamp()))
        print(f"Best asset found : {best_asset['id']} | {best_asset['start_time']} | {best_asset['cloud_cover']}% of cloud cover")
        
        return best_asset["id"]

    def fetch_image_data(self, asset_id, region, bands, width=800, height=800):
        #todo : add input argument validators
        url = f'{EE_PUBLIC}/assets/{asset_id}:getPixels'
        body = {
            'fileFormat': 'GEO_TIFF',
            'bandIds': bands,
            'region': region,
            'grid': {
                'dimensions': {'width': width, 'height': height},
            },
        }
        response = self.session.post(url, json=body)
        return response.content

    def process_image_data(self, image_content, band):
        #todo : add input argument validators
        geotiff = rasterio.MemoryFile(io.BytesIO(image_content)).open()
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
        return geotiff, pixels

    def calculate_indices(self, pixels):
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

    def compute_TCT(self, pixels):
        #todo : add input argument validators
        B2_band = pixels['B2'].astype(float)
        B3_band = pixels['B3'].astype(float)
        B4_band = pixels['B4'].astype(float)
        B5_band = pixels['B5'].astype(float)
        B6_band = pixels['B6'].astype(float)
        B7_band = pixels['B7'].astype(float)

        Br = 0.3029*B2_band + 0.2786*B3_band + 0.4733*B4_band + 0.5599*B5_band + 0.5080*B6_band + 0.1872*B7_band
        Gr = -0.2941*B2_band - 0.243*B3_band - 0.5424*B4_band + 0.7276*B5_band + 0.071*B6_band - 0.1608*B7_band
        We = 0.1511*B2_band + 0.1973*B3_band + 0.3283*B4_band + 0.3407*B5_band - 0.7117*B6_band - 0.4559*B7_band
        return Br, Gr, We

    def display_TCT(self, r_in, g_in, w_in):
        #todo : add input argument validators
        fig, axes = plt.subplots(1, 3, figsize=(12, 6))
        ax_iter = axes.flat
        rasterio.plot.show(r_in, ax=next(ax_iter), cmap='gray', title='TCT Brightness')
        rasterio.plot.show(g_in, ax=next(ax_iter), cmap='gray', title='TCT Greenness')
        rasterio.plot.show(w_in, ax=next(ax_iter), cmap='gray', title='TCT Wetness')
        plt.show()

    def compute_sma(self, in_pixels):
        #todo : add input argument validators
        red = in_pixels['B4'].astype(np.float32)
        nir = in_pixels['B8'].astype(np.float32)
        swir = in_pixels['B11'].astype(np.float32)
        endmembers = np.array([
            [0.3, 0.4, 0.2],
            [0.05, 0.6, 0.1],
            [0.02, 0.03, 0.5]
        ])
        rows, cols = in_pixels['B4'].shape
        new_pixels = np.stack([red.ravel(), nir.ravel(), swir.ravel()], axis=1)
        fractions = np.array([nnls(endmembers.T, in_pixel)[0] for in_pixel in new_pixels])
        fractions = fractions.reshape(rows, cols, -1)
        return fractions

    def save_sma(self, in_fractions, outfile_tif):
        #todo : add input argument validators
        meta = {
            'driver': 'GTiff',
            'count': 3,
            'dtype': 'float32',
            'crs': 'EPSG:4326',
            'transform': rasterio.transform.from_origin(west=0, north=10, xsize=1, ysize=1),
            'width': 800,
            'height': 800
        }
        with rasterio.open(outfile_tif, 'w', **meta) as dst:
            dst.write(in_fractions[:, :, 0], 1)
            dst.write(in_fractions[:, :, 1], 2)
            dst.write(in_fractions[:, :, 2], 3)
        return True

    def display_tif_picture_with_3_bands(self, input_file):
        #todo : add input argument validators
        with rasterio.open(input_file) as src:
            band1 = src.read(1)
            band2 = src.read(2)
            band3 = src.read(3)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        axes[0].imshow(band1, cmap='gray')
        axes[0].set_title('Bande 1')
        axes[1].imshow(band2, cmap='gray')
        axes[1].set_title('Bande 2')
        axes[2].imshow(band3, cmap='gray')
        axes[2].set_title('Bande 3')
        plt.tight_layout()
        plt.show()

    def compute_VIF(self, in_pixels):
        #todo : add input argument validators
        red, nir, swir = 'B4', 'B8', 'B11'
        ndvi = (in_pixels[nir].astype(float) - in_pixels[red].astype(float)) / (in_pixels[nir] + in_pixels[red])
        swir_band = in_pixels['B11'].astype(float)
        nir_band = in_pixels['B8'].astype(float)
        red_band = in_pixels['B4'].astype(float)
        ndbi = (swir_band - nir_band) / (swir_band + nir_band)
        ndvi = (nir_band - red_band) / (nir_band + red_band)
        vif = ndvi - ndbi
        return vif, ndvi, ndbi

    def computeCenterMapCenter(self, aoi):
        #todo : add input argument validators
        return list(reversed(aoi['coordinates']))

    def show_aoi_map(self, ASSET_ID, map_center):
        #todo : add input argument validators
        url = f'{EE_PUBLIC}/assets/{ASSET_ID}'
        response = self.session.get(url)
        asset = response.json()
        folium_fig, folium_map = self.setup_folium(location=map_center, zoom=15, width=800, height=800)
        folium.Marker(map_center, popup='FIIT STU').add_to(folium_map)
        folium.GeoJson(asset['geometry'], name='Asset geometry').add_to(folium_map)
        map_path = "map.html"
        folium_map.save(map_path)
        webbrowser.open(map_path)

    def setup_folium(self, location, zoom=None, width=None, height=None):
        #todo : add input argument validators
        fig = folium.Figure(width=width, height=height)
        map = folium.Map(location=location, zoom_start=zoom)
        map.add_to(fig)
        return fig, map

def main():
    input_from_geoson = "./input_from_geoson_io_with_aoi.json"
    processor = GeoProcessor(KEY)
    coordinates_aio, coordinates_polygon = processor.extract_coordinates(input_from_geoson)
    polygon_json = processor.create_json_tile(coordinates_polygon)
    AOI_JSON = processor.create_json_aoi(coordinates_aio)
    MAP_CENTER = processor.computeCenterMapCenter(AOI_JSON)
    print(MAP_CENTER)
    assets = processor.fetch_assets(AOI_JSON)
    if "error" in assets:
        print(f"Error: {assets['error']['message']}")
        return
    ASSET_ID = processor.filtered_assets(assets)
    assert(ASSET_ID)
    image_content = processor.fetch_image_data(ASSET_ID, polygon_json, ['B2','B3','B4','B5','B6','B7','B8','B11'])
    geotiff, pixels = processor.process_image_data(image_content,['B2','B3','B4','B5','B6','B7','B8','B11'])
    ndvi, ndbi, mndwi, bai, rationir_swir, ratioswir_nir = processor.calculate_indices(pixels)
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