import random

from shapely.geometry import Polygon, Point, LineString, MultiPoint
import matplotlib.pyplot as plt
import requests
import csv
import geopandas as gpd
import pandas as pd
from shapely.wkt import loads
from shapely.geometry import MultiPolygon


def get_response(lat, lng):
    headers = {
        'X-Auth-Token': 'sajY_zirdieU8d8WfsCn',
        'X-Auth-Email': 'finn.51@gmail.com',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) '
                      'Chrome/115.0.0.0 Safari/537.36 Edg/115.0.1901.183',
    }
    url = f'https://green.parcels.id.land/parcels/parcels/by_location.json?lng={lng}&lat={lat}'
    response = requests.get(url, headers=headers)
    return response.json()

new= get_response(35.17299611176094, -101.81457284471622)
print(new)


def find_polygon_corners(polygon):
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')
    for x, y in polygon.exterior.coords:
        min_x = min(min_x, x)
        min_y = min(min_y, y)
        max_x = max(max_x, x)
        max_y = max(max_y, y)
    corners = [
        (min_x, min_y),  # Bottom-left corner
        (max_x, min_y),  # Bottom-right corner
        (max_x, max_y),  # Top-right corner
        (min_x, max_y),  # Top-left corner
    ]
    return corners


def find_midpoints(corners):
    midpoints = []
    num_corners = len(corners)
    for i in range(num_corners):
        x_avg = (corners[i][0] + corners[(i + 1) % num_corners][0]) / 2
        y_avg = (corners[i][1] + corners[(i + 1) % num_corners][1]) / 2
        midpoints.append((x_avg, y_avg))
    return midpoints


def find_points_on_corners(polygon, distance, num_points=12):
    # Create buffer polygon
    buffer_polygon = polygon.buffer(distance)
    buffer_points = []
    # Calculate points along the boundary of the buffer polygon at equal intervals
    for i in range(0, num_points):
        point_on_boundary = buffer_polygon.boundary.interpolate(i / num_points, normalized=True)
        buffer_points.append((point_on_boundary.x, point_on_boundary.y))
    return buffer_points


# Polygon vertices
# polygon_vertices = [
#    (-101.830037876095, 35.178138626423),
#    (-101.829703132988, 35.1781256793867),
#    (-101.829701398771, 35.1778606290531),
#    (-101.829699606065, 35.1775865816554),
#    (-101.829697900717, 35.1773259302299),
#    (-101.829714544802, 35.1773121251595),
#    (-101.830035143174, 35.1773108352077),
#    (-101.830368101361, 35.177309372873),
#    (-101.830654799592, 35.1773081063686),
#    (-101.830671402242, 35.1773216327745),
#    (-101.830672006827, 35.1775825218913),
#    (-101.830672892886, 35.1781104163247),
#    (-101.830640462973, 35.1781371436679)

# ]

polygon_vertices = [
    (-101.813784875329, 35.172456066025),
    (-101.813679269762, 35.1724562068822),
    (-101.813661139675, 35.1724751238405),
    (-101.813644457815, 35.1724878396618),
    (-101.813646279369, 35.1727859681255),
    (-101.813786818382, 35.1727857777297),
    (-101.813784875329, 35.172456066025)
 ]

# polygon_vertices = [
#     (-101.79553987431, 35.3090191916404),
#     (-101.795486993006, 35.2945159538594),
#     (-101.782048905939, 35.2944826309803),
#     (-101.777676876069, 35.2944716474824),
#     (-101.777576416348, 35.2944713938042),
#     (-101.777636595966, 35.3090300262766),
#     (-101.795381030267, 35.3090284543452),
#     (-101.795539840892, 35.3090284222647),
#     (-101.79553987431, 35.3090191916404)
# ]

# Define your polygon and buffer distances
polygon = Polygon(polygon_vertices)
buffer_distances = [0.00010, 0.00020, 0.00100]


def create_buffer_points(polygon, buffer_distance, num_points=12):
    buffer_points = find_points_on_corners(polygon, buffer_distance, num_points)
    if not buffer_points:
        buffer_points = []  # Empty list if no points are found
    return pd.DataFrame(buffer_points, columns=['lng', 'lat'])


def scrape_data_within_buffer(buffer_points_df):
    parcel_responses = []
    for index, row in buffer_points_df.iterrows():
        lat = row['lat']
        lng = row['lng']
        response = get_response(lat, lng)
        parcel_responses.append(response)
    return parcel_responses


def create_and_scrape_buffer_zones(polygon, buffer_distances, num_points=12):
    for buffer_distance in buffer_distances:
        buffer_points_df = create_buffer_points(polygon, buffer_distance, num_points)
        parcel_responses = scrape_data_within_buffer(buffer_points_df)

        # Check if any response is not an error or empty
        valid_responses = [response for response in parcel_responses if response and ("errors" not in response)]

        if valid_responses:
            return buffer_points_df, parcel_responses
        else:
            print(
                f"No valid response found in buffer zone with distance {buffer_distance}. Trying the next buffer zone.")

    return None, []


# Create and scrape buffer zones
buffer_points_df, parcel_responses = create_and_scrape_buffer_zones(polygon, buffer_distances)

buffer_points_df['parcel_data'] = parcel_responses

# Create a GeoDataFrame with a geometry column containing Point objects
geometry = [Point(xy) for xy in zip(buffer_points_df['lng'], buffer_points_df['lat'])]
buffer_points_gdf = gpd.GeoDataFrame(buffer_points_df, geometry=geometry, crs='EPSG:4326')

if not buffer_points_gdf.empty:
    # Save the GeoDataFrame to a new CSV file with the parcel data and geometry
    buffer_points_gdf.to_csv('buffer_points_parcel_data1.csv', index=False)
else:
    print("No valid response found in any buffer zone. No CSV file will be created.")


def generate_12_buffer_points(random_point_geometry, buffer_distance):
    buffer_point_df = create_buffer_points(random_point_geometry, buffer_distance)
    parcel_responses = scrape_data_within_buffer(buffer_point_df)

    valid_responses = [response for response in parcel_responses if response and ("errors" not in response)]

    if valid_responses:
        return buffer_point_df.head(12).values.tolist(), parcel_responses[:12]
    else:
        print(f"No valid response found in buffer zone with distance {buffer_distance}.")
        return [], []

def extract_geometry_from_parcel_data(parcel_data):
    try:
        parcels = eval(parcel_data)['parcels']
        if parcels and 'parcel_data' in parcels[0]:
            geom_as_wkt = parcels[0]['parcel_data']['geom_as_wkt']
            return loads(geom_as_wkt)
        else:
            return None
    except Exception as e:
        print(f"Error extracting geometry: {e}")
        return None


def process_buffer_points(random_point_csv, buffer_distances, selected_buffer_distance):
    random_point_df = pd.read_csv(random_point_csv)

    latitude = random_point_df['lat'][0]
    longitude = random_point_df['lng'][0]

    random_point_geometry = Point(longitude, latitude)

    result_buffer_points, result_parcel_responses = generate_12_buffer_points(
        random_point_geometry, selected_buffer_distance)

    # Create a DataFrame from the result buffer points
    result_buffer_points_df = pd.DataFrame(result_buffer_points, columns=['lng', 'lat'])

    # Add the parcel responses list as a new column in the DataFrame
    result_buffer_points_df['parcel_data'] = result_parcel_responses

    # Create a GeoDataFrame with a geometry column containing Point objects
    geometry = [Point(xy) for xy in zip(result_buffer_points_df['lng'], result_buffer_points_df['lat'])]
    result_buffer_points_gdf = gpd.GeoDataFrame(result_buffer_points_df, geometry=geometry, crs='EPSG:4326')

    # Save the GeoDataFrame to a new CSV file with the parcel data
    result_buffer_points_gdf.to_csv('new_buffer_points.csv', index=False)

    new_buffer_points_df = pd.read_csv('new_buffer_points.csv')

    # Extract the new buffer points' coordinates as a set of tuples
    new_buffer_points_set = set(new_buffer_points_df[['lng', 'lat']].apply(tuple, axis=1))

    # Load the old buffer points from the CSV file
    old_buffer_points_df = pd.read_csv('buffer_points_parcel_data1.csv')

    # Extract the old buffer points' coordinates as a set of tuples
    old_buffer_points_set = set(old_buffer_points_df[['lng', 'lat']].apply(tuple, axis=1))

    plt.scatter(new_buffer_points_df['lng'], new_buffer_points_df['lat'], color='blue', label='Buffer Points')

    # Create a scatter plot for the old buffer points
    plt.scatter(old_buffer_points_df['lng'], old_buffer_points_df['lat'], color='red', label='Old Buffer Points')

    # Customize the plot if needed
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Buffer Points Comparison')
    plt.legend()
    plt.show()

def filter_new_buffer_points(result_buffer_points_df, old_buffer_points_set):

    old_buffer_points_df = pd.read_csv('buffer_points_parcel_data1.csv')
    old_buffer_points_df['geometry'] = old_buffer_points_df['parcel_data'].apply(extract_geometry_from_parcel_data)
    old_buffer_gdf = gpd.GeoDataFrame(old_buffer_points_df, geometry='geometry', crs='EPSG:4326')
    old_buffer_gdf = old_buffer_gdf.dropna(subset=['geometry'])  # Filter out None geometries

    new_buffer_points_df = pd.read_csv('new_buffer_points.csv')
    new_buffer_points_df['geometry'] = new_buffer_points_df['parcel_data'].apply(extract_geometry_from_parcel_data)
    new_buffer_gdf = gpd.GeoDataFrame(new_buffer_points_df, geometry='geometry', crs='EPSG:4326')
    new_buffer_gdf = new_buffer_gdf.dropna(subset=['geometry'])  # Filter out None geometries

# Filter out points with missing or invalid geometries
    new_buffer_points_filtered = new_buffer_gdf.dropna(subset=['geometry'])

# Use .within() to check if new buffer points are within old buffer points' geometries
    new_buffer_points_filtered = new_buffer_points_filtered[
        ~new_buffer_points_filtered['geometry'].apply(lambda geom: any(geom.within(old_geom) for old_geom in old_buffer_gdf['geometry']))]

# Save the filtered new buffer points without the 'geometry' column
    new_buffer_points_filtered.to_csv('new_buffer_points.csv', columns=['lng', 'lat', 'parcel_data'], index=False)

    duplicated_points = new_buffer_gdf[
        new_buffer_gdf['geometry'].apply(lambda geom: any(geom.within(old_geom) for old_geom in old_buffer_gdf['geometry']))]


# Save the non-duplicated points to a CSV file
    non_duplicated_points = new_buffer_gdf[
        ~new_buffer_gdf['geometry'].apply(lambda geom: any(geom.within(old_geom) for old_geom in old_buffer_gdf['geometry']))]

# Save the non-duplicated points to a CSV file
    non_duplicated_points.to_csv('non_duplicated_buffer_points.csv', columns=['lng', 'lat', 'parcel_data'], index=False)

# Save the duplicated points to a CSV file
    duplicated_points.to_csv('duplicated_buffer_points.csv', columns=['lng', 'lat', 'parcel_data'], index=False)


def generate_buffer_points_and_responses(selected_buffer_distance):
    # Load random point data
    random_point_df = pd.read_csv('new_point_and_response.csv')
    latitude = random_point_df['lat'][0]
    longitude = random_point_df['lng'][0]
    random_point_geometry = Point(longitude, latitude)

    # Generate buffer points and responses
    buffer_points, parcel_responses = generate_12_buffer_points(random_point_geometry, selected_buffer_distance)
    result_buffer_points_df = pd.DataFrame(buffer_points, columns=['lng', 'lat'])
    result_buffer_points_df['parcel_data'] = parcel_responses
    return result_buffer_points_df

def load_old_buffer_points():
    old_buffer_points_df = pd.read_csv('buffer_points_parcel_data1.csv')
    old_buffer_points_set = set(old_buffer_points_df[['lng', 'lat']].apply(tuple, axis=1))
    return old_buffer_points_df, old_buffer_points_set


def main():
    selected_buffer_distance = 0.00010
    # Generate buffer points and responses
    result_buffer_points_df = generate_buffer_points_and_responses(selected_buffer_distance)

    # Load old buffer points for comparison
    old_buffer_points_df, old_buffer_points_set = load_old_buffer_points()

    # Filter new buffer points that are not in the old set
    filter_new_buffer_points(result_buffer_points_df, old_buffer_points_set)

    # Call the function to process and visualize the buffer points
    process_buffer_points('new_point_and_response.csv', [0.00010, 0.00020, 0.00100], selected_buffer_distance)

if __name__ == "__main__":
    main()
