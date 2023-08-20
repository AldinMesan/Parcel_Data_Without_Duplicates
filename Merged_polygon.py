import math

from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import requests
import pandas as pd
import geopandas as gpd
from shapely.ops import nearest_points
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

# Create a shapely Polygon object
polygon = Polygon(polygon_vertices)
# Desired buffer distance
buffer_distance = 0.0001
# Find the four corners of the polygon
corners = find_polygon_corners(polygon)
# Find points on corners of the buffer polygon
buffer_points = find_points_on_corners(polygon, buffer_distance, num_points=12)
# print("Buffer Points:", buffer_points)
# Convert buffer points to a DataFrame
df_buffer_points = pd.DataFrame(buffer_points, columns=['lng', 'lat'])
# Save the DataFrame to a CSV file
df_buffer_points.to_csv('buffer_points.csv', index=False)
# Read the CSV file into a DataFrame
df_coordinates = pd.read_csv('buffer_points.csv')
# Initialize lists to store the scraped data
parcel_responses = []
for index, row in df_coordinates.iterrows():
    lat = row['lat']
    lng = row['lng']
    response = get_response(lat, lng)
    parcel_responses.append(response)
# Add the parcel_responses list as a new column in the DataFrame
df_coordinates['parcel_data'] = parcel_responses
# Save the DataFrame to a new CSV file with the parcel data
df_coordinates.to_csv('buffer_points_parcel_data.csv', index=False)

# Initialize an empty MultiPolygon to store the merged result
merged_polygon = MultiPolygon([polygon])
# Iterate through the DataFrame and merge each polygon with the original one
for index, row in df_coordinates.iterrows():
    lat = row['lat']
    lng = row['lng']
    response = get_response(lat, lng)
    if 'parcels' in response and len(response['parcels']) > 0:
        parcel_data = response['parcels'][0]['parcel_data']
        parcel_polygon_wkt = parcel_data['geom_as_wkt']
        parcel_polygon = loads(parcel_polygon_wkt)  # Convert WKT to Shapely Polygon
        merged_polygon = merged_polygon.union(parcel_polygon)
    else:
        # print(f"No parcels found for lat={lat}, lng={lng}")
        continue

polygon_strings = []
print(merged_polygon)

# Iterate through all polygons in the MultiPolygon
for polygon in merged_polygon.geoms:
    # Convert each polygon to a string format and append to the list
    polygon_string = ", ".join([f"{x} {y}" for x, y in polygon.exterior.coords])
    polygon_strings.append(f"(({polygon_string}))")

# Create a DataFrame with the polygon strings
polygon_df = pd.DataFrame({'Polygon': polygon_strings})


# Save the DataFrame to the CSV file
polygon_df.to_csv("multipolygon.csv", index=False)

def visualize_polygon(polygon):
    # Create a GeoDataFrame with the Polygon
    gdf = gpd.GeoDataFrame(geometry=[polygon], crs='EPSG:4326')

    # Plot the Polygon
    gdf.plot()
def find_new_polygon(multipolygon, selected_polygon):
    # Check if the input is a MultiPolygon
    if not isinstance(multipolygon, MultiPolygon):
        raise ValueError("Input must be a MultiPolygon.")

    # Check if the selected polygon is valid
    if not selected_polygon.is_valid:
        raise ValueError("Selected polygon is not valid.")

    # Get the boundary of the selected polygon
    selected_boundary = selected_polygon.boundary

    # Find the adjacent polygons that intersect with the selected boundary
    new_polygons = []
    for polygon in multipolygon.geoms:
        if polygon.intersects(selected_boundary):
            new_polygons.append(polygon)

    return new_polygons


selected_polygon_index = 0
selected_polygon = merged_polygon.geoms[selected_polygon_index]
print("Is selected polygon valid?", selected_polygon.is_valid)

# Find the adjacent polygons and verify their validity
adjacent_polygons = find_new_polygon(merged_polygon, selected_polygon)

# Visualize the selected polygon and each adjacent polygon separately
print("Visualizing the selected polygon (index 0):")
visualize_polygon(selected_polygon)

print("Visualizing adjacent polygons:")
for i, adjacent_polygon in enumerate(adjacent_polygons):
    print(f"Adjacent polygon {i + 1}:")
    print(adjacent_polygon)
    print("Is adjacent polygon valid?", adjacent_polygon.is_valid)
    visualize_polygon(adjacent_polygon)
# Plot the original polygon in blue
x, y = polygon.exterior.xy
plt.plot(x, y, color='blue')

# Check if merged_polygon is a MultiPolygon or a single Polygon
if isinstance(merged_polygon, MultiPolygon):
    for polygon in merged_polygon.geoms:
        x, y = polygon.exterior.xy
        plt.plot(x, y, color='red')
else:
    x, y = merged_polygon.exterior.xy
    plt.plot(x, y, color='red')

plt.show()

