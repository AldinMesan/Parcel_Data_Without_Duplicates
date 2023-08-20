from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import requests
import pandas as pd
import geopandas as gpd
import time
import random
from shapely.wkt import loads


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


def get_response_with_delay(lat, lng):
    response = get_response(lat, lng)
    time.sleep(random.uniform(0.5, 1))
    return response


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
    buffer_polygon = polygon.buffer(distance)

    buffer_points = []

    # Calculate points along the boundary of the buffer polygon at equal intervals
    for i in range(0, num_points):
        point_on_boundary = buffer_polygon.boundary.interpolate(i / num_points, normalized=True)
        buffer_points.append((point_on_boundary.x, point_on_boundary.y))

    return buffer_points


# Polygon vertices
polygon_vertices = [
    (-101.79553987431, 35.3090191916404),
    (-101.795486993006, 35.2945159538594),
    (-101.782048905939, 35.2944826309803),
    (-101.777676876069, 35.2944716474824),
    (-101.777576416348, 35.2944713938042),
    (-101.777636595966, 35.3090300262766),
    (-101.795381030267, 35.3090284543452),
    (-101.795539840892, 35.3090284222647),
    (-101.79553987431, 35.3090191916404)
]

# Create a shapely Polygon object
polygon = Polygon(polygon_vertices)

buffer_distance = 0.0001
# Find the four corners of the polygon
corners = find_polygon_corners(polygon)
# Find points on corners of the buffer polygon
buffer_points = find_points_on_corners(polygon, buffer_distance, num_points=12)
print("Buffer Points:", buffer_points)

df_buffer_points = pd.DataFrame(buffer_points, columns=['lng', 'lat'])

df_buffer_points.to_csv('first_buffer_points.csv', index=False)

df_coordinates = pd.read_csv('first_buffer_points.csv')

# Initialize lists to store the scraped data
# Initialize lists to store the scraped data for the first buffer zone
first_buffer_parcel_responses = []
for index, row in df_buffer_points.iterrows():
    lat = row['lat']
    lng = row['lng']
    response = get_response_with_delay(lat, lng)
    first_buffer_parcel_responses.append(response)

# Add the parcel data response to the df_buffer_points DataFrame
df_buffer_points['parcel_data'] = first_buffer_parcel_responses

# Save the DataFrame with first buffer points and parcel data to a CSV file
df_buffer_points.to_csv('first_buffer_points_parcel_data.csv', index=False)

second_buffer_distance = 0.0002

# Find points on corners of the second buffer polygon
second_buffer_points = find_points_on_corners(polygon, second_buffer_distance, num_points=12)
print("Second Buffer Points:", second_buffer_points)

# Create a DataFrame for the second buffer points
df_second_buffer_points = pd.DataFrame(second_buffer_points, columns=['lng', 'lat'])

# Initialize lists to store the scraped data for the second buffer zone
second_buffer_parcel_responses = []
for index, row in df_second_buffer_points.iterrows():
    lat = row['lat']
    lng = row['lng']
    response = get_response_with_delay(lat, lng)
    second_buffer_parcel_responses.append(response)

# Add the parcel data response to the df_second_buffer_points DataFrame
df_second_buffer_points['parcel_data'] = second_buffer_parcel_responses

# Save the DataFrame with second buffer points and parcel data to a CSV file
df_second_buffer_points.to_csv('second_buffer_points_parcel_data.csv', index=False)

# Extract x and y coordinates for scatter plot if points are available
if buffer_points:
    buffer_x, buffer_y = zip(*buffer_points)
else:
    buffer_x, buffer_y = [], []

# Calculate the midpoints between consecutive corners
midpoints = find_midpoints(corners)
print("Midpoints:", midpoints)
# Extract x and y coordinates for midpoints scatter plot if points are available
if midpoints:
    mid_x, mid_y = zip(*midpoints)
else:
    mid_x, mid_y = [], []

# Extract x and y coordinates for the second buffer zone points if available
if second_buffer_points:
    second_buffer_x, second_buffer_y = zip(*second_buffer_points)
else:
    second_buffer_x, second_buffer_y = [], []
    # Assuming you have calculated and obtained the 'second_buffer_points'

    import geopandas as gpd
    from shapely.geometry import Point

    # Create a DataFrame to store the second buffer points
    df_second_buffer = pd.DataFrame({'x': second_buffer_x, 'y': second_buffer_y})

    # Convert the DataFrame to a GeoDataFrame with Point geometries
    geometry = [Point(xy) for xy in zip(df_second_buffer['x'], df_second_buffer['y'])]
    gdf_second_buffer = gpd.GeoDataFrame(df_second_buffer, geometry=geometry, crs='EPSG:4326')

    # Perform the spatial join using 'within'
    points_within_first_buffer = gdf_second_buffer[gdf_second_buffer.geometry.within(polygon)]

    # Get the parcel data for points within the first buffer zone
    second_buffer_parcel_responses = []
    for index, row in points_within_first_buffer.iterrows():
        lat = row.geometry.y
        lng = row.geometry.x
        response = get_response_with_delay(lat, lng)
        second_buffer_parcel_responses.append(response)

# Visualize the polygons and points for both buffer zones
x, y = polygon.exterior.xy
x_buf, y_buf = polygon.buffer(buffer_distance).exterior.xy

# Create a new subplot for the visualizations
plt.figure(figsize=(12, 6))

# First subplot: Original and First Buffer Zone Visualization
plt.subplot(1, 2, 1)
plt.plot(x, y, label='Original Polygon')
plt.plot(x_buf, y_buf, label='Buffer Polygon')
plt.scatter(buffer_x, buffer_y, color='red', label='Points on Buffer Corners')
plt.scatter(mid_x, mid_y, color='green', label='Midpoints')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Original and Buffer Polygon Visualization')
plt.legend()
plt.grid(True)

# Second subplot: Second Buffer Zone Visualization
plt.subplot(1, 2, 2)
# Visualize the second buffer polygon
x_second_buf, y_second_buf = polygon.buffer(second_buffer_distance).exterior.xy
plt.plot(x, y, label='Original Polygon')
plt.plot(x_second_buf, y_second_buf, label='Second Buffer Polygon',
         linestyle='dashed')  # Dashed line for differentiation
plt.scatter(buffer_x, buffer_y, color='red', label='Points on Buffer Corners (1st Buffer)')
plt.scatter(mid_x, mid_y, color='green', label='Midpoints (1st Buffer)')
plt.scatter(second_buffer_x, second_buffer_y, color='blue', label='Points on Buffer Corners (2nd Buffer)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Second Buffer Polygon Visualization')
plt.legend()
plt.grid(True)

plt.tight_layout()  # To adjust spacing between subplots


# plt.show()
def extract_lat_lon(wkt):
    polygon = loads(wkt)
    # Assuming the WKT is in the form 'POLYGON((lng1 lat1, lng2 lat2, ..., lngN latN))'
    coords = polygon.exterior.coords
    return coords


def second_within_first(gdf_first_buffer, gdf_second_buffer, matching_gdff, non_matching_gdff):
    """
Find geometries from the second buffer that are within the first buffer and separate them into matching and
    non-matching GeoDataFrames.

    :param gdf_first_buffer: str
        File path to the first buffer GeoDataFrame in CSV format.
    :param gdf_second_buffer: str
        File path to the second buffer GeoDataFrame in CSV format.
    :param matching_gdff: str
        File path where the matching geometries will be saved in CSV format.
    :param non_matching_gdff: str
        File path where the non-matching geometries will be saved in CSV format.
    :return: tuple
        A tuple containing two GeoDataFrames. The first one contains the geometries from the second buffer that
        are within the first buffer, and the second one contains the non-matching geometries.
    """

    gdf_first_buffer = gpd.read_file(gdf_first_buffer)
    gdf_second_buffer = gpd.read_file(gdf_second_buffer)
    # gdf_first_buffer = gpd.read_file('combined_data_first_buffer_no_duplicates.csv')
    # gdf_second_buffer = gpd.read_file('combined_data_second_buffer_no_duplicates.csv')

    # Step 2: Convert 'geom_as_wkt' to valid geometries using Shapely
    gdf_first_buffer['geometry'] = gdf_first_buffer['parcel_data'].apply(
        lambda x: loads(eval(x)['parcels'][0]['parcel_data']['geom_as_wkt']))
    gdf_second_buffer['geometry'] = gdf_second_buffer['parcel_data'].apply(
        lambda x: loads(eval(x)['parcels'][0]['parcel_data']['geom_as_wkt']))

    # Step 3: Drop the 'parcel_data' column, as it is no longer needed
    gdf_first_buffer.drop(columns=['parcel_data'], inplace=True)
    gdf_second_buffer.drop(columns=['parcel_data'], inplace=True)

    # Step 4: Drop rows with null geometries
    gdf_first_buffer = gdf_first_buffer[gdf_first_buffer.geometry.notnull()]
    gdf_second_buffer = gdf_second_buffer[gdf_second_buffer.geometry.notnull()]
    print(gdf_first_buffer)
    print(gdf_second_buffer)

    matching_parcel_indices = []
    for idx, first_geometry in gdf_first_buffer.iterrows():
        is_within = gdf_second_buffer.geometry.within(first_geometry['geometry'])
        if is_within.any():
            matching_indices = gdf_second_buffer.index[is_within].tolist()
            matching_parcel_indices.extend([(idx, second_idx) for second_idx in matching_indices])

            matching_gdf = gdf_first_buffer.iloc[[idx for idx, _ in matching_parcel_indices]]
            non_matching_gdf = gdf_first_buffer[
                ~gdf_first_buffer.index.isin([idx for idx, _ in matching_parcel_indices])]

            matching_gdf.to_csv(matching_gdff, index=False)
            non_matching_gdf.to_csv(non_matching_gdff, index=False)

            matching_parcel_gdf = gpd.GeoDataFrame(
                gdf_second_buffer.loc[[index for _, index in matching_parcel_indices]],
                crs=gdf_second_buffer.crs
            )
            
            
            fig, ax = plt.subplots(figsize=(10, 10))
            gdf_first_buffer.plot(ax=ax, color='blue', label='First Buffer Parcels')
            gdf_second_buffer.plot(ax=ax, color='green', label='Second Buffer Parcels')
            matching_parcel_gdf.plot(ax=ax, color='red', label='Matching Parcels')

            # Add labels to the points for better identification
            for x, y, label in zip(gdf_first_buffer.geometry.centroid.x, gdf_first_buffer.geometry.centroid.y,
                                   gdf_first_buffer.index):
                ax.text(x, y, label, color='blue', fontsize=8, ha='center', va='center')

            for x, y, label in zip(gdf_second_buffer.geometry.centroid.x, gdf_second_buffer.geometry.centroid.y,
                                   gdf_second_buffer.index):
                ax.text(x, y, label, color='green', fontsize=8, ha='center', va='center')

            for x, y, label in zip(matching_parcel_gdf.geometry.centroid.x, matching_parcel_gdf.geometry.centroid.y,
                                   matching_parcel_gdf.index):
                ax.text(x, y, label, color='red', fontsize=8, ha='center', va='center')

            # Create custom legend labels using proxy artists
            blue_patch = plt.Rectangle((0, 0), 1, 1, color='blue', label='First Buffer Parcels')
            green_patch = plt.Rectangle((0, 0), 1, 1, color='green', label='Second Buffer Parcels')
            red_patch = plt.Rectangle((0, 0), 1, 1, color='red', label='Matching Parcels')
            ax.legend(handles=[blue_patch, green_patch, red_patch], loc='upper right')

            # Set plot title
            ax.set_title('Matching Parcels Visualization')
            plt.show()


def extract_parcel_id(input_csv, output_csv):
    df = pd.read_csv(input_csv, delimiter=",")

    all_data = []

    for _, row in df.iterrows():
        # Replace single quotes with double quotes in the JSON-like string
        parcel_data_str = row['parcel_data'].replace("'", "\"")

        # Find the start and end index of the 'parcel_id' field
        start_index = parcel_data_str.find('"parcel_id":')
        if start_index != -1:
            end_index = parcel_data_str.find(',', start_index)
            if end_index == -1:
                end_index = parcel_data_str.find('}', start_index)
            # Extract the 'parcel_id' value from the string
            parcel_id = parcel_data_str[start_index + len('"parcel_id":'):end_index].strip(' "')
            # Add the parcel_id to the list
            all_data.append({'parcel_id': parcel_id})
        else:
            print(f"No parcel_id found in: {parcel_data_str}")

    # Create a new DataFrame from the list of parcel IDs
    merged_df = pd.DataFrame(all_data)

    # Save the DataFrame to the output CSV file
    merged_df.to_csv(output_csv, index=False)


def combine_data(input_csv, parcel_id_csv, output_csv):
    """
    Combine 'lat', 'lng', 'parcel_data', and 'parcel_id' from two CSV files and save them to a new CSV file.

    :param input_csv: Path to the input CSV file containing 'lat', 'lng', and 'parcel_data' columns.
    :param parcel_id_csv: Path to the CSV file containing the previously extracted 'parcel_id'.
    :param output_csv: Path to the output CSV file where the combined data will be saved.
    :return: None
    """

    df_input = pd.read_csv(input_csv)
    df_parcel_id = pd.read_csv(parcel_id_csv)
    df_combined = pd.concat([df_input, df_parcel_id], axis=1)
    df_combined.to_csv(output_csv, index=False)


def delete_duplicates(input_csv, output_csv):
    """
    Delete duplicates from the input CSV file and save the result to a new CSV file.

    :param input_csv: Path to the input CSV file.
    :param output_csv: Path to the output CSV file where the data without duplicates will be saved.
    :return: None
    """
    df = pd.read_csv(input_csv)
    df_no_duplicates = df.drop_duplicates(subset='parcel_id', keep='first')
    df_no_duplicates.to_csv(output_csv, index=False)


# Usage example:

if __name__ == "__main__":

    extract_parcel_id('first_buffer_points_parcel_data.csv', 'first_buffer_points_parcel_ids.csv')
    extract_parcel_id('second_buffer_points_parcel_data.csv', 'second_buffer_points_parcel_ids.csv')

    combine_data('first_buffer_points_parcel_data.csv', 'first_buffer_points_parcel_ids.csv',
                 'combined_data_first_buffer.csv')
    combine_data('second_buffer_points_parcel_data.csv', 'second_buffer_points_parcel_ids.csv',
                 'combined_data_second_buffer.csv')

    delete_duplicates('combined_data_first_buffer.csv', 'combined_data_first_buffer_no_duplicates.csv')
    delete_duplicates('combined_data_second_buffer.csv', 'combined_data_second_buffer_no_duplicates.csv')

    second_within_first('combined_data_first_buffer_no_duplicates.csv', 'combined_data_second_buffer_no_duplicates.csv',
                        'matching_parcels.csv', 'non_matching_parcels.csv')

