import requests
import json
from sgp4.api import Satrec, jday
from datetime import datetime, timedelta
import math

# -------------------------------
# Part 1: Download and Save TLE Data
# -------------------------------

def download_tle_data(url):
    try:
        print(f"Fetching data from {url}...")
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching the data: {e}")
        return None

# URL for Starlink TLE data from CelesTrak
tle_url = "https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle"
tle_text = download_tle_data(tle_url)

if tle_text is None:
    exit("Failed to download TLE data.")

# Optionally, save the raw TLE to a file (if needed for debugging)
tle_file_path = "starlink_tle.txt"
with open(tle_file_path, "w") as f:
    f.write(tle_text)
    print(f"TLE data saved to {tle_file_path}.")

# -------------------------------
# Part 2: Parse TLE Data into Satellite Entries
# -------------------------------

def read_tle_data(tle_text):
    satellites = []
    lines = [line.strip() for line in tle_text.splitlines() if line.strip()]
    for i in range(0, len(lines), 3):
        try:
            satellite_name = lines[i]
            line1 = lines[i + 1]
            line2 = lines[i + 2]
            if line1.startswith("1 ") and line2.startswith("2 "):
                satellites.append((satellite_name, line1, line2))
            else:
                print(f"Skipping malformed TLE structure: {satellite_name}")
        except IndexError:
            print(f"Skipping incomplete TLE entry starting at line {i + 1}")
    return satellites

satellites = read_tle_data(tle_text)

# -------------------------------
# Helper Functions for Orbit Calculation
# -------------------------------

def get_launch_group(line1):
    try:
        international_designator = line1[9:17]
        # Handle two-digit year: if > "57", assume 1900s, else 2000s.
        launch_year = "19" + international_designator[:2] if international_designator[:2] > "57" else "20" + international_designator[:2]
        launch_number = international_designator[2:5]
        return f"{launch_year}-Launch-{launch_number}"
    except Exception as e:
        print(f"Error extracting launch group: {line1}, Error: {e}")
        return "Unknown"

def is_tle_valid(line1):
    if len(line1) < 32:
        print(f"Invalid TLE format: {line1}")
        return False
    try:
        epoch_year = int(line1[18:20])
        epoch_day = float(line1[20:32])
        epoch_year += 2000 if epoch_year < 57 else 1900
        epoch_date = datetime(epoch_year, 1, 1) + timedelta(days=epoch_day - 1)
        days_diff = (datetime.utcnow() - epoch_date).days
        return days_diff < 30  # Consider TLE valid if epoch is within the last 30 days
    except (ValueError, IndexError) as e:
        print(f"Error parsing TLE line: {line1}, Error: {e}")
        return False

def eci_to_latlonalt(eci_position):
    EARTH_RADIUS = 6371.0  # in kilometers
    x, y, z = eci_position
    lon = math.degrees(math.atan2(y, x))
    lat = math.degrees(math.atan2(z, math.sqrt(x**2 + y**2)))
    alt = math.sqrt(x**2 + y**2 + z**2) - EARTH_RADIUS
    # Normalize longitude to [-180, 180]
    lon = (lon + 180) % 360 - 180
    # Return altitude in meters
    return lat, lon, alt * 1000

def is_valid_point(lat, lon):
    return -90 <= lat <= 90 and -180 <= lon <= 180

def fix_idl_crossings(trajectory):
    fixed_trajectory = []
    prev_lon = None
    for point in trajectory:
        lon, lat, alt = point
        if prev_lon is not None and abs(lon - prev_lon) > 180:
            lon += -360 if lon > prev_lon else 360
        fixed_trajectory.append([lon, lat, alt])
        prev_lon = lon
    # Re-normalize longitudes
    for point in fixed_trajectory:
        point[0] = (point[0] + 180) % 360 - 180
    return fixed_trajectory

def generate_trajectory_with_timestamps(satellite, start_time, duration_hours, time_step_minutes):
    trajectory = []
    timestamps = []
    current_time = start_time
    total_steps = int((duration_hours * 60) / time_step_minutes)
    
    for _ in range(total_steps):
        jd, fr = jday(current_time.year, current_time.month, current_time.day,
                      current_time.hour, current_time.minute, current_time.second)
        error_code, position, velocity = satellite.sgp4(jd, fr)
        
        if error_code == 0:
            lat, lon, alt = eci_to_latlonalt(position)
            if is_valid_point(lat, lon):
                trajectory.append([lon, lat, alt])
                timestamps.append(current_time.isoformat())
        else:
            raise RuntimeError(f"SGP4 error code {error_code} for satellite")
        
        current_time += timedelta(minutes=time_step_minutes)
    
    trajectory = fix_idl_crossings(trajectory)
    return trajectory, timestamps

# -------------------------------
# Part 3: Generate GeoJSON from TLE Data
# -------------------------------

output_geojson_file = "starlink_trajectories.geojson"
features = []
skipped_satellites = []

start_time = datetime.utcnow()
duration_hours = 24       # Orbit calculation is set for the next 24 hours.
time_step_minutes = 5     # Time step between trajectory points

for satellite_name, line1, line2 in satellites:
    if not is_tle_valid(line1):
        print(f"Skipping {satellite_name}: TLE is outdated or invalid.")
        skipped_satellites.append(satellite_name)
        continue

    try:
        satellite = Satrec.twoline2rv(line1, line2)
        trajectory, timestamps = generate_trajectory_with_timestamps(
            satellite, start_time, duration_hours, time_step_minutes
        )
        launch_group = get_launch_group(line1)
        features.append({
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": trajectory
            },
            "properties": {
                "name": satellite_name,
                "timestamps": timestamps,
                "launch_group": launch_group
            }
        })
    except Exception as e:
        print(f"Skipping {satellite_name} due to error: {e}")
        skipped_satellites.append(satellite_name)

geojson = {
    "type": "FeatureCollection",
    "features": features
}

with open(output_geojson_file, "w") as f:
    json.dump(geojson, f, indent=2)
    print(f"\nGeoJSON trajectories saved to {output_geojson_file}.")

if skipped_satellites:
    print("\nSkipped Satellites:")
    for skipped in skipped_satellites:
        print(f"- {skipped}")
