import requests
import json
from sgp4.api import Satrec, jday
from datetime import datetime, timedelta
import math

# -------------------------------
# Part 1: Download and Save TLE Data
# -------------------------------

def download_tle_data(url, timeout_seconds=30):
    try:
        print(f"Fetching data from {url}...")
        response = requests.get(url, timeout=timeout_seconds)
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

WGS84_A = 6378.137  # Semi-major axis in kilometers
WGS84_B = 6356.7523142  # Semi-minor axis in kilometers
WGS84_E2 = 1 - (WGS84_B ** 2) / (WGS84_A ** 2)
WGS84_EP2 = (WGS84_A ** 2) / (WGS84_B ** 2) - 1


def gstime(jdut1):
    tut1 = (jdut1 - 2451545.0) / 36525.0
    temp = (-6.2e-6 * tut1**3) + (0.093104 * tut1**2) + ((876600.0 * 3600 + 8640184.812866) * tut1) + 67310.54841
    temp = (temp % 86400.0) / 240.0  # Convert seconds to degrees
    gmst = math.radians(temp % 360.0)
    return gmst


def teme_to_ecef(teme_position_km, jd_ut1):
    gmst = gstime(jd_ut1)
    cos_gmst = math.cos(gmst)
    sin_gmst = math.sin(gmst)
    x_teme, y_teme, z_teme = teme_position_km
    x_ecef = cos_gmst * x_teme + sin_gmst * y_teme
    y_ecef = -sin_gmst * x_teme + cos_gmst * y_teme
    z_ecef = z_teme
    return x_ecef, y_ecef, z_ecef


def ecef_to_geodetic(x_km, y_km, z_km):
    p = math.sqrt(x_km**2 + y_km**2)
    if p == 0.0:
        lon = 0.0
    else:
        lon = math.atan2(y_km, x_km)

    theta = math.atan2(z_km * WGS84_A, p * WGS84_B)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)

    lat = math.atan2(
        z_km + WGS84_EP2 * WGS84_B * sin_theta**3,
        p - WGS84_E2 * WGS84_A * cos_theta**3
    )

    sin_lat = math.sin(lat)
    N = WGS84_A / math.sqrt(1 - WGS84_E2 * sin_lat**2)
    alt_km = (p / math.cos(lat)) - N

    lon_deg = (math.degrees(lon) + 180.0) % 360.0 - 180.0
    lat_deg = math.degrees(lat)
    alt_m = alt_km * 1000.0
    return lat_deg, lon_deg, alt_m


def teme_to_latlonalt(position_km, jd, fr):
    jd_ut1 = jd + fr
    x_ecef, y_ecef, z_ecef = teme_to_ecef(position_km, jd_ut1)
    return ecef_to_geodetic(x_ecef, y_ecef, z_ecef)

def is_valid_point(lat, lon):
    return -90 <= lat <= 90 and -180 <= lon <= 180

def fix_idl_crossings(trajectory):
    if not trajectory:
        return []

    fixed = []
    offset = 0.0
    prev_lon = None

    for lon, lat, alt in trajectory:
        if prev_lon is not None:
            diff = (lon + offset) - prev_lon
            if diff > 180.0:
                offset -= 360.0
            elif diff < -180.0:
                offset += 360.0
        adjusted_lon = lon + offset
        fixed.append([adjusted_lon, lat, alt])
        prev_lon = adjusted_lon

    return fixed

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
            lat, lon, alt = teme_to_latlonalt(position, jd, fr)
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
