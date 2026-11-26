#!/usr/bin/env python3
"""
geography.py

Turn what a BioSample says about where a mosquito came from into a province.

Used by: tables/make_china_province_table.py, hosts/lookup_host_species.py,
tables/summarise_bioproject.py
"""

import json
import re
import time
import urllib.request

# Set the coordinate pattern, the Nominatim settings + the China province lookups
# Nominatim allows one request a second, zoom 12 is city level, and the four municipalities
# name a district rather than a province
LAT_LON = re.compile(r"\s*([\d.]+)\s*([NS])[\s,]+([\d.]+)\s*([EW])")

NOMINATIM = "https://nominatim.openstreetmap.org/reverse"
NOMINATIM_PAUSE_SECONDS = 1.1
NOMINATIM_RETRIES = 3
NOMINATIM_TIMEOUT_SECONDS = 60
USER_AGENT = "biof3001-project7-summary/1.0"

NOMINATIM_ZOOM = 12

MUNICIPALITY_DISTRICTS = {
    "Dongcheng": "Beijing", "Xicheng": "Beijing", "Chaoyang": "Beijing",
    "Haidian": "Beijing", "Fengtai": "Beijing",
    "Pudong": "Shanghai", "Minhang": "Shanghai", "Jiading": "Shanghai",
    "Songjiang": "Shanghai",
}


CHINESE_PROVINCES = {
    "北京": "Beijing", "上海": "Shanghai", "云南": "Yunnan", "内蒙古": "Inner Mongolia",
    "四川": "Sichuan", "宁夏": "Ningxia", "安徽": "Anhui", "山东": "Shandong",
    "山西": "Shanxi", "广东": "Guangdong", "广西": "Guangxi", "新疆": "Xinjiang",
    "江西": "Jiangxi", "河南": "Henan", "浙江": "Zhejiang", "海南": "Hainan",
    "湖北": "Hubei", "湖南": "Hunan", "甘肃": "Gansu", "西藏": "Xizang",
    "贵州": "Guizhou", "陕西": "Shaanxi", "黑龙江": "Heilongjiang",
    "黑龙江省": "Heilongjiang",
}

PROVINCE_SUFFIX = re.compile(r"\s+(Province|City|Municipality|Autonomous Region|"
                             r"Zhuang Autonomous Region|Hui Autonomous Region|"
                             r"Uygur Autonomous Region|District)$")


def split_geo_loc_name(value):
    """Split "China: Yunnan, Mengla" into ("Yunnan", "Mengla"), or ("", "") if country only."""
    if not value or ":" not in str(value):
        return "", ""

    location = str(value).split(":", 1)[1].strip()
    if "," not in location:
        return location.strip(), ""

    province, site = location.split(",", 1)
    return province.strip(), site.strip()


def parse_lat_lon(value):
    """Parse "31.61 N 104.82 E" into (31.61, 104.82), or None if it isn't coordinates."""
    match = LAT_LON.match(str(value))
    if not match:
        return None

    latitude = float(match.group(1)) * (1 if match.group(2) == "N" else -1)
    longitude = float(match.group(3)) * (1 if match.group(4) == "E" else -1)
    return latitude, longitude


def reverse_geocode(latitude, longitude):
    """Return the province, city and district for one coordinate, or {} if the lookup failed."""
    url = (f"{NOMINATIM}?format=jsonv2&zoom={NOMINATIM_ZOOM}&accept-language=en"
           f"&lat={latitude}&lon={longitude}")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    # Step 1: ask Nominatim, retrying a few times before giving up on the coordinate
    for attempt in range(NOMINATIM_RETRIES):
        try:
            with urllib.request.urlopen(request, timeout=NOMINATIM_TIMEOUT_SECONDS) as response:
                address = json.load(response).get("address", {})
            break
        except Exception as error:
            print(f"  geocoder request failed ({error}), retrying...")
            time.sleep(5)
    else:
        return {}

    # Step 2: read the names out, falling back to the city as a municipality has no "state"
    province = address.get("state") or address.get("province") or address.get("city", "")
    city = address.get("city") or address.get("state_district") or address.get("town", "")
    district = (address.get("city_district") or address.get("county")
                or address.get("district") or address.get("suburb", ""))

    return {"province": province, "city": city, "district": district}


def normalise_province(province):
    """Return one spelling for a province, handling Chinese names, suffixes + districts."""
    name = (province or "").strip()
    if not name:
        return ""

    for chinese, english in CHINESE_PROVINCES.items():
        if chinese in name:
            return english

    name = PROVINCE_SUFFIX.sub("", name).strip()
    return MUNICIPALITY_DISTRICTS.get(name, name)
