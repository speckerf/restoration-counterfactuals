"""Annual boundary means and donor-pixel outcomes on each dataset's native grid."""

import re

import ee
import pandas as pd


def annual_bands(collection):
    """Use source annual image IDs, as in the existing toBands extraction."""
    image = collection.toBands()
    names = image.bandNames().getInfo()
    years = []
    for name in names:
        matches = re.findall(r"(?:19|20)\d{2}", name)
        if len(matches) != 1:
            raise ValueError(f"Cannot identify a unique annual year in {name!r}.")
        years.append(int(matches[0]))
    if not years or len(set(years)) != len(years):
        raise ValueError("Expected exactly one outcome band per year.")
    return image.rename([f"y{year}" for year in years]), years


def _long_table(features, id_column, years, variable):
    rows = []
    for feature in features["features"]:
        properties = feature["properties"]
        for year in years:
            rows.append(
                {
                    id_column: properties[id_column],
                    "year": year,
                    "variable": variable,
                    "value": properties.get(f"y{year}"),
                }
            )
    return pd.DataFrame(rows, columns=[id_column, "year", "variable", "value"])


def extract_project_means(site, collection, projection, variable):
    """Area-weighted annual raster mean directly inside the project polygon."""
    image, years = annual_bands(collection)
    grid = projection.getInfo()
    values = image.reduceRegions(
        collection=ee.FeatureCollection([site]),
        reducer=ee.Reducer.mean(),
        crs=grid["crs"],
        crsTransform=grid["transform"],
    ).getInfo()
    return _long_table(values, "site_id", years, variable)


def extract_donor_values(donors, collection, projection, variable):
    """Retain masked annual values as missing, rather than dropping whole donors."""
    image, years = annual_bands(collection)
    grid = projection.getInfo()
    values = image.reduceRegions(
        collection=donors,
        reducer=ee.Reducer.first(),
        crs=grid["crs"],
        crsTransform=grid["transform"],
    ).getInfo()
    return _long_table(values, "donor_id", years, variable)
