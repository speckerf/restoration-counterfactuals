import ee
from loguru import logger


def load_alpha_earth(year: int, bounds: ee.Geometry) -> ee.Image:
    """Load annual AlphaEarth/Satellite Embedding image for a specific year."""
    assert 2017 <= year <= 2025, "Year must be between 2017 and 2025 inclusive."
    return (
        ee.ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")
        .filterDate(f"{year}-01-01", f"{year + 1}-01-01")
        .filterBounds(bounds)
        .mosaic()
    )


def load_hansen_disturbance(
    gfc_asset: str = "UMD/hansen/global_forest_change_2025_v1_13",
    intervention_year: int | None = None,
) -> ee.Image:
    """
    Load Hansen loss year and mask post intervention disturbances if intervention_year is provided.

    Hansen lossyear:
      0 = no loss
      1 = 2001
      ...
      25 = 2025


    """
    loss_year = ee.Image(gfc_asset).select("lossyear")
    loss_binary = ee.Image(gfc_asset).select("loss")

    if intervention_year is not None:
        # Keep only losses strictly before intervention year.
        # Example: intervention_year = 2011 -> keep lossyear <= 10.
        max_pre_lossyear = intervention_year - 2000 - 1

        pre_intervention_mask = loss_binary.eq(0).Or(
            loss_year.lte(ee.Image.constant(max_pre_lossyear))
        )

        loss_year = loss_year.updateMask(pre_intervention_mask)

    return loss_year.toInt()  # explicit cast from short to int


def load_human_modification(
    year: int,
    asset: str = "projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M",
) -> ee.Image:
    """Load Global Human Modification overall layer for a specific year."""
    year = min(year, 2020)
    assert 1990 <= year <= 2020, "Year must be between 1990 and 2020 inclusive."

    if year % 5 != 0:
        logger.trace(
            f"Interpolating Human Modification for year {year} since it is not a multiple of 5."
        )
        hmi_before = (
            ee.ImageCollection(asset)
            .filter(ee.Filter.eq("year", (year // 5) * 5))
            .first()
        )
        hmi_after = (
            ee.ImageCollection(asset)
            .filter(ee.Filter.eq("year", ((year // 5) + 1) * 5))
            .first()
        )
        hmi_img = hmi_before.add(
            hmi_after.subtract(hmi_before).multiply((year % 5) / 5)
        )

    else:
        logger.trace(f"Loading Human Modification for year {year}.")
        hmi_img = ee.ImageCollection(asset).filter(ee.Filter.eq("year", year)).first()

    return hmi_img.rename("human_modification")


def _remap_lulc_classes(lulc_img: ee.Image) -> ee.Image:

    terra_firma_true_desert = list(range(0, 2))
    terra_firma_semi_arid = list(range(2, 19))
    terra_firma_dense_short_vegetation = list(range(19, 25))
    terra_firma_treecover_beloweq10m = list(range(25, 33))
    terra_firma_treecover_above10m = list(range(33, 49))

    wetland_salt_pan = list(range(100, 102))
    wetland_sparse_vegetation = list(range(102, 119))
    wetland_dense_short_vegetation = list(range(119, 125))
    wetland_treecover_beloweq10m = list(range(125, 133))
    wetland_treecover_above10m = list(range(133, 149))

    open_surface_water = list(range(200, 208))
    snow_ice = [241]
    cropland = [244]
    built_up = [250]
    ocean = [254]
    no_data = [255]

    # Define target classes
    CLASS_MAP = {
        1: terra_firma_true_desert,
        2: terra_firma_semi_arid,
        3: terra_firma_dense_short_vegetation,
        4: terra_firma_treecover_beloweq10m,
        5: terra_firma_treecover_above10m,
        6: wetland_salt_pan,
        7: wetland_sparse_vegetation,
        8: wetland_dense_short_vegetation,
        9: wetland_treecover_beloweq10m,
        10: wetland_treecover_above10m,
        11: open_surface_water,
        12: snow_ice,
        13: cropland,
        14: built_up,
        15: ocean,
        255: no_data,
    }

    # Build remap lists
    from_values = []
    to_values = []

    for new_class, old_classes in CLASS_MAP.items():
        from_values.extend(old_classes)
        to_values.extend([new_class] * len(old_classes))

    # Remap
    lulc_remapped = lulc_img.remap(from_values, to_values, defaultValue=255).rename(
        "lulc"
    )

    return lulc_remapped


def load_lulc(
    year: int,
) -> ee.Image:
    """Load Global Land Analysis and Discovery (GLAD) Land Use/Land Cover dataset for the year of intervention. Load from 5 year image preintervention. Remap to main classes (e.g. forest, agriculture, urban, etc.)"""
    assert year >= 2000, "Intervention year must be greater than 2000."
    available_years = [2000, 2005, 2010, 2015, 2020]
    closest_year = max([y for y in available_years if y <= year])
    logger.debug(
        f"Loading LULC data for year {year}. Closest available year in GLAD dataset is {closest_year}."
    )
    lulc_img = ee.Image(f"projects/glad/GLCLU2020/v2/LCLUC_{closest_year}")
    lulc_remapped = _remap_lulc_classes(lulc_img)
    return lulc_remapped


def load_svh_imgc() -> tuple[ee.ImageCollection, ee.Projection]:
    path = "projects/global-pasture-watch/assets/gsvh-30m/v1/short-veg-height_m"
    imgc = ee.ImageCollection(path)

    # scale by 0.1:
    imgc = imgc.map(lambda img: img.multiply(0.1).copyProperties(img, img.propertyNames()))

    projection = imgc.first().projection()
    return imgc, projection


def load_agb_imgc() -> tuple[ee.ImageCollection, ee.Projection]:
    path = "projects/sat-io/open-datasets/CTREES-GLOBAL-AGB-100M"
    imgc = ee.ImageCollection(path).select("agb")

    # scale by 0.1 and mask out < 0:
    def rescale_agb(image: ee.Image) -> ee.Image:
        scaled = image.multiply(0.1)
        return scaled.updateMask(scaled.gte(0)).copyProperties(
            image, image.propertyNames()
        )

    imgc = imgc.map(lambda img: rescale_agb(img))

    projection = imgc.first().projection()
    return imgc, projection
