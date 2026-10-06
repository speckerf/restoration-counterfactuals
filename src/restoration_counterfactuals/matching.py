"""Single-site matching, adapted from 01_find_candidate_donors.py."""

import ee
from loguru import logger

from .datasets import (
    load_alpha_earth,
    load_hansen_disturbance,
    load_human_modification,
    load_lulc,
)


def scale_by_sqrt_area(area_ha, minimum, maximum):
    fraction = (area_ha**0.5 - 1) / (1000**0.5 - 1)
    return float(min(max(minimum + fraction * (maximum - minimum), minimum), maximum))


def donor_search_area(
    geometry: ee.Geometry | ee.Feature,
    area_ha: int,
    settings: dict,
    other_projects: ee.FeatureCollection | None,
) -> ee.Geometry:
    """Exclude spillover and buffered local restoration boundaries (required input)."""
    outer = settings["INCLUSION_BUFFER"]
    inner = settings["SPILLOVER_EXCLUDE"]

    inclusion_buffer = scale_by_sqrt_area(area_ha, outer["MIN"], outer["MAX"])
    logger.debug(f"Inclusion buffer: {inclusion_buffer} m")
    region = geometry.buffer(inclusion_buffer, 10)

    exclusion_buffer = scale_by_sqrt_area(area_ha, inner["MIN"], inner["MAX"])
    logger.debug(f"Exclusion (spillover) buffer: {exclusion_buffer} m")
    region = region.difference(geometry.buffer(exclusion_buffer, 10), 10)

    # Buffer each polygon first, including projects just outside the search perimeter.
    buffered = other_projects.map(
        lambda f: f.buffer(settings["BUFFER_AROUND_OTHER_PROJECTS_M"], 10)
    ).filterBounds(region)

    return ee.Geometry(
        ee.Algorithms.If(
            buffered.size().gt(0), region.difference(buffered.geometry(10), 10), region
        )
    )


def matching_layers(intervention_year: int, bounds: ee.Geometry):
    if intervention_year < 2001:
        raise ValueError(
            "Pre-intervention GLAD land cover requires an intervention after 2000."
        )
    loss = load_hansen_disturbance(intervention_year=intervention_year)
    # The source matcher uses 9999 for no pre-intervention loss, including later loss.
    loss = loss.where(loss.eq(0), 9999).unmask(9999)
    hmi = load_human_modification(intervention_year - 1)

    # post 2017: use alpha earth by default
    mode = "alphaearth" if intervention_year >= 2018 else "landcover"
    lulc = (
        load_alpha_earth(intervention_year - 1, bounds)
        if mode == "alphaearth"
        else load_lulc(intervention_year - 1)
    )
    # otherwise use GLAD LULC simplified land cover map.
    if mode == "landcover":
        lulc = lulc.updateMask(lulc.neq(255))
    return {"disturbance": loss, "hmi": hmi, "lulc": lulc, "mode": mode}


def sample_treatment_pixels(site, layers, n_pixels, scale, seed=1111):
    image = (
        ee.Image.constant(1)
        .rename("valid_treatment")
        .addBands(ee.Image.cat([layers["disturbance"], layers["hmi"], layers["lulc"]]))
    )
    sampled = image.stratifiedSample(
        numPoints=n_pixels,
        classBand="valid_treatment",
        region=site.geometry(),
        scale=scale,
        seed=seed,
        dropNulls=True,
        tileScale=4,
        geometries=True,
    )
    return sampled.map(lambda f: f.set("treatment_id", f.id()))


def similarity_masks(pixel, layers, intervention_year, settings):
    """Return separate and combined masks relative to one treatment pixel."""
    hmi = (
        layers["hmi"]
        .subtract(pixel.getNumber("human_modification"))
        .abs()
        .lte(settings["HMI_MAX_DELTA"])
    )
    loss = layers["disturbance"]
    target = pixel.getNumber("lossyear")
    disturbed = target.gt(0).And(target.lt(intervention_year - 2000))
    disturbance = ee.Image(
        ee.Algorithms.If(
            disturbed,
            loss.subtract(target)
            .abs()
            .lte(settings["DISTURBANCE_TOLERANCE"])
            .And(loss.lt(intervention_year - 2000)),
            loss.eq(9999),
        )
    )
    if layers["mode"] == "alphaearth":
        bands = ee.List([f"A{i:02d}" for i in range(64)])
        vector = ee.Image.constant(bands.map(lambda b: pixel.getNumber(b))).rename(
            bands
        )
        distance = ee.Image.constant(1).subtract(
            layers["lulc"].multiply(vector).reduce(ee.Reducer.sum())
        )
        lulc = distance.lte(settings["AEF_MAX_DELTA"])
    else:
        lulc = layers["lulc"].eq(pixel.getNumber("lulc"))
    return {
        "human_modification": hmi,
        "disturbance": disturbance,
        "ecological_similarity": lulc,
        "eligible": hmi.And(disturbance).And(lulc).rename("eligible"),
    }


def select_donors(
    treatments: ee.FeatureCollection,
    region: ee.Geometry,
    layers: dict,
    intervention_year: int,
    settings: dict,
    seed=1234,
) -> ee.FeatureCollection:
    """Return unique donor locations; duplicate locations do not get extra ridge weight."""

    def sample(pixel):
        mask = similarity_masks(pixel, layers, intervention_year, settings)[
            "eligible"
        ].selfMask()
        points = mask.stratifiedSample(
            numPoints=settings["N_SELECTED_DONOR_PIXELS_PER_TREATMENT"],
            classBand="eligible",
            region=region,
            scale=settings["SAMPLE_SCALE"],
            seed=ee.Number.parse(pixel.get("treatment_id")).add(seed).toInt(),
            dropNulls=True,
            tileScale=4,
            geometries=True,
        )
        return points.map(lambda f: f.set("treatment_id", pixel.get("treatment_id")))

    donors = ee.FeatureCollection(treatments.map(sample)).flatten()

    def identify(f):
        coords = f.geometry().coordinates()
        key = (
            ee.Number(coords.get(0))
            .format("%.8f")
            .cat("_")
            .cat(ee.Number(coords.get(1)).format("%.8f"))
        )
        return f.set("donor_id", key)

    return donors.map(identify).distinct(["donor_id"])
