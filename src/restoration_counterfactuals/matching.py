"""Single-site matching, adapted from 01_find_candidate_donors.py."""
import ee
from .datasets import load_alpha_earth, load_hansen_disturbance, load_human_modification, load_lulc


def scale_by_sqrt_area(area_ha, minimum, maximum):
    fraction = (area_ha ** 0.5 - 1) / (1000 ** 0.5 - 1)
    return float(min(max(minimum + fraction * (maximum - minimum), minimum), maximum))


def donor_search_area(geometry, area_ha, other_projects, settings):
    """Exclude spillover and buffered local restoration boundaries (required input)."""
    outer = settings['INCLUSION_BUFFER']
    inner = settings['SPILLOVER_EXCLUDE']
    region = geometry.buffer(scale_by_sqrt_area(area_ha, outer['MIN'], outer['MAX']), 10)
    region = region.difference(geometry.buffer(scale_by_sqrt_area(area_ha, inner['MIN'], inner['MAX']), 10), 10)
    # Buffer each polygon first, including projects just outside the search perimeter.
    buffered = other_projects.map(lambda f: f.buffer(settings['BUFFER_AROUND_OTHER_PROJECTS_M'], 10)).filterBounds(region)
    return ee.Geometry(ee.Algorithms.If(buffered.size().gt(0), region.difference(buffered.geometry(10), 10), region))


def matching_layers(intervention_year, bounds):
    if intervention_year < 2001:
        raise ValueError('Pre-intervention GLAD land cover requires an intervention after 2000.')
    loss = load_hansen_disturbance(intervention_year=intervention_year)
    # The source matcher uses 9999 for no pre-intervention loss, including later loss.
    loss = loss.where(loss.eq(0), 9999).unmask(9999)
    hmi = load_human_modification(intervention_year - 1)
    mode = 'alphaearth' if intervention_year >= 2018 else 'landcover'
    ecology = (load_alpha_earth(intervention_year - 1, bounds) if mode == 'alphaearth'
               else load_lulc(intervention_year - 1))
    if mode == 'landcover':
        ecology = ecology.updateMask(ecology.neq(255))
    return dict(disturbance=loss, hmi=hmi, ecology=ecology, mode=mode)


def sample_treatment_pixels(site, layers, n_pixels, scale, seed=1234):
    image = ee.Image.constant(1).rename('valid_treatment').addBands(
        ee.Image.cat([layers['disturbance'], layers['hmi'], layers['ecology']]))
    sampled = image.stratifiedSample(numPoints=n_pixels, classBand='valid_treatment',
        region=site.geometry(), scale=scale, seed=seed, dropNulls=True, tileScale=4, geometries=True)
    return sampled.map(lambda f: f.set('treatment_id', f.id()))


def similarity_masks(pixel, layers, intervention_year, settings):
    """Return separate and combined masks relative to one treatment pixel."""
    hmi = layers['hmi'].subtract(pixel.getNumber('human_modification')).abs().lte(settings['HMI_MAX_DELTA'])
    loss = layers['disturbance']
    target = pixel.getNumber('lossyear')
    disturbed = target.gt(0).And(target.lt(intervention_year - 2000))
    disturbance = ee.Image(ee.Algorithms.If(disturbed,
        loss.subtract(target).abs().lte(settings['DISTURBANCE_TOLERANCE']).And(loss.lt(intervention_year - 2000)),
        loss.eq(9999)))
    if layers['mode'] == 'alphaearth':
        bands = ee.List([f'A{i:02d}' for i in range(64)])
        vector = ee.Image.constant(bands.map(lambda b: pixel.getNumber(b))).rename(bands)
        distance = ee.Image.constant(1).subtract(layers['ecology'].multiply(vector).reduce(ee.Reducer.sum()))
        ecology = distance.lte(settings['AEF_MAX_DELTA'])
    else:
        ecology = layers['ecology'].eq(pixel.getNumber('lulc'))
    return dict(human_modification=hmi, disturbance=disturbance, ecological_similarity=ecology,
                eligible=hmi.And(disturbance).And(ecology).rename('eligible'))


def select_donors(treatments, region, layers, intervention_year, settings, seed=1234):
    """Return unique donor locations; duplicate locations do not get extra ridge weight."""
    def sample(pixel):
        mask = similarity_masks(pixel, layers, intervention_year, settings)['eligible'].selfMask()
        points = mask.stratifiedSample(
            numPoints=settings['N_SELECTED_DONOR_PIXELS_PER_TREATMENT'], classBand='eligible',
            region=region, scale=settings['SAMPLE_SCALE'],
            seed=ee.Number.parse(pixel.get('treatment_id')).add(seed).toInt(),
            dropNulls=True, tileScale=4, geometries=True)
        return points.map(lambda f: f.set('treatment_id', pixel.get('treatment_id')))
    donors = ee.FeatureCollection(treatments.map(sample)).flatten()
    def identify(f):
        coords = f.geometry().coordinates()
        key = ee.Number(coords.get(0)).format('%.8f').cat('_').cat(ee.Number(coords.get(1)).format('%.8f'))
        return f.set('donor_id', key)
    return donors.map(identify).distinct(['donor_id'])
