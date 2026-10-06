"""Return maps and figures; callers decide how to display and save them."""

import ee
import geemap
import glasbey
import ipywidgets as widgets
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

# def matching_map(
#     site: ee.Feature,
#     region: ee.Geometry,
#     treatments: ee.FeatureCollection,
#     donors: ee.FeatureCollection,
#     masks: dict | None = None,
#     exclusions: ee.FeatureCollection | None = None,
# ) -> geemap.Map:

#     treatment_ids = treatments.aggregate_array("treatment_id").distinct().sort()
#     n = treatment_ids.size().getInfo()

#     palette = ee.List(
#         [color.lstrip("#") for color in glasbey.create_palette(palette_size=n)]
#     )

#     def style_treatment(feature):
#         index = treatment_ids.indexOf(feature.get("treatment_id"))
#         color = palette.get(ee.Number(index).mod(palette.size()))
#         return feature.set(
#             "style",
#             {
#                 "color": color,
#                 "fillColor": color,
#                 "pointSize": 6,
#                 "width": 1,
#             },
#         )

#     m = geemap.Map(basemap="HYBRID")
#     m.centerObject(site, 12)
#     m.addLayer(region, {"color": "888888"}, "Donor search area")

#     if exclusions is not None:
#         m.addLayer(exclusions, {"color": "ff8800"}, "Other project exclusions")

#     m.addLayer(site, {"color": "ff0000"}, "Project boundary")
#     m.addLayer(
#         treatments.map(style_treatment).style(styleProperty="style"),
#         {},
#         "Treatment samples",
#     )
#     m.addLayer(
#         donors.map(style_treatment).style(styleProperty="style"),
#         {},
#         "Selected donors",
#     )

#     if masks:
#         selector = widgets.Dropdown(
#             options=[(str(tid), tid) for tid in masks],
#             description="Treatment:",
#         )
#         output = widgets.Output()
#         mask_layers = []

#         @output.capture(clear_output=True)
#         def show_masks(treatment_id):
#             for layer in mask_layers:
#                 m.remove_layer(layer)
#             mask_layers.clear()

#             for i, (name, mask) in enumerate(masks[treatment_id].items()):
#                 layer_name = f"{name} — treatment {treatment_id}"
#                 m.addLayer(
#                     mask.selfMask().clip(region),
#                     {"palette": ["66bb66"]},
#                     layer_name,
#                     i == 0,  # Show the first mask; others remain toggleable.
#                     0.6,
#                 )
#                 mask_layers.append(m.find_layer(layer_name))

#         def on_change(change):
#             show_masks(change["new"])

#         selector.observe(on_change, names="value")
#         m.add_widget(
#             widgets.VBox([selector, output]),
#             position="bottomleft",
#         )
#         show_masks(selector.value)

#     return m


def matching_map(
    site: ee.Feature,
    region: ee.Geometry,
    treatments: ee.FeatureCollection,
    donors: ee.FeatureCollection,
    masks: dict | None = None,
    exclusions: ee.FeatureCollection | None = None,
) -> geemap.Map:

    treatment_ids = treatments.aggregate_array("treatment_id").distinct().sort()
    ids = treatment_ids.getInfo()

    colors = (
        [color.lstrip("#") for color in glasbey.create_palette(palette_size=len(ids))]
        if ids
        else []
    )

    palette = ee.List(colors)
    treatment_colors = dict(zip(ids, colors))

    def style_treatment(feature):
        index = treatment_ids.indexOf(feature.get("treatment_id"))
        color = palette.get(index)
        return feature.set(
            "style",
            {
                "color": color,
                "fillColor": color,
                "pointSize": 6,
                "width": 1,
            },
        )

    m = geemap.Map(basemap="HYBRID")
    m.centerObject(site, 12)
    m.addLayer(region, {"color": "888888"}, "Donor search area", opacity=0.6)

    if exclusions is not None:
        m.addLayer(
            exclusions, {"color": "ff8800"}, "Other project exclusions", opacity=0.6
        )

    m.addLayer(site, {"color": "ff0000"}, "Project boundary", opacity=0.6)

    selector = widgets.Dropdown(
        options=[("All treatments", None)] + [(str(tid), tid) for tid in ids],
        value=None,
        description="Treatment:",
        layout=widgets.Layout(width="250px"),
    )
    swatch = widgets.HTML(layout=widgets.Layout(width="24px"))
    output = widgets.Output()

    layer_cache = {}
    active_key = None

    def add_cached_layer(image, vis, name, shown=True, opacity=1):
        m.addLayer(image, vis, name, shown, opacity)
        return m.find_layer(name)

    def build_layers(treatment_id):
        """Build once per selection; subsequent visits reuse the layers."""
        layers = []

        if treatment_id is None:
            selected_treatments = treatments
            selected_donors = donors
            suffix = "all"
        else:
            selection = ee.Filter.eq("treatment_id", treatment_id)
            selected_treatments = treatments.filter(selection)
            selected_donors = donors.filter(selection)
            suffix = str(treatment_id)

            for i, (name, mask) in enumerate(
                (masks or {}).get(treatment_id, {}).items()
            ):
                layer = add_cached_layer(
                    mask.selfMask().clip(region),
                    {"palette": ["66bb66"]},
                    f"{name} — treatment {suffix}",
                    shown=i == 0 and "human" not in name.lower(),
                    opacity=0.6,
                )
                layers.append(layer)

        for collection, name in [
            (selected_treatments, "Treatment samples"),
            (selected_donors, "Selected donors"),
        ]:
            layer = add_cached_layer(
                collection.map(style_treatment).style(styleProperty="style"),
                {},
                f"{name} — {suffix}",
            )
            layers.append(layer)

        return layers

    @output.capture(clear_output=True)
    def show_selection(treatment_id):
        nonlocal active_key

        # Remember visibility and remove the previous selection's layers.
        if active_key in layer_cache:
            entry = layer_cache[active_key]
            entry["visibility"] = [layer.visible for layer in entry["layers"]]

            for layer in entry["layers"]:
                m.remove_layer(layer)

            entry["layers"] = []

        # Restore previously selected visibility, if available.
        previous_visibility = layer_cache.get(treatment_id, {}).get("visibility")

        layers = build_layers(treatment_id)

        if previous_visibility is not None:
            for layer, visible in zip(layers, previous_visibility):
                layer.visible = visible

        layer_cache[treatment_id] = {
            "layers": layers,
            "visibility": [layer.visible for layer in layers],
        }
        active_key = treatment_id

        if treatment_id is None:
            swatch.value = ""
        else:
            color = treatment_colors[treatment_id]
            swatch.value = (
                '<span style="display:inline-block;'
                "width:16px;height:16px;border-radius:50%;"
                f'background:#{color};border:1px solid #555;"></span>'
            )

    def on_change(change):
        show_selection(change["new"])

    selector.observe(on_change, names="value")

    m.add_widget(
        widgets.VBox(
            [
                widgets.HBox(
                    [swatch, selector],
                    layout=widgets.Layout(align_items="center"),
                ),
                output,
            ]
        ),
        position="bottomleft",
    )

    show_selection(selector.value)
    # show_selection(selector.value)

    return m


def impact_figure_with_donors(
    trajectories,
    donor_timeseries,
    weights,
    intervention_year,
    variable,
    zero_weight_tol=1e-10,
):
    """Plot the restoration site, synthetic control, donors, and estimated effect."""

    units = {"svh": "m", "agb": "Mg/ha"}
    unit = units.get(variable, "")

    donors = donor_timeseries.copy()
    if "variable" in donors.columns:
        donors = donors.loc[donors["variable"] == variable]

    donors = donors.merge(
        weights[["donor_id", "weight"]],
        on="donor_id",
        how="inner",
        validate="many_to_one",
    )

    zero_donors = donors.loc[donors["weight"] <= zero_weight_tol]
    weighted_donors = donors.loc[donors["weight"] > zero_weight_tol]

    # Light-to-dark blue; blue is reserved exclusively for donor trajectories.
    donor_cmap = LinearSegmentedColormap.from_list(
        "donor_blues",
        ["#C6DBEF", "#6BAED6", "#2171B5"],
    )

    max_weight = weighted_donors["weight"].max() if not weighted_donors.empty else 1.0
    weight_norm = Normalize(vmin=0, vmax=max_weight)

    colors = {
        "site": "#D55E00",  # vermilion
        "synthetic": "#222222",  # charcoal
        "zero_donor": "#B8C0CC",  # neutral gray-blue
        "effect": "#009E73",  # green
        "intervention": "#555555",
        "post_period": "#F5F5F5",
        "grid": "#E5E7EB",
    }

    fig, (ax, effect_ax) = plt.subplots(
        2,
        1,
        figsize=(11, 7.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]},
    )
    fig.subplots_adjust(
        top=0.84,
        bottom=0.10,
        left=0.10,
        right=0.97,
        hspace=0.08,
    )

    # Draw zero-weight donors first so they form quiet background context.
    for _, donor in zero_donors.groupby("donor_id"):
        donor = donor.sort_values("year")
        ax.plot(
            donor["year"],
            donor["value"],
            color=colors["zero_donor"],
            linewidth=0.65,
            alpha=0.12,
            zorder=1,
        )

    # Draw low-weight donors first and high-weight donors last.
    weighted_ids = (
        weighted_donors[["donor_id", "weight"]].drop_duplicates().sort_values("weight")
    )

    for row in weighted_ids.itertuples(index=False):
        donor = weighted_donors.loc[
            weighted_donors["donor_id"] == row.donor_id
        ].sort_values("year")

        relative_weight = row.weight / max_weight

        ax.plot(
            donor["year"],
            donor["value"],
            color=donor_cmap(weight_norm(row.weight)),
            # linewidth=0.8 + 2.2 * np.sqrt(relative_weight),
            linewidth=1,
            alpha=0.35 + 0.50 * np.sqrt(relative_weight),
            zorder=2,
        )

    # A subtle white outline separates the two focal trajectories from donors.
    foreground_effect = [
        pe.Stroke(linewidth=4.8, foreground="white"),
        pe.Normal(),
    ]

    (site_line,) = ax.plot(
        trajectories["year"],
        trajectories["observed"],
        color=colors["site"],
        marker="o",
        markersize=6,
        markerfacecolor="white",
        markeredgewidth=1.5,
        linewidth=2.7,
        label="Restoration site",
        path_effects=foreground_effect,
        zorder=5,
    )

    (synthetic_line,) = ax.plot(
        trajectories["year"],
        trajectories["synthetic"],
        color=colors["synthetic"],
        marker="s",
        markersize=5.5,
        markerfacecolor="white",
        markeredgewidth=1.4,
        linestyle="--",
        linewidth=2.5,
        label="Synthetic control",
        path_effects=foreground_effect,
        zorder=4,
    )

    effect_ax.plot(
        trajectories["year"],
        trajectories["effect"],
        color=colors["effect"],
        marker="o",
        markersize=5.5,
        markerfacecolor="white",
        markeredgewidth=1.4,
        linewidth=2.4,
        zorder=3,
    )
    effect_ax.axhline(0, color="#777777", linewidth=1.1, zorder=1)

    xmin = trajectories["year"].min()
    xmax = trajectories["year"].max()

    for current_ax in (ax, effect_ax):
        # Lightly distinguish the post-intervention period.
        current_ax.axvspan(
            intervention_year,
            xmax + 0.5,
            color=colors["post_period"],
            zorder=0,
        )
        current_ax.axvline(
            intervention_year,
            color=colors["intervention"],
            linestyle="--",
            linewidth=1.3,
            zorder=3,
        )
        current_ax.grid(
            axis="y",
            color=colors["grid"],
            linewidth=0.8,
            zorder=0,
        )
        current_ax.spines[["top", "right"]].set_visible(False)
        current_ax.tick_params(labelsize=10)

    ax.text(
        intervention_year + 0.15,
        0.97,
        f"Intervention: {intervention_year}",
        transform=ax.get_xaxis_transform(),
        va="top",
        ha="left",
        fontsize=9,
        color=colors["intervention"],
    )

    # A symmetric effect scale makes departures from zero easier to compare.
    finite_effects = trajectories["effect"].dropna()
    if not finite_effects.empty:
        effect_limit = max(abs(finite_effects.min()), abs(finite_effects.max()))
        effect_limit = max(effect_limit * 1.10, 0.1)
        effect_ax.set_ylim(-effect_limit, effect_limit)

    ax.set_ylabel(f"{variable.upper()} ({unit})", fontsize=11)
    effect_ax.set_ylabel(f"Site − synthetic ({unit})", fontsize=11)
    effect_ax.set_xlabel("Year", fontsize=11)
    effect_ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    fig.suptitle(
        f"{variable.upper()} restoration impact",
        x=0.10,
        y=0.965,
        ha="left",
        fontsize=15,
    )

    legend_handles = [
        site_line,
        synthetic_line,
        Line2D(
            [0],
            [0],
            color="#2171B5",
            linewidth=2.5,
            label="Weighted donors",
        ),
        Line2D(
            [0],
            [0],
            color=colors["zero_donor"],
            linewidth=1.5,
            alpha=0.5,
            label="Zero-weight donors",
        ),
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.10, 0.925),
        ncol=4,
        frameon=False,
        fontsize=10,
        handlelength=2.8,
        columnspacing=1.8,
    )

    fig.text(
        0.97,
        0.925,
        "Weighted donors: darker = greater weight",
        ha="right",
        va="center",
        fontsize=9,
        color="#52606D",
    )

    ax.set_xlim(xmin - 0.5, xmax + 0.5)

    return fig


def impact_figure(trajectories, intervention_year, variable):
    units = {"svh": "m", "agb": "Mg/ha"}
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
    axes[0].plot(
        trajectories.year, trajectories.observed, "o-", label="Observed boundary mean"
    )
    axes[0].plot(
        trajectories.year, trajectories.synthetic, "o-", label="Synthetic control"
    )
    axes[0].legend()
    axes[0].set_ylabel(f"{variable.upper()} ({units[variable]})")
    axes[1].plot(trajectories.year, trajectories.effect, "o-", color="tab:green")
    axes[1].axhline(0, color="gray", linewidth=1)
    axes[1].set_ylabel(f"Observed − synthetic ({units[variable]})")
    axes[1].set_xlabel("Year")
    for ax in axes:
        ax.axvline(
            intervention_year, color="black", linestyle="--", label="Intervention"
        )
    return fig


def weights_figure(weights):
    selected = weights.loc[weights.weight > 0].sort_values("weight")
    fig, ax = plt.subplots(
        figsize=(8, max(3, 0.22 * len(selected))), layout="constrained"
    )
    ax.barh(selected.donor_id.astype(str), selected.weight)
    ax.set(xlabel="Donor weight", ylabel="Donor ID", title="Positive donor weights")
    return fig
