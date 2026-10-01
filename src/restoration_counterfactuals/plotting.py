"""Return maps and figures; callers decide how to display and save them."""

import matplotlib.pyplot as plt


def matching_map(site, region, treatments, donors, masks, exclusions):
    import geemap

    m = geemap.Map(basemap="HYBRID")
    m.centerObject(site, 12)
    m.addLayer(region, {"color": "888888"}, "Donor search area")
    m.addLayer(exclusions, {"color": "ff8800"}, "Other project exclusions")
    for name, mask in masks.items():
        m.addLayer(
            mask.selfMask().clip(region),
            {"palette": ["66bb66"]},
            name + " (first treatment pixel)",
            False,
        )
    m.addLayer(site, {"color": "ff0000"}, "Project boundary")
    m.addLayer(treatments, {"color": "ffff00"}, "Treatment samples")
    m.addLayer(donors, {"color": "0000ff"}, "Selected donors")
    return m


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
