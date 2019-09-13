import matplotlib.pyplot as plt
import seaborn as sns


def demasculinization(data, ax=None, title=None, legend=False, **kwargs):
    """Stacked barplot common for display of demasculinization.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame organized with chomosome as index, columns as Up/NS/Down,
        and values as proportions.
    ax : plt.Axes, optional
        Alternative axes to draw the plot, by default None
    title : str, optional
        Title to add to the plot, by default None
    legend : bool, optional
        If to keep the legend, by default False

    Returns
    -------
    plt.Axes
        Matplotlib axes with the stacked bar plot.

    Example
    -------
    >>> data = pd.DataFrame({"Up": [.1, .3], "NS": [.7, .4], "Down": [.2, .3]}, index=["X", "A"])
    >>> demasculinization(data)
    """
    if ax is None:
        fig, ax = plt.subplots()

    plot_defaults = dict(
        stacked=True,
        color=["red", "lightgrey", "blue"],
        width=0.9,
        edgecolor="k",
        linewidth=0.2,
    )
    plot_defaults.update(kwargs)

    data.plot.bar(ax=ax, **plot_defaults)

    # Clean up X
    plt.setp(ax.get_xticklabels(), rotation=0)
    ax.set_xlabel("")

    # Clean up Y
    ax.set_ylim(0, 1)
    ax.set_ylabel("")

    # Clean up other parts
    if title is not None:
        ax.set_title(title)

    if not legend:
        ax.legend_ = None

    sns.despine(ax=ax, left=True, bottom=True)

    return ax
