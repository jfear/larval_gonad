
def pval_to_string(pvalue):
    if pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "NS"


def format_pval(x, y, pvalue, ax, **kwargs):
    pval_string = pval_to_string(pvalue)
    if pval_string == "NS":
        return ax

    ax.text(x, y, pval_string, ha="center", **kwargs)
    return ax


def add_pvals(x, y, pval, ax, **kwargs):
    """Loops over and adds formatted p-value to plot.
    
    Parameters
    ----------
    x : array-like
        A list of x locations.
    y : array-like
        A list of y locations.
    pval : array-like
        A list of p-values or q-values for plotting.
    ax : plt.Axes
        The axes to add formatted p-values.
    """
    for x_i, y_i, pval_i in zip(x, y, pval):
        format_pval(x_i, y_i, pval_i, ax, **kwargs)


def plot_statsmodels_results(file: str, results: str):
    plt.text(0.01, 0.05, results, {'fontsize': 10}, fontproperties = 'monospace') # approach improved by OP -> monospace!
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(file, dpi=200)
