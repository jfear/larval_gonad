from functools import wraps

import matplotlib as mpl
import matplotlib.pyplot as plt


def add_styles(dirname):
    mpl.style.core.USER_LIBRARY_PATHS.append(dirname)
    mpl.style.core.update_user_library(mpl.style.library)
    mpl.style.reload_library()


def make_ax(*args, **kwargs):
    fig, ax = plt.subplots(*args, **kwargs)
    return ax


def figure_element(func):
    """Decorate figure elements to create an Axes if none is given.

    I define a figure element as a function that imports, munges, and plots. A
    figure panel can be built by combining multiple figure elements using
    matplotlib.Axes. However, some times I will only want to plot the element,
    so this wrapper will generate a new Axes and pass it to the function for
    me.

    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        ax = kwargs.get('ax', None)
        if ax is None:
            fig, ax = plt.subplots()
            kwargs.update({'ax': ax})
        return func(*args, **kwargs)


def make_figs(fname=None, styles=None, formats=None, layout=True,
              kws_layout=None):
    if isinstance(formats, str):
        formats = [formats]
    elif formats is None:
        formats = ['png', 'eps']

    if isinstance(styles, str):
        styles = [styles]
    elif styles is None:
        styles = ['notebook']

    if kws_layout is None:
        kws_layout = {}

    def _plot_all(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            def plot_style(style, formats):
                with plt.style.context(['common', style]):
                    func(*args, **kwargs)
                    if layout:
                        plt.tight_layout(**kws_layout)
                    if ('notebook' not in style) & (fname is not None):
                        fn = fname + '_' + style
                        for f in formats:
                            plt.savefig('{}.{}'.format(fn, f))
                        plt.close()

            for style in styles:
                plot_style(style, formats)

        return wrapper
    return _plot_all


def plot_barcode_rank(umi, selected=None, title=None, **kwargs):
    """Plot Barcode Rank Plot."""
    options = {
        'kind': 'scatter',
        's': 3,
        'logx': True,
        'logy': True
    }
    options.update(kwargs)
    try:
        ax = options.pop('ax')
    except KeyError:
        ax = make_ax()

    dat = umi.to_frame()
    dat.columns = ['umi']

    dat = dat.sort_values('umi', ascending=False)
    dat['cell_num'] = list(range(1, dat.shape[0] + 1))

    dat.plot('cell_num', 'umi', c='lightgrey', ax=ax, **options)

    if selected is not None:
        dat.loc[dat.index.isin(selected), :].plot('cell_num', 'umi', c='g', ax=ax, **options)

    if title is not None:
        ax.set_title(title)
