from functools import wraps

def make_figs(fname=None, styles=None, formats=None):
    def _plot_all(func, styles=styles, formats=formats):
        @wraps(func)
        def wrapper(*args, styles=styles, formats=formats, **kwargs):
            def plot_style(style, formats):
                if isinstance(formats, str):
                    formats = [formats]
                elif formats is None:
                    formats = ['png', 'eps']

                with plt.style.context('seaborn-' + style):
                    fn = fname + '_' + style
                    func(*args,  **kwargs)
                    plt.tight_layout()
                    if style != 'notebook':
                        for f in formats:
                            plt.savefig('{}.{}'.format(fn, f))
                        plt.close()

            if isinstance(styles, str):
                styles = [styles]
            elif styles is None:
                styles = ['notebook', 'talk', 'poster', 'paper']

            for style in styles:
                plot_style(style, formats)

        return wrapper
    return _plot_all
