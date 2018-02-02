"""Set of helper functions for the notebook."""
import os
from pathlib import Path
from datetime import datetime
from yaml import load

import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
from IPython import get_ipython

from .plotting import add_styles


class Nb(object):
    def __init__(self, nb_name=None, project_dir=None, config_dir=None,
                 ref_dir=None, fig_dir=None, table_dir=None,
                 subproject_dir=None, cache=None, formats=None, styles=None,
                 styles_wide=None, watermark=None, **kwargs):
        """Helper method for working consistently in notebook.

        Stores a set a bunch of useful attributes. Turns on a bunch of commonly
        used notebook magics. If matplotlib stylelib exists in the config_dir
        then it imports user defined styles.

        Parameters
        ----------
        nb_name : str
            Name of the current notebook.
        project_dir : str
            Name of the project directory.
        config_dir : str
            Name of the config directory.
        ref_dir : str
            Name of the references directory.
        subproject_dir : str
            Name of the subproject directory for placing output.
        fig_dir : str
            Name of the figures directory.
        table_dir : str
            Name of the tables directory.
        cache : str
            Name of the cache directory.
        formats : str or list
            Default list of formats to use for plotting. For example 'png' or
            ['png', 'svg'].
        styles : str or list
            Default list of matplotlib.style.library to use for plotting. For
            example 'seaborn-notebook' or ['seaborn-notebook',
            'seaborn-paper'].
        watermark : bool
            If true turn on watermarking.
        **kwargs
            Additional arguments that are stored as attributes

        Attributes
        ----------
        nb_name : str
            Name of the current notebook.
        project_dir : str
            Name of the project directory.
        config_dir : str
            Name of the config directory.
        ref_dir : str
            Name of the references directory.
        fig_dir : str
            Name of the figures directory.
        table_dir : str
            Name of the tables directory.
        cache : str
            Name of the cache directory.
        formats : str or list
            Default list of formats to use for plotting. For example 'png' or
            ['png', 'svg'].
        styles : str or list
            Default list of matplotlib.style.library to use for plotting. For
            example 'seaborn-notebook' or ['seaborn-notebook',
            'seaborn-paper'].
        styles_wide : str or list
            Default list of matplotlib.style.library to use for plotting wide
            (two column) images. For example 'seaborn-notebook' or
            ['seaborn-notebook', 'seaborn-paper'].
        date : str
            Current date, generated upon creation.
        fasta : str
            Path to fasta file.
        chromsizes : str
            Path to chromsizes file.
        gtf : str
            Path to gtf file.
        gtf_db : str
            Path to gtf_db file.
        annot : str
            Path to annot file.
        syn : str
            Path to syn file.
        seurat : Seurat
            Useful Seurat paths.

        """
        self.nb_name = nb_name
        self.project_dir = project_dir
        self.config_dir = config_dir
        self.ref_dir = ref_dir
        self.fig_dir = fig_dir
        self.table_dir = table_dir
        self.cache = cache
        self.formats = formats
        self.styles = styles
        self.styles_wide = styles_wide
        self.date = datetime.now().strftime("%Y-%m-%d")

        # Add useful paths
        assembly = kwargs['assembly']
        tag = kwargs['tag']
        self.fasta = os.path.join(self.ref_dir, assembly, tag, 'fasta',
                                  f'{assembly}_{tag}.fasta')

        self.chromsizes = os.path.join(self.ref_dir, assembly, tag, 'fasta',
                                       f'{assembly}_{tag}.chromsizes')

        self.gtf = os.path.join(self.ref_dir, assembly, tag, 'gtf',
                                f'{assembly}_{tag}.gtf')

        self.gtf_db = os.path.join(self.ref_dir, assembly, tag, 'gtf',
                                   f'{assembly}_{tag}.gtf.db')

        self.annot = os.path.join(self.ref_dir, assembly, tag,
                                  'fb_annotation',
                                  f'{assembly}_{tag}.fb_annotation')

        self.syn = os.path.join(self.ref_dir, assembly, tag,
                                'fb_synonym',
                                f'{assembly}_{tag}.fb_synonym')

        self.seurat = Seurat(subproject_dir)

        # Add useful mappers
        _annot = pd.read_csv(self.annot, sep='\t', index_col=1)
        self.fbgn2symbol = _annot['gene_symbol'].to_dict()
        self.symbol2fbgn = {v: k for k, v in self.fbgn2symbol.items()}

        self.fbgn2chrom = pd.read_csv(
            os.path.join(self.project_dir, 'output/fbgn2chrom.tsv'),
            sep='\t', index_col=0)

        # Add Colors
        self.colors = sns.color_palette('Paired', n_colors=12)
        sns.set_palette(self.colors)

        self.color_female = sns.xkcd_rgb['rose pink']
        self.color_male = sns.xkcd_rgb['dodger blue']
        self.colors_sex = sns.color_palette([self.color_female,
                                             self.color_male])

        self.color_c1 = sns.xkcd_rgb['dusty purple']
        self.color_c2 = sns.xkcd_rgb['golden rod']

        # Add any key word args
        self._config_attrs = kwargs.keys()
        for k, v in kwargs.items():
            setattr(self, k, v)

        # turn on magics
        self._start_magics(watermark=watermark)

        # Set up plotting
        self._setup_plotting()

        # Turn off scientific notation
        np.set_printoptions(precision=5, suppress=True)

    def _start_magics(self, watermark=None):
        """Start up the notebook magics I commonly use."""
        mgc = get_ipython().magic

        # Activate the autoreload extension for easy reloading of external
        # packages
        mgc('reload_ext autoreload')
        mgc('autoreload 2')

        # Trun on the water mark
        if watermark:
            mgc('reload_ext watermark')
            mgc('watermark -u -d -g')

        # Plot inline
        mgc('matplotlib inline')

    def _setup_plotting(self):
        styles = os.path.join(self.config_dir, 'stylelib')
        if os.path.exists(styles):
            add_styles(styles)

        #sns.set_context('notebook')
        mpl.style.use(['default', 'notebook'])

    @classmethod
    def setup_notebook(cls, nb_name=None, config_name='common.yml',
                       watermark=True, **kwargs):
        """Helper function to consistently setup notebooks.

        Functions detects working folder and sets up a larval_gonad.notebook.Nb
        object with sane defaults.

        Parameters
        ----------
        nb_name : str
            Name of the current notebook.
        config_name : str
            Name of the a YAML config file.
        watermark : bool
            If truen then output watermark information.
        kwargs
            Additional arguments to pass to Nb.

        """
        if nb_name is None:
            cache_dir = 'cache'
        else:
            cache_dir = os.path.join('cache', nb_name)

        # Figure out current project, config, and references folder
        prj = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '../../')
        )
        cfg = os.path.join(prj, 'config')
        ref = os.environ.get('REFERENCES_DIR', None)

        # set defaults
        defaults = {
            'nb_name': nb_name,
            'project_dir': prj,
            'config_dir': cfg,
            'ref_dir': ref,
            'fig_dir': './figures',
            'table_dir': './tables',
            'cache': os.path.join(prj, 'output', cache_dir),
            'formats': ['png', 'pdf', 'svg'],
            'styles': ['notebook', 'paper', 'talk', 'poster'],
            'styles_wide': ['notebook-wide', 'paper-wide', 'talk-wide',
                            'poster-wide'],
            'watermark': watermark
        }

        defaults.update(kwargs)

        # Import external config
        fname = os.path.join(cfg, config_name)
        with open(fname) as fh:
            config = load(fh)
        defaults.update(config)

        return cls(**defaults)

    def fig_name(self, fname):
        if self.nb_name is not None:
            fname = '_'.join([self.nb_name, fname])

        return os.path.join(self.fig_dir, fname)

    def table_name(self, fname):
        if self.nb_name is not None:
            fname = '_'.join([self.nb_name, fname])

        return os.path.join(self.table_dir, fname)

    def __repr__(self):
        return str(self)

    def __str__(self):
        keys = ['nb_name', 'project_dir', 'config_dir', 'fig_dir', 'table_dir',
                'cache', 'formats', 'styles', 'styles_wide', 'date']
        keys.extend(self._config_attrs)
        res = []
        for key in keys:
            value = self.__dict__[key]
            if value is None:
                value = 'None'
            res.append('{}:\t{}'.format(key, value))

        return '\n'.join(res)


class Seurat(object):
    """Class that stores basic paths for files from Seurat analysis."""
    def __init__(self, path=None):
        """Create a list of Seurat paths.
        raw : str
            Path to raw.
        scaled : str
            Path to scaled data.
        dispersion : str
            Path to dispersion estimates.
        normalized_read_counts : str
            Path to normalized read counts.
        principal_components_cell : str
            Path to PCA cell loadings.
        principal_components_gene : str
            Path to PCA gene loadings.
        tsne : str
            Path to tsne.
        biomarkers : str
            Path to biomarkers list.
        clusters : str
            Path to clusters identities.
        robj : str
            Path to seurat R object.

        """

        if path is None:
            self.raw = None
            self.scaled = None
            self.dispersion = None
            self.normalized_read_counts = None
            self.principal_components_cell = None
            self.principal_components_gene = None
            self.tsne = None
            self.biomarkers = None
            self.clusters = None
            self.robj = None
        else:
            Path(path).mkdir(exist_ok=True)
            self.raw = os.path.join(path, 'raw.tsv')
            self.scaled = os.path.join(path, 'scaled.tsv')
            self.dispersion = os.path.join(path, 'dispersion.tsv')
            self.normalized_read_counts = os.path.join(
                path, 'normalized_read_counts.tsv'
            )
            self.principal_components_cell = os.path.join(
                path, 'principal_components_cell.tsv'
            )
            self.principal_components_gene = os.path.join(
                path, 'principal_components_gene.tsv'
            )
            self.tsne = os.path.join(path, 'tsne.tsv')
            self.biomarkers = os.path.join(path, 'biomarkers.tsv')
            self.clusters = os.path.join(path, 'clusters.tsv')
            self.robj = os.path.join(path, 'seurat.Robj')

    def get_raw(self):
        df = pd.read_csv(self.raw, sep='\t')
        df.index.name = 'FBgn'
        df.columns.name = 'cell_id'
        return df

    def get_scaled(self):
        df = pd.read_csv(self.scaled, sep='\t')
        df.index.name = 'FBgn'
        df.columns.name = 'cell_id'
        return df

    def get_dispersion(self):
        df = pd.read_csv(self.dispersion, sep='\t')
        df.index.name = 'FBgn'
        return df

    def get_normalized_read_counts(self):
        df = pd.read_csv(self.normalized_read_counts, sep='\t')
        df.index.name = 'FBgn'
        return df

    def get_principal_components_cell(self):
        df = pd.read_csv(self.principal_components_cell, sep='\t')
        df.index.name = 'cell_id'
        df.columns.name = 'PC'
        return df


    def get_principal_components_gene(self):
        df = pd.read_csv(self.principal_components_gene, sep='\t')
        df.index.name = 'FBgn'
        df.columns.name = 'PC'
        return df

    def get_tsne(self):
        df = pd.read_csv(self.tsne, sep='\t')
        df.index.name = 'FBgn'
        return df

    def get_biomarkers(self):
        df = pd.read_csv(self.biomarkers, sep='\t', index_col='gene')
        df.index.name = 'FBgn'
        return df

    def get_clusters(self):
        df = pd.read_csv(self.clusters, sep='\t')
        df.index.name = 'cell_id'
        df.columns = ['cluster']
        return df.cluster
