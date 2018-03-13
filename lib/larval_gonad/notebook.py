"""Set of helper functions for the notebook."""
import os
from pathlib import Path
from datetime import datetime
from subprocess import check_output

import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
from IPython import get_ipython

from .config import config, PROJECT_DIR, CONFIG_DIR, CACHE_DIR, REFERENCES_DIR
from .plotting import add_styles


class Nb(object):
    def __init__(self, nb_name=None, project_dir=None, subproject_dir=None,
                 seurat_dir=None, config_dir=None, ref_dir=None, fig_dir=None,
                 formats=None, styles=None, styles_wide=None,
                 styles_full=None, watermark=None, **kwargs):
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
        seurat_dir : str
            Name of the  directory with seurat output.
        fig_dir : str
            Name of the figures directory.
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
        subproject_dir : str
            Directory to save outputs from this subproject.
        seurat_dir : str
            Location of Seurat output. Will default to the subproject_dir.
        config_dir : str
            Name of the config directory.
        ref_dir : str
            Name of the references directory.
        fig_dir : str
            Name of the figures directory.
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
        styles_full : str or list
            Default list of matplotlib.style.library to use for plotting wide
            (two column) images. For example 'seaborn-notebook' or
            ['seaborn-notebook', 'seaborn-paper'].
        date : str
            Current date, generated upon creation.
        conda_env : str
            Name of the current conda environment location.
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
        self.subproject_dir = subproject_dir
        self.seurat_dir = seurat_dir
        self.config_dir = config_dir
        self.ref_dir = ref_dir
        self.fig_dir = fig_dir
        self.formats = formats
        self.styles = styles
        self.styles_wide = styles_wide
        self.styles_full = styles_full
        self.date = datetime.now().strftime("%Y-%m-%d")
        self.conda_env = self.get_conda()

        # Add useful reference paths
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

        if (seurat_dir is None):
            self.seurat = None
        else:
            self.seurat = Seurat(seurat_dir)

        # Add useful mappers
        _annot = pd.read_csv(self.annot, sep='\t', index_col=1)
        self.fbgn2symbol = _annot['gene_symbol'].to_dict()
        self.symbol2fbgn = {v: k for k, v in self.fbgn2symbol.items()}

        try:
            self.fbgn2chrom = pd.read_csv(
                os.path.join(self.project_dir, 'output/fbgn2chrom.tsv'),
                sep='\t', index_col=0)
        except Exception:
            print(
                'Please check output/fbgn2chrom.tsv. '
                'If it does not exist, run bin/fbgn2chrom.py'
            )

        # Add Colors
        colors = kwargs.pop('colors')
        self.colors = colors['cycle']

        self.color_chrom = colors['chrom']
        self.color_xa = [colors['chrom'][5], colors['chrom'][0]]

        self.color_female = colors['female']
        self.color_male = colors['male']
        self.colors_sex = sns.color_palette([self.color_female,
                                             self.color_male])

        self.color_c1 = colors['c1']
        self.color_c2 = colors['c2']
        self.color_c3 = colors['c3']

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

        mpl.style.use(['common', 'notebook'])
        sns.set_palette(self.colors)
        mpl.rcParams['pdf.fonttype'] = 42

    def get_conda(self):
        conda_info = check_output(['conda', 'info']).decode('utf-8')
        for x in conda_info.split('\n'):
            if 'envs directories' in x:
                return x.split(':')[1].strip()

    @classmethod
    def setup_notebook(cls, nb_name=None, subproject_dir=None, seurat_dir=None,
                       watermark=True, **kwargs):
        """Helper function to consistently setup notebooks.

        Functions detects working folder and sets up a larval_gonad.notebook.Nb
        object with sane defaults.

        Parameters
        ----------
        nb_name : str
            Name of the current notebook.
        subproject_dir : str
            Directory to save outputs from this subproject.
        seurat_dir : str
            Location of Seurat output. Will default to the subproject_dir.
        watermark : bool
            If truen then output watermark information.
        kwargs
            Additional arguments to pass to Nb.

        """
        # Set seurat_dir to subproject_dir if it was None.
        if subproject_dir is None:
            subproject_dir = Path(PROJECT_DIR, 'output').as_posix()

        fig_dir = Path(subproject_dir, 'figures')
        fig_dir.mkdir(parents=True, exist_ok=True)

        # set defaults
        defaults = {
            'nb_name': nb_name,
            'project_dir': PROJECT_DIR,
            'subproject_dir': subproject_dir,
            'seurat_dir': seurat_dir,
            'config_dir': CONFIG_DIR,
            'ref_dir': REFERENCES_DIR,
            'fig_dir': fig_dir.as_posix(),
            'formats': ['png', 'pdf'],
            'styles': ['notebook', 'talk'],
            'watermark': watermark
        }

        defaults.update(kwargs)

        # Import external config
        defaults.update(config)

        # Add wide and full styles
        _styles = defaults['styles']
        defaults['styles_wide'] = [x + '-wide' for x in _styles],
        defaults['styles_full'] = [x + '-full' for x in _styles],

        return cls(**defaults)

    def fig_name(self, fname):
        if self.nb_name is not None:
            fname = '_'.join([self.nb_name, fname])

        return os.path.join(self.fig_dir, fname)

    def table_name(self, fname):
        if self.nb_name is not None:
            fname = '_'.join([self.nb_name, fname])

        return os.path.join(self.subproject_dir, fname)

    def __repr__(self):
        return str(self)

    def __str__(self):
        keys = ['nb_name', 'project_dir', 'config_dir', 'fig_dir',
                'subproject_dir', 'seurat_dir', 'formats', 'styles',
                'styles_wide', 'styles_full', 'date']
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
        metadata : str
            Path to metadata.
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

        if (path is None) | Path(path).exists():
            self.raw = None
            self.metadata = None
            self.scaled = None
            self.dispersion = None
            self.normalized_read_counts = None
            self.principal_components_cell = None
            self.principal_components_gene = None
            self.principal_components_stdev = None
            self.tsne = None
            self.biomarkers = None
            self.clusters = None
            self.robj = None
        else:
            self.raw = os.path.join(path, 'raw.tsv')
            self.metadata = os.path.join(path, 'metadata.tsv')
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
            self.principal_components_stdev = os.path.join(
                path, 'principal_components_stdev.tsv'
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

    def get_metadata(self):
        df = pd.read_csv(self.metadata, sep='\t')
        df.index.name = 'cell_id'
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

    def get_principal_components_stdev(self):
        df = pd.read_csv(self.principal_components_stdev, sep='\t')
        df.index.name = 'PC'
        df.columns.name = 'stdev'
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
