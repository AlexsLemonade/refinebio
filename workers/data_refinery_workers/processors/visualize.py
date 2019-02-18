import copy
import time

import datashader as ds
import pandas as pd
import numpy as np
import holoviews as hv
import matplotlib as mpl
from holoviews import opts
from datashader import transfer_functions as tf
from dask import dataframe as daskdf

import multiprocessing

from holoviews.operation.datashader import datashade, shade, dynspread, rasterize
from holoviews.operation import decimate
from holoviews import opts

from data_refinery_common.logging import get_and_configure_logger

hv.extension('bokeh')
logger = get_and_configure_logger(__name__)


def visualize(input_frame, output_path, width=2000, height=2000, logz=True, backend='bokeh'):
    """
    Generate high resolution visualizations of supplied Pandas DataFrames

    """
    try:
        dask_df = daskdf.from_pandas(input_frame, npartitions=multiprocessing.cpu_count()).persist()
        logger.info("Converted to Dask Frame..")

        num_genes, num_samples = input_frame.shape
        da = dask_df.to_dask_array(True).persist()
        logger.info("Converted to Dask Array..")

        # Visualize
        if backend == 'matplotlib':
            img = hv.Image((np.arange(num_samples), np.arange(num_genes), da))
            rasterized_img = rasterize(img)
            viridis_bad_black = copy.copy(mpl.cm.get_cmap('viridis')) # copy the default cmap
            viridis_bad_black.set_bad((0,0,0))
            rasterized_img.opts(cmap=viridis_bad_black, dpi=400, logz=True)
            hv.save(rasterized_img, output_path, backend=backend)
        else:
            img = hv.Image((np.arange(num_samples), np.arange(num_genes), da))
            rasterized_img = rasterize(img, width=width, height=height)

            viridis_bad_black = copy.copy(mpl.cm.get_cmap('viridis')) # copy the default cmap
            viridis_bad_black.set_bad((1,0,1))

            rasterized_img.opts(width=width, height=height, cmap=viridis_bad_black, logz=True)
            hv.save(rasterized_img, output_path, backend='bokeh')

        logger.info("Output visualization!", output_path=output_path)
        return output_path
    except Exception as e:
        logger.exception("Unable to visualize dataframe!",
            output_path=output_path,
            width=width,
            height=height,
            logz=logz,
            backend=backend
        )
        return None
