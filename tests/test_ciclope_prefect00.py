import numpy as np
import prefect
from prefect import task, Flow, Parameter
from ciclope import recon_utils as ru
from skimage.filters import threshold_otsu, gaussian
from scipy import ndimage, misc

# DATA_FILE = Parameter("DATA_FILE", default="/home/gianthk/Data/2019.001.coop_Tuberlin_simulierte_Mensch.iorig/trabecular_sample_mini3/2000L_crop_imgaussfilt_60micron_uint8_0000.tif")
input_file = Parameter("DATA_FILE", default='/home/gianthk/Data/2019.001.coop_TUberlin_simulierte_Mensch.iorig/trabecular_sample_mini3/2000L_crop_imgaussfilt_60micron_uint8_0000.tif')

def script_path(filename):
    """
    A convenience function to get the absolute path to a file in this
    tutorial's directory. This allows the tutorial to be launched from any
    directory.
    """
    import os

    filepath = os.path.join(os.path.dirname(__file__))
    return os.path.join(filepath, filename)

@task
def loop(master):
    """
    The start step:
    1) Loads master with metadata into pandas dataframe.
    2) Loop through ciclopes
    3) Launch ciclopes marked for run.
    """
    import pandas
    from io import StringIO

    logger = prefect.context.get("logger")
    logger.info("Hello world!")

    # Load the data set into a pandas dataframe.
    df = pandas.read_csv(StringIO(master))

@task
def load(filein: str):
    logger = prefect.context.get("logger")
    logger.info("Reading image: %s" % filein)
    print("load")
    return ru.read_tiff_stack(filein)

@task
def smooth(I):
    logger = prefect.context.get("logger")
    logger.info("Apply gaussian smooth")
    return gaussian(I, sigma=1, preserve_range=True)

with Flow("hello-flow") as flow:
    #master = script_path("test_ciclope_flow.csv")

    I = load(filein=input_file)

    I = smooth(I=I)

    print("here")
    # # An optional parameter "people", with a default list of names
    # people = Parameter("people", default=["Arthur", "Ford", "Marvin"])
    # # Map `say_hello` across the list of names
    # say_hello.map(people)

state = flow.run()
type(state._result.value)
task_ref = flow.get_tasks()[2]
# ru.plot_midplanes(state.result[task_ref]._result.value)


# # gaussian smooth
# if self.dataframe["smooth"][self.id] == 1:
#     print("Apply gaussian smooth")
