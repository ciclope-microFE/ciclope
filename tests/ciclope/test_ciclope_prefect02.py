import numpy as np
import pandas as pd
import prefect
from prefect import task, Flow, Parameter
from ciclope import recon_utils as ru
from skimage.filters import threshold_otsu, gaussian
from scipy import ndimage, misc
import pandas
from io import StringIO

# input_file = Parameter("DATA_FILE", default='/home/gianthk/Data/2019.001.coop_TUberlin_simulierte_Mensch.iorig/trabecular_sample_mini3/2000L_crop_imgaussfilt_60micron_uint8_0000.tif')
master_file = 'test_ciclope_flow.csv'
# files = ['/home/gianthk/Data/2019.001.coop_TUberlin_simulierte_Mensch.iorig/trabecular_sample_mini3/2000L_crop_imgaussfilt_60micron_uint8_0000.tif', '/home/gianthk/Data/StefanFly_test/test_00__rec/test_00__rec_mini2/recon_0000.tif']


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
def read_master(filename: str) -> pd.DataFrame:
    # Load the data set into a pandas dataframe.
    logger = prefect.context.get("logger")
    logger.info("Loading master file: %s" % filename)
    return pandas.read_csv(StringIO(filename))

@task
def load(filein: str) -> np.ndarray:
    logger = prefect.context.get("logger")
    logger.info("Reading image: %s" % filein)
    return ru.read_tiff_stack(filein)

@task
def smooth(image: np.ndarray, smooth=0) -> np.ndarray:
    if smooth == 1:
        logger = prefect.context.get("logger")
        logger.info("Apply gaussian smooth")
        image = gaussian(image, sigma=4, preserve_range=True)

    return image

def ciclope_flow():
    with Flow("ciclope-flow") as flow:
        """
        2) Map and Loop through ciclopes selected for run.
        """

        filename = Parameter('filename', required=True)
        s = Parameter('smooth')
        # vs = Parameter('vs')
        # r = Parameter('r')
        # t = Parameter('t')

        # I = load(filein=input_file)
        # I = load.map(files)
        I = load.map(filename)
        I = smooth.map(image=I, smooth=s)

        print("here")
        # # An optional parameter "people", with a default list of names
        # people = Parameter("people", default=["Arthur", "Ford", "Marvin"])
        # # Map `say_hello` across the list of names
        # say_hello.map(people)

    return flow

if __name__ == '__main__':
    df = pd.read_csv(script_path(master_file))
    df = df.where(pd.notnull(df), None)

    flow = ciclope_flow()

    parameters = {'filename': df['filein'][df['run'] == 1].tolist(),
                  'smooth': df['smooth'][df['run'] == 1].tolist(),
                  # 'vs': df['vs'][df['run'] == 1].tolist(),
                  # 'r': df['r'][df['run'] == 1].tolist(),
                  # 't': df['t'][df['run'] == 1].tolist()
                  }

    state = flow.run(parameters)

    type(state._result.value)
    task_ref = flow.get_tasks()[1]
    ru.plot_midplanes(state.result[task_ref]._result.value[0])



# class ciclopeFlow(FlowSpec):
#     """
#     Flow of ciclope runs.
#     The flow performs the following steps:
#     1) Ingests a CSV into a Pandas Dataframe.
#     2) Select ciclopes for execution.
#     3) Launch ciclope with given args.
#     4) Save something..
#     """
#
#     master = IncludeFile(
#         "master_table",
#         help="Path to master table.",
#         default=script_path("test_ciclope_flow.csv"),
#     )
#
#     @step
#     def start(self):
#         """
#         The start step:
#         1) Loads master with metadata into pandas dataframe.
#         2) Select ciclopes marked for run.
#         3) Run all ciclopes.
#         """
#         import pandas
#         from io import StringIO
#
#         # Load the data set into a pandas dataframe.
#         self.dataframe = pandas.read_csv(StringIO(self.master))
#
#         # The column 'run' is a run flag. Get only ciclopes marked for run.
#         # self.dataframe = self.dataframe[self.dataframe['run'] == 1]
#         self.torun = self.dataframe.index[self.dataframe['run'] == 1]
#
#         # We want to compute some statistics for each genre. The 'foreach'
#         # keyword argument allows us to compute the statistics for each genre in
#         # parallel (i.e. a fan-out).
#         self.next(self.run_ciclope, foreach="torun")
#
#     @step
#     def run_ciclope(self):
#         """
#         Run one ciclope step (imresize).
#         """
#         # The ciclope currently being processed is a class property called
#         # 'input'.
#         self.id = self.input
#         print("Running ciclope n. %s" % self.id)
#
#         # read data
#         print("Reading image: %s" % self.dataframe["filein"][self.id])
#         I = ru.read_tiff_stack(self.dataframe["filein"][self.id])
#
#         # gaussian smooth
#         if self.dataframe["smooth"][self.id] == 1:
#             print("Apply gaussian smooth")
#             I = gaussian(I, sigma=1, preserve_range=True)
#
#         if self.dataframe["r"][self.id] != 0:
#             I = ndimage.zoom(I, 1 / self.dataframe["r"][self.id], output=None, order=2)
#
#         self.dataframe["note"][self.id] = "processed"
#         self.next(self.join)
#
#     @step
#     def join(self, inputs):
#         self.results = [input.id for input in inputs]
#         self.next(self.end)
#
#     @step
#     def end(self):
#         print('\n'.join(self.results))
