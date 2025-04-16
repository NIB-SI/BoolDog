''' API to download and use SBML-qual models from bio.models directly. '''

import logging
import re
from pathlib import Path
import requests


try:
    from biomodels_restful_api_client import services as bmservices
    _BIOMODELS_AVAILABLE = True
except ImportError as e:
    raise ImportError(
        'biomodels-restful-api-client '
        '(https://bitbucket.org/biomodels/biomodels-resftful-api-client/) '
        'is needed to fetch models through the biomodels API. '
        'We suggest you install it using pip. ') from e

logger = logging.getLogger(__name__)

MAMO_ACCESSIONS = [
    "MAMO_0000030",
    "MAMO_0000053"
]

EXAMPLE_MODEL_ID = "BIOMD0000000562"


def fetch_model(model_id, sbml_file=None, local_file=None, check_modelling_approach=False):
    '''Fetch SBML from BioModels by model id

    Parameters
    ----------

    model_id: str
        BioModels model identifier (e.g. 'BIOMD0000000562')

    sbml_file: str
        Name of the file containing the model in BioModels. Optional, if not
        given, will be determined from the model info.


    local_file: str or path-like
        Name of path to save the downloaded model to. Optional, if not given
        will save the model to the a file in the cwd with the name as on the
        remote server.

    check_modelling_approach: Bool
        Whether to check modelling annotation "modellingApproach" falls into
        Boolean or logical modelling.

    Returns
    -------

    local_file: str
        Name of the local file containing the download.

    Notes
    --------
    This first collects the model info/metadata, and uses the 'files' --> 'main' attribute
    to find the first file with 'SBML'/'sbml' in the beginning of the description.

    Will always overwrite an existing file.

    '''

    model_info = bmservices.get_model_info(model_id)


    if check_modelling_approach:
        modelling_approach = model_info["modellingApproach"]["accession"]
        if not modelling_approach in MAMO_ACCESSIONS:
            raise ValueError(
                f'Model approach {modelling_approach} is not in logical or Boolean (model info).')

    if not model_info['format']['name'].lower() == 'sbml':
        raise ValueError(
            f'Model {model_id} is not in SBML format (model info).')


    sbml_file = None

    for f in model_info['files']['main']:
        if re.match('^sbml', f['description'].lower()) is not None:
            sbml_file = f['name']

    if sbml_file is None:
        raise ValueError(
            f'Model {model_id} is not in SBML format (file description).')

    local_file = _download(model_id, filename=sbml_file, local_file=local_file)

    logger.info("Saved model %s from biomodels to %s.", model_id, local_file)

    return local_file


# GET /model/download/{model_id}
def _download(model_id, filename=None, local_file=None):
    '''Downloads a file from a model in BioModels

    Parameters
    ----------

    model_id: str
        BioModels model identifier (e.g. )

    Returns
    -------

    local_file: str
        Name of the local file containing the download.

    '''
    download_url = bmservices.API_URL + "/model/download/" + model_id

    if filename is not None:
        response = requests.get(download_url + "?filename=" + filename,
                                timeout=10)
    else:
        response = requests.get(download_url, timeout=10)

    # Determine local file name if not given
    if local_file is None:
        if filename is not None:
            local_file = Path(filename)
        else:
            # make up a name for the entire archive
            local_file = Path(f"{model_id}.omex")

    # Save the file data to the local file
    with open(local_file, 'wb') as file:
        file.write(response.content)

    return local_file
