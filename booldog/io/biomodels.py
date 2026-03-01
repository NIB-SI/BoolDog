''' API to download and use SBML-qual models from bio.models directly. '''

import logging
import re
from pathlib import Path
import requests

logger = logging.getLogger(__name__)

BIOMODELS_BASE_URL = "https://www.biomodels.org/"

MAMO_ACCESSIONS = ["MAMO_0000030", "MAMO_0000053"]
'''list: MAMO accessions for models that are compatible with BoolDog.

The list contains the following accessions:
* MAMO_0000030: logical model (http://identifiers.org/mamo/MAMO_0000030)
* MAMO_0000053: Boolean model (http://identifiers.org/mamo/MAMO_0000053)
'''

EXAMPLE_MODEL_ID = "BIOMD0000000562"
'''str: BioModels model identifier for an example SBML-qual model (Chaouiya2013
- EGF and TNFalpha mediated signalling pathway). See
https://www.ebi.ac.uk/biomodels/BIOMD0000000562 for more details.
'''

def fetch_model(model_id,
                sbml_file=None,
                local_file=None,
                check_modelling_approach=False):
    '''Fetch SBML from BioModels by model id

    Parameters
    ----------

    model_id: str
        BioModels model identifier (e.g. 'BIOMD0000000562')

    sbml_file: str
        Name of the file containing the model in the BioModels registry.
        Optional, if not given, will be determined from the model info.

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

    model_info = fetch_model_info(model_id)

    if check_modelling_approach:
        modelling_approach = model_info["modellingApproach"]["accession"]
        if not modelling_approach in MAMO_ACCESSIONS:
            raise ValueError(
                f'Model approach {modelling_approach} is not in logical or Boolean (model info).'
            )

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


# get model info
def fetch_model_info(model_id):
    '''Fetch model info/metadata from BioModels by model id.

    Parameters
    ----------
    model_id: str
        BioModels model identifier (e.g. 'BIOMD0000000562')
    Returns
    -------
    dict
        Model info/metadata as a dictionary.
    '''
    response = requests.get(f"{BIOMODELS_BASE_URL}/{model_id}?format=json", timeout=30)
    return response.json()

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
    download_url = f"{BIOMODELS_BASE_URL}/model/download/{model_id}"
    if filename is not None:
        download_url += f"?filename={filename}"

    response = requests.get(download_url, timeout=30)

    if not response.ok:
        raise ValueError(
            f"Failed to download model {model_id} from BioModels. "
            f"HTTP status code: {response.status_code}, URL: {download_url}"
        )

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
