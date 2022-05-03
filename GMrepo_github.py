import requests
import json
from pandas.core.frame import DataFrame


Health = 'D006262'


BA = 1680


# BA Health-------------------------------------------------------------------

data_query = {'mesh_id': Health, "ncbi_taxon_id": BA}  ## -- to get statistics on MeSH ID D006262
url = 'https://gmrepo.humangut.info/api/getMicrobeAbundancesByPhenotypeMeshIDAndNCBITaxonID'
data = requests.post(url, data=json.dumps(data_query))

## --get DataFrames
abundance_and_meta_data = DataFrame(data.json().get('abundance_and_meta_data'))
list(abundance_and_meta_data)

abundance_and_meta_data.to_csv("AKK_Health_abundance_and_meta_data.csv", index_label="index_label")


