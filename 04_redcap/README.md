# Import data to REDCap for chart review & export following data entry for analysis

## Import to REDCap 

We first design the data entry fields on REDCap, then download the data dictionary in a CSV format. (`FrequencyAndResultsOfGeneticTe_ImportTemplate_2022-05-31.csv` - the one where features are in columns.) We then use this data dictionary to clean up our data at hand to not only match the field names but also to encode multiple-choice fields to integer mappings. 

## Export from REDCap 

Following the completion of manual data entry on REDCap, we download the completed dataset from REDCap in a CSV format. In order for further analysis, we will construct a JSON-format map describing the integer mappings of multiple-choice entries on REDCap. (`redcap_datadict_to_json.py`) In order to do this, we download a **different** type of data dictionary from REDCap, one where the features are in rows. 
