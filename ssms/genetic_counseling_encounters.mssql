use NM_BI;

/***
The purpose of this query is to find all encounters that were genetic counseling.
I am trying to do this by joining the relevant table with the provider table,
and only selecting rows where provider.title='CGC'.
***/

SELECT COUNT(DISTINCT eo.patient_key) AS n_patients, COUNT(DISTINCT eo.encounter_outpatient_key) as n_encounters
FROM fact.vw_encounter_outpatient eo
	INNER JOIN dim.vw_provider prov
		ON eo.visit_provider_key = prov.provider_key
		--AND prov.is_current=1
		AND prov.title='CGC'

