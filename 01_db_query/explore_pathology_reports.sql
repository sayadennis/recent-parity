use NM_BI;

DECLARE @cohort_id AS INT = 5096;

IF object_id('tempdb..#PATHOLREPORTS') IS NOT NULL
	DROP TABLE #PATHOLREPORTS;

SELECT DISTINCT 
	p.ir_id, p.cerner_central_id, p.clarity_central_id, p.west_mrn, p.first_name, p.last_name, p.birth_date, 
	pc.case_collect_date_key, pcr.report_description, pcrs.section_description, pcrs.note_text AS pathol_note_text
INTO #PATHOLREPORTS
FROM dim.vw_patient_current p
	INNER JOIN fsm_analytics.cohort.cohort_patients cp 
		ON cp.source_ir_id = p.ir_id
            AND cp.is_dltd_flg = 0 
            AND cp.cohort_id = @cohort_id
	INNER JOIN [FSM_Analytics].[pathology_fact].[pathology_case] pc
		ON p.ir_id = pc.patient_ir_id
			AND YEAR(pc.case_collect_date_key)>=2010
			AND YEAR(pc.case_collect_date_key)<=2021
			AND pc.is_canceled = 0
	INNER JOIN [FSM_Analytics].[pathology_dim].[pathology_case_group] pcg
		ON pc.pathology_case_group_key = pcg.pathology_case_group_key
			AND pcg.is_current = 1
			AND pcg.group_desc = 'Surgical Pathology'
	INNER JOIN [FSM_Analytics].[pathology_fact].[pathology_case_report] pcr
		ON pc.pathology_case_key = pcr.pathology_case_key
			AND ((pcr.report_description = 'Surg Path Final Report') 
				OR (pcr.report_description = 'Surg Path Addendum Report'))
	INNER JOIN [FSM_Analytics].[pathology_fact].[pathology_case_report_section] pcrs
		ON pcrs.pathology_case_report_key = pcr.pathology_case_report_key
			AND ((pcrs.section_description = 'Final Diagnosis')
				OR (pcrs.section_description = 'Addendum'))
			AND pcrs.note_text LIKE '%breast%'
			AND (
				(pcrs.note_text LIKE '%tnm%') OR (pcrs.note_text LIKE '%marker%')
			)
	--INNER JOIN fact.vw_note n
	--	ON n.note_key = pcrs.note_key
	--	AND n.source_system_table = 'hno_note_text'
	--	AND note_type IN ('Progress Notes', 'Letter')
	--	and n.note_text LIKE '%breast%'
;

select *
from #PATHOLREPORTS
order by ir_id;

-- 1246 patients in cohort 