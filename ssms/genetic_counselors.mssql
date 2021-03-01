use NM_BI;

--This is the list of genetic counselors (title CGC)
--I think it includes counselors that rarely see breast cancer patients too.
--The list that Theresa gave me had 7 counselors:
--Sciaraffa, Park, Fallen, Szymaniak, Hyde, DeGreef, Williams

select ir_id, full_name, first_name, last_name, specialty_fpsc, primary_department
from dim.vw_provider
where title='CGC'
	and is_current=1;
