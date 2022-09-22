import re

# genelist = ['brca', 'brca1', 'brca2', 'palb2', 'palb', 'tp53', 'pten', 'cdh1', 'cdh', 'stk11', 'stk', 'chek2', 'atm']
genelist = ['brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']

def clean_string(notestring):
    """
    This function takes the raw string of the genetic test result section as input
    and cleans it by:
    - removing any remnants of RTF formats (this is no longer necessary after I refined my SQL query)
    - changing all sequences of whitespaces to a single space 
    - stripping whitespace at the edges
    """
    # cleaned_string = re.sub('\\\\par+', '', notestring)
    # cleaned_string = re.sub('\\\\plain\\\\f\d\\\\fs22|\\\\lang1033\\\\hich\\\\f\d\\\\dbch\\\\f\d\\\\loch\\\\f\d\\\\cf\d\\\\fs22|\\\\[ib]\d*', '', cleaned_string)
    cleaned_string = re.sub('\s+', ' ', notestring) # replace all whitespaces with single space
    cleaned_string = cleaned_string.strip() # strip outer whitespace 
    return cleaned_string


def string_to_dict(ir_id, rawstring):
    resultsdict, regiondict = {}, {}
    resultsdict['ir_id'] = ir_id
    regiondict['ir_id'] = ir_id
    ## Find keywords: positive, variants of unknown significance, or negative 
    # notestring = re.sub('\s', '', rawstring).lower() # remove all whitespace because weird spaces are included sometimes
    notestring = rawstring
    matches = re.findall(
        r'positive|likely pathogenic|(?!no )mutations? found|(?!no )mutation identified|negative|no mutations? found|no mutations? detected|no mutations? identified|no pathogenic mutations? found|no brca mutation found|of uncertain significance|hypermethylation', 
        notestring, re.IGNORECASE # brca1(c.68_69del,c.5266dupc)andbrca2(c.5946del)|
    )
    matches_span = [x.span() for x in re.finditer(
        r'positive|likely pathogenic|(?!no )mutations? found|(?!no )mutation identified|negative|no mutations? found|no mutations? detected|no mutations? identified|no pathogenic mutations? found|no brca mutation found|of uncertain significance|hypermethylation', 
        notestring, re.IGNORECASE # brca1(c.68_69del,c.5266dupc)andbrca2(c.5946del)|
    )]
    #### Troubleshooting code ####
    if len(matches)==0:
        print('No known pattern found in section: %s' % notestring)
    ##############################
    # create a list of tuples that are ('match', [span]) 
    match_tuples = []
    for i in range(len(matches)):
        # match_tuples.append((matches[i], matches_span[i]))
        if i < len(matches)-1:
            match_tuples.append((matches[i], notestring[matches_span[i][1]:matches_span[i+1][0]-1]))
        else:
            match_tuples.append((matches[i], notestring[matches_span[i][1]:len(notestring)]))
    # end up with a tuples that looks like [('positive' : ['content1', 'content2']) ('negative' : [content]), ('uncertain' : ['content'])]
    for item in match_tuples:
        if len(item)==0:
            print('Could not extract test results: %s' % rawstring)
        elif len(item)==2:
            if len(item[1])==0: # when gene names not specified 
                if re.match(r'positive|likely pathogenic|(?!no )mutation identified', item[0], re.IGNORECASE): # |brca1(c.68_69del,c.5266dupc)andbrca2(c.5946del)
                    resultsdict['general'] = 'positive'
                elif re.match(r'negative|no mutations? found|(?!no )mutations? found|no mutations? identified|no mutations? detected|no pathogenic mutations? found|no brca mutation found', item[0], re.IGNORECASE):
                    resultsdict['general'] = 'negative'
                elif re.match(r'of uncertain significance', item[0], re.IGNORECASE):
                    resultsdict['general'] = 'variant'
                else:
                    print('Tuple key does not match known pattern: (%s, %s)' % item)
            else:
                for genename in genelist:
                    if genename in item[1].lower():
                        if re.match(r'positive|likely pathogenic|(?!no )mutation identified', item[0], re.IGNORECASE): # |brca1(c.68_69del,c.5266dupc)andbrca2(c.5946del)
                            resultsdict[genename] = 'positive'
                            region = re.findall(r'[ ,(][cp]\.[a-zA-Z0-9\-_>]+[ ;,)(]', item[1])
                            if len(region)>0:
                                regiondict[genename] = str(region)
                        elif re.match(r'negative|no mutations? found|no mutations? identified|(?!no )mutations? found|no mutations? detected|no pathogenic mutations? found|no brca mutation found', item[0], re.IGNORECASE):
                            resultsdict[genename] = 'negative'
                        elif re.match(r'of uncertain significance', item[0], re.IGNORECASE):
                            resultsdict[genename] = 'variant'
                            region = re.findall(r'[ ,(][cp]\.[a-zA-Z0-9\-_>]+[ ;,)(]', item[1])
                            if len(region)>0:
                                regiondict[genename] = str(region)
                        else:
                            print('Tuple key does not match known pattern: (%s, %s)' % item)
    return resultsdict, regiondict
