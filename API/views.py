# from snippets.models import Snippet
# from snippets.serializers import SnippetSerializer
from django.http import Http404
from rest_framework.response import Response
from rest_framework.decorators import api_view
from rest_framework import status
import cx_Oracle
import pandas as pd
import json
import copy
from django.views.decorators.csrf import csrf_exempt
from .serializers import FindingSerializer, InitFindingSerializer
from API.utils import extract

def get_stats(group):
    return {'min': group.min(), 'max': group.max()}

# onto_df = pd.read_pickle("API/static/data/ontology.pkl")
# Get the dose range for each study
info_df = pd.read_pickle("API/static/data/animals_per_group_per_sex.pkl")
info_df.study_id = info_df.study_id.astype(int).astype(str)
range_df = info_df[info_df.dose > 0]
range_df = range_df.groupby(('study_id')).dose.apply(get_stats).unstack().reset_index()[['study_id','max','min']]
range_df.columns = ['study_id','dose_max','dose_min']

# Load dataframes with information on studies, compounds and findings
compound_df = pd.read_pickle("API/static/data/compound.pkl")
findings_df = pd.read_pickle("API/static/data/findings.pkl.gz", compression='gzip')
findings_df.study_id = findings_df.study_id.astype(int).astype(str)
study_df = pd.read_pickle("API/static/data/study.pkl")
study_df.study_id = study_df.study_id.astype(int).astype(str)
organ_onto_df = pd.read_pickle("API/static/data/organ_ontology.pkl")
observation_onto_df = pd.read_pickle("API/static/data/observation_ontology.pkl")

study_df = pd.merge(study_df, range_df,
                    how='left', on='study_id', left_index=False, right_index=False,
                    sort=False)
study_cmpd_df = pd.merge(study_df[['study_id', 'subst_id', 'normalised_administration_route',
                            'normalised_species', 'normalised_strain', 'source_company',
                            'exposure_period_days', 'report_number']],
                         compound_df[['subst_id', 'smiles', 'status', 'common_name',
                            'cas_number','pharmacological_action']],
                         how='left', on='subst_id', left_index=False, right_index=False,
                         sort=False)
all_df = pd.merge(study_cmpd_df,
                  findings_df[['study_id','source', 'observation_normalised', 'grade', 
                                'organ_normalised', 'dose', 'relevance', 'normalised_sex']],
                  how='left', on='study_id', left_index=False, right_index=False,
                  sort=False)

# Get the summary of data per single substance
study_count_df = all_df[['subst_id', 'normalised_species', 'study_id']].groupby(['subst_id', 'normalised_species']).study_id.nunique().reset_index()
study_count_df.columns = ['subst_id', 'normalised_species', 'study_count']
study_count_df['count'] = study_count_df.normalised_species + ': ' + study_count_df.study_count.astype(int).astype(str)
study_count_df = study_count_df[['subst_id', 'count']].groupby('subst_id').agg(lambda x : '\n'.join(x)).reset_index()

# Merge target/action
withAction = compound_df[compound_df.action.notnull()]
withAction['targetAction'] = withAction.target+' : '+withAction.action
noAction = compound_df[~compound_df.action.notnull()]
noAction['targetAction'] = noAction.target
target_df = pd.concat([withAction, noAction], ignore_index=True)[['subst_id', 'targetAction']]
target_df = target_df.groupby('subst_id').agg(lambda x : '\n'.join(x)).reset_index()

# Generate output dataframe
subst_info_df = pd.merge(target_df, study_count_df, 
                how='left', on='subst_id', left_index=False, right_index=False, 
                sort=False)
subst_info_df = subst_info_df.drop_duplicates()


@api_view(['GET'])
def source(request):
    global find_df

    host = ''
    port = ''
    sid = ''
    user = ''
    password = ''

    """conn = connectDB(host, port, sid, user, password)
    cursor = conn.cursor()"""

    # Generate normalised study dataframe
    """cmd = "SELECT DISTINCT PARENT_LUID AS study_id, RELEVANCE, \
        STANDARDISED_PATHOLOGY, STANDARDISED_ORGAN, DOSE \
        FROM HISTOPATHOLOGICALFI \
        WHERE STANDARDISED_PATHOLOGY IS NOT NULL \
        AND STANDARDISED_ORGAN IS NOT NULL"
    cursor.execute(cmd)
    results = cursor.fetchall()
    tmp_table = []
    for (study_id, relevance, observation, organ, dose) in results:
        tmp_table.append([study_id, relevance, observation, organ, dose])
        print (str(study_id)+","+str(relevance)+","+str(observation)+","+str(organ)+","+str(dose))
    find_df = pd.DataFrame(tmp_table, columns=['study_id', 'relevance',
                                               'observation_normalised', 'organ_normalised', 'dose'])"""

    pass

@api_view(['GET'])
def initFindings(request):
    global all_df, compound_df, subst_info_df

    num_studies = len(all_df.study_id.unique().tolist())
    num_structures = len(all_df.subst_id.unique().tolist())

    optionsDict = {}

    optionsDict['sources'] = all_df.source.dropna().unique().tolist()
    optionsDict['sources'].sort()

    # optionsDict['grades'] = all_df.grade.dropna().unique().tolist()
    # optionsDict['grades'].sort()

    optionsDict['routes'] = all_df.normalised_administration_route.dropna().unique().tolist()
    optionsDict['routes'].sort()

    optionsDict['sex'] = all_df.normalised_sex.dropna().unique().tolist()
    optionsDict['sex'].sort()

    optionsDict['species'] = all_df.normalised_species.dropna().unique().tolist()
    optionsDict['species'].sort()

    optionsDict['pharmacological_action'] = all_df.pharmacological_action.dropna().unique().tolist()
    optionsDict['pharmacological_action'].sort()

    optionsDict['compound_name'] = all_df.common_name.dropna().unique().tolist()
    optionsDict['compound_name'].sort()

    optionsDict['cas_number'] = all_df.cas_number.dropna().unique().tolist()
    optionsDict['cas_number'].sort()

    exposure_range = all_df.exposure_period_days.dropna().unique().tolist()
    exposure_range.sort()
    optionsDict['exposure_min'] = int(exposure_range[0])
    optionsDict['exposure_max'] = int(exposure_range[-1])

    optionsDict['organs'] = {}
    optionsDict['observations'] = {}
    for source in optionsDict['sources']:
        organs = all_df[all_df.source.str.lower() == source.lower()].organ_normalised.dropna().unique().tolist()
        # Create nested dictionary for angular treeviews
        organs_df = organ_onto_df[organ_onto_df.child_term.str.lower().isin([x.lower() for x in organs])]
        organs_df = getValuesForTree(organs_df,organ_onto_df)
        relations = organs_df.groupby(by='parent_term')['child_term'].apply(list).to_dict()
        parents = set(relations.keys()) & set(organ_onto_df[organ_onto_df.level == 1].child_term.tolist())
        optionsDict['organs'][source] = create_dictionary(relations, parents)

        observations = all_df[all_df.source.str.lower() == source.lower()].observation_normalised.dropna().unique().tolist()
        # Create nested dictionary for angular treeviews
        observations_df = observation_onto_df[observation_onto_df.child_term.str.lower().isin([x.lower() for x in observations])]
        observations_df = getValuesForTree(observations_df,observation_onto_df)
        relations = observations_df.groupby(by='parent_term')['child_term'].apply(list).to_dict()
        parents = set(relations.keys()) & set(observation_onto_df[observation_onto_df.level == 1].child_term.tolist())
        optionsDict['observations'][source] = create_dictionary(relations, parents)

    ##############
    # Pagination #
    ##############
    page = 1
    init = 0
    total = len(compound_df)
    if total < 5:
        end = total
        num_pages = total
    else:
        end = 5
        num_pages = int(total / 5)

    # Range of pages to show
    previous = 0
    nexts = 7
    range_pages = range(1, num_pages + 1)[previous:nexts]
    previous_page = page - 1
    next_page = page + 1

    if (next_page >= num_pages):
        next_page = 0

    output = pd.merge(compound_df, subst_info_df, 
                    how='left', on='subst_id', left_index=False, right_index=False, 
                    sort=False)
    output.drop(['target', 'action'], axis=1, inplace=True)
    output = output.drop_duplicates()
    output.common_name = output.common_name.str.replace(', ', '\n')
    output.pharmacological_action = output.pharmacological_action.str.replace(', ', '\n')
    output = output[init:end].fillna(value="-").to_dict('records')

    results = {
        'data': output,
        'allOptions': optionsDict,
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page,
        'num_studies': num_studies,
        'num_structures': num_structures
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def findings(request):

    global all_df, study_df

    filtered_tmp = all_df[:]

    # Relevancy
    relevant = request.GET.get("treatmentRelated")
    if relevant:
        filtered_tmp = filtered_tmp[filtered_tmp.relevance == 'Treatment related']

    # Sex
    sex = request.GET.get("sex")
    if sex:
        filtered_tmp = filtered_tmp[filtered_tmp.normalised_sex == sex]

    # Exposure
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    if min_exposure and max_exposure:
        filtered_tmp = filtered_tmp[(filtered_tmp.exposure_period_days >= int(min_exposure)) &
                    (filtered_tmp.exposure_period_days <= int(max_exposure))]

    queryDict = {}
    # Pharmacological action
    all_pharm = request.GET.getlist("pharmacological_action")
    if len(all_pharm) > 0:
        queryDict['pharmacological_action'] = 'pharmacological_action == @all_pharm'

    # Pharmacological action
    all_compound_name = request.GET.getlist("compound_name")
    if len(all_compound_name) > 0:
        queryDict['compound_name'] = 'common_name == @all_compound_name'

    # CAS number
    all_cas_number = request.GET.getlist("cas_number")
    if len(all_cas_number) > 0:
        queryDict['cas_number'] = 'cas_number == @all_cas_number'

    # Administration route
    all_routes = request.GET.getlist("routes")
    if len(all_routes) > 0:
        queryDict['routes'] = 'normalised_administration_route == @all_routes'

    # Species
    all_species = request.GET.getlist("species")
    if len(all_species) > 0:
        queryDict['species'] = 'normalised_species == @all_species'

    ##
    ## Filter organs, observations and grades by category
    ##

    # Organs
    all_organs = request.GET.getlist("organs")
    if len(all_organs) > 0:
        all_organs = all_organs[0].split(', ')
        tmp_dict = {}
        for v in all_organs:
            category, val = v.split(' | ')
            if category not in tmp_dict:
                tmp_dict[category] = [val]
            else:
                tmp_dict[category].append(val)
        queryList = []
        for category in tmp_dict:
            tmp_list = '[%s]' %(', '.join(['\'%s\'' %x.strip() for x in tmp_dict[category]]))
            queryList.append('(source == \'%s\' and organ_normalised == %s)' %(category.strip(), tmp_list))
        queryDict['organs'] = ' and '.join(list(queryList))

    # Observations
    all_observations = request.GET.getlist("observations")
    if len(all_observations) > 0:
        all_observations = all_observations[0].split(', ')
        tmp_dict = {}
        for v in all_observations:
            category, val = v.split(' | ')
            if category not in tmp_dict:
                tmp_dict[category] = [val]
            else:
                tmp_dict[category].append(val)
        queryList = []
        for category in tmp_dict:
            tmp_list = '[%s]' %(', '.join(['\'%s\'' %x.strip() for x in tmp_dict[category]]))
            queryList.append('(source == \'%s\' and observation_normalised == %s)' %(category.strip(), tmp_list))
        queryDict['observations'] = ' and '.join(list(queryList))

    # Grade
    # all_grades = request.GET.getlist("grade")
    # if len(all_grades) > 0:
    #     all_grades = all_grades[0].split(', ')
    #     tmp_dict = {}
    #     for v in all_grades:
    #         category, val = v.split(' | ')
    #         if category not in tmp_dict:
    #             tmp_dict[category] = [val]
    #         else:
    #             tmp_dict[category].append(val)
    #     queryList = []
    #     for category in tmp_dict:
    #         tmp_list = '[%s]' %(', '.join(['\'%s\'' %x.strip() for x in tmp_dict[category]]))
    #         queryList.append('(source == \'%s\' and grade == %s)' %(category.strip(), tmp_list))
    #     queryDict['grades'] = ' and '.join(list(queryList))

    #####################
    # Apply all filters #
    #####################
    query_string = ''
    if queryDict != {}:
        query_string = ' and '.join(list(queryDict.values()))
        filtered = filtered_tmp.query(query_string)
    else:
        filtered = filtered_tmp[:]

    sources = filtered.source.dropna().unique().tolist()

    optionsDict = {}
    if not filtered.empty:
        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('organs', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['organs'] = {}
        for source in sources:
            organs = tmp_df[tmp_df.source.str.lower() == source.lower()].organ_normalised.dropna().unique().tolist()
            # Create nested dictionary for angular treeviews
            organs_df = organ_onto_df[organ_onto_df.child_term.str.lower().isin([x.lower() for x in organs])]
            organs_df = getValuesForTree(organs_df, organ_onto_df)
            relations = organs_df.groupby(by='parent_term')['child_term'].apply(list).to_dict()
            parents = set(relations.keys()) & set(organ_onto_df[organ_onto_df.level == 1].child_term.tolist())
            optionsDict['organs'][source] = create_dictionary(relations, parents)

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('observations', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['observations'] = {}
        for source in sources:
            observations = tmp_df[tmp_df.source.str.lower() == source.lower()].observation_normalised.dropna().unique().tolist()
            # Create nested dictionary for angular treeviews
            observations_df = observation_onto_df[observation_onto_df.child_term.str.lower().isin([x.lower() for x in observations])]
            observations_df = getValuesForTree(observations_df, observation_onto_df)
            relations = observations_df.groupby(by='parent_term')['child_term'].apply(list).to_dict()
            parents = set(relations.keys()) & set(observation_onto_df[observation_onto_df.level == 1].child_term.tolist())
            optionsDict['observations'][source] = create_dictionary(relations, parents)

        # tmp_dict = copy.deepcopy(queryDict)
        # tmp_dict.pop('grade', None)
        # valuesL = list(tmp_dict.values())
        # if len(valuesL) > 0:
        #     query_string = ' and '.join(valuesL)
        #     tmp_df = filtered_tmp.query(query_string)
        # else:
        #     tmp_df = filtered_tmp
        # optionsDict['grade'] = tmp_df.grade.dropna().unique().tolist()
        # optionsDict['grade'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('pharmacological_action', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['pharmacological_action'] = tmp_df.pharmacological_action.dropna().unique().tolist()
        optionsDict['pharmacological_action'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('cas_number', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['cas_number'] = tmp_df.cas_number.dropna().unique().tolist()
        optionsDict['cas_number'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('compound_name', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['compound_name'] = tmp_df.common_name.dropna().unique().tolist()
        optionsDict['compound_name'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('routes', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['routes'] = tmp_df.normalised_administration_route.dropna().unique().tolist()
        optionsDict['routes'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('species', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['species'] = tmp_df.normalised_species.dropna().unique().tolist()
        optionsDict['species'].sort()

        exposure_range = filtered.exposure_period_days.dropna().unique().tolist()
        exposure_range.sort()
        optionsDict['exposure_min'] = int(exposure_range[0])
        optionsDict['exposure_max'] = int(exposure_range[-1])

        optionsDict['sex'] = all_df.normalised_sex.dropna().unique().tolist()
        optionsDict['sex'].sort()

    #############
    # Aggregate #
    #############

    num_studies = filtered.study_id.nunique()
    num_structures = filtered.subst_id.nunique()
    
    output = pd.merge(filtered[['subst_id']], compound_df, 
                    how='left', on='subst_id', left_index=False, right_index=False, 
                    sort=False)
    output = pd.merge(output, subst_info_df, 
                    how='left', on='subst_id', left_index=False, right_index=False, 
                    sort=False)
    output.drop(['target', 'action'], axis=1, inplace=True)
    output = output.drop_duplicates()
    output.common_name = output.common_name.str.replace(', ', '\n')
    output.pharmacological_action = output.pharmacological_action.str.replace(', ', '\n')

    ##############
    # Pagination #
    ##############

    page = int(request.GET.get("page"))

    if page != 0:
        init = (int(page) - 1) * 5
        end = init + 5
    else:
        init = 0
        end = len(filtered)

    num_pages = int(len(output) / 5)

    # Range of pages to show
    if page < 4:
        previous = 0
        nexts = 7
    elif page > (num_pages - 4):

        previous = num_pages - 7
        nexts = num_pages
    else:
        previous = page - 4
        nexts = page + 3
    range_pages = range(1, num_pages + 1)[previous:nexts]
    previous_page = page - 1
    next_page = page + 1

    if (next_page >= num_pages):
        next_page = 0

    # output = filtered[init:end].fillna(value="-").to_dict('records')
    output = output[init:end].fillna(value="-").to_dict('records')

    results = {
        'data': output,
        'allOptions': optionsDict,
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page,
        'num_studies': num_studies,
        'num_structures': num_structures
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def qualitative(request):
    global findings_df, study_df

    all_organs = request.GET.getlist("organs")
    all_observation = request.GET.getlist("observations")
    all_species = request.GET.getlist("species")
    all_routes = request.GET.getlist("routes")
    sex = request.GET.getlist("sex")
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    relevant = request.GET.get("treatmentRelated")

    #############
    # Filter    #
    #############

    filter_dict = {}
    if min_exposure:
        filter_dict['min_exposure'] = int(min_exposure[0])
    if max_exposure:
        filter_dict['max_exposure'] = int(max_exposure[0])
    if len(all_routes) > 0:
        filter_dict['route'] = all_routes
    if len(all_species) > 0:
        filter_dict['species'] = all_species

    filtered = extract.filter_study(filter_dict, study_df)
    filtered = findings_df[findings_df.study_id.isin(filtered.study_id)]

    filtered = pd.merge(filtered[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                  'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days', 'report_number']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    if sex:
        filtered = filtered[filtered.normalised_sex.str.lower() == sex[0].lower()]
    if relevant:
        filtered = filtered[filtered.relevance == 'treatment related']
    if len(all_organs) > 0:
        filter_dict['organs'] = all_organs
        filtered = filtered[filtered['organ_normalised'].isin(all_organs)]
    if len(all_observation) > 0:
        filter_dict['observations'] = all_observation
        filtered = filtered[filtered['observation_normalised'].isin(all_observation)]

    ###################################
    # Get stats for relevant findings #
    ###################################

    # Get the number of studies per substance
    count_df = filtered.groupby(('subst_id')).study_id.nunique().to_frame().reset_index()
    # Get the global dose range per substance
    range_df = filtered[filtered.dose > 0]
    range_df = range_df.groupby(('subst_id')).dose.apply(extract.get_stats).unstack().reset_index()
    # Get all stats into a single dataframe
    stats_df = pd.merge(count_df, range_df, how='inner', on='subst_id',
                        left_index=False, right_index=False, sort=False)
    stats_df.columns = ['subst_id', 'study_count', 'dose_max', 'dose_min']

    # Expand
    filtered_piv = extract.expand(filtered, filter_dict, onto_df)

    ######################################
    # Aggragate by substance and finding #
    ######################################
    # Define finding as organ+observation
    filtered_piv.organ_normalised = filtered_piv.organ_normalised.fillna('NA')
    filtered_piv.observation_normalised = filtered_piv.observation_normalised.fillna('NA')
    filtered_piv['finding'] = filtered_piv.apply(
        lambda row: row.organ_normalised + '_' + row.observation_normalised, axis=1)
    filtered_piv = filtered_piv[['subst_id', 'finding', 'dose']]

    # Aggregate by substance and finding (as defined above), keeping the minimum dose
    # for each substance/finding instance
    group_df = filtered_piv.groupby(('subst_id', 'finding')).min().add_prefix('min_').reset_index()

    #######################################
    # Pivot so that each finding is a row #
    #######################################
    ### Qualitative
    group_df.loc[:, 'min_dose'] = 1
    pivotted_df = group_df.pivot_table(index='subst_id', columns='finding', values='min_dose').reset_index()
    pivotted_df['is_active'] = 'True'
    qualitative_df = pd.merge(stats_df, pivotted_df, how='left', on='subst_id',
                              left_index=False, right_index=False, sort=False)
    qualitative_df.is_active = qualitative_df.is_active.fillna('False')
    # Reorder columns
    cols = qualitative_df.columns.tolist()
    cols = cols[0:4] + [cols[-1]] + cols[4:-1]
    qualitative_df = qualitative_df[cols]

    results = {
        "data": qualitative_df.fillna(value=0)
    }

    return Response(results)

@api_view(['GET'])
def quantitative(request):
    global findings_df, study_df

    all_organs = request.GET.getlist("organs")
    all_observation = request.GET.getlist("observations")
    all_species = request.GET.getlist("species")
    all_routes = request.GET.getlist("routes")
    sex = request.GET.getlist("sex")
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    relevant = request.GET.get("treatmentRelated")

    global findings_df, study_df

    #############
    # Filter    #
    #############

    filter_dict = {}
    if min_exposure:
        filter_dict['min_exposure'] = int(min_exposure[0])
    if max_exposure:
        filter_dict['max_exposure'] = int(max_exposure[0])
    if len(all_routes) > 0:
        filter_dict['route'] = all_routes
    if len(all_species) > 0:
        filter_dict['species'] = all_species

    filtered = extract.filter_study(filter_dict, study_df)
    filtered = findings_df[findings_df.study_id.isin(filtered.study_id)]

    filtered = pd.merge(filtered[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                  'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days', 'report_number']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    if sex:
        filtered = filtered[filtered.normalised_sex.str.lower() == sex[0].lower()]
    if relevant:
        filtered = filtered[filtered.relevance == 'treatment related']
    if len(all_organs) > 0:
        filter_dict['organs'] = all_organs
        filtered = filtered[filtered['organ_normalised'].isin(all_organs)]
    if len(all_observation) > 0:
        filter_dict['observations'] = all_observation
        filtered = filtered[filtered['observation_normalised'].isin(all_observation)]

    ###################################
    # Get stats for relevant findings #
    ###################################

    # Get the number of studies per substance
    count_df = filtered.groupby(('subst_id')).study_id.nunique().to_frame().reset_index()
    # Get the global dose range per substance
    range_df = filtered[filtered.dose > 0]
    range_df = range_df.groupby(('subst_id')).dose.apply(extract.get_stats).unstack().reset_index()
    # Get all stats into a single dataframe
    stats_df = pd.merge(count_df, range_df, how='inner', on='subst_id',
                        left_index=False, right_index=False, sort=False)
    stats_df.columns = ['subst_id', 'study_count', 'dose_max', 'dose_min']

    # Expand
    filtered_piv = extract.expand(filtered, filter_dict, onto_df)

    ######################################
    # Aggragate by substance and finding #
    ######################################
    # Define finding as organ+observation
    filtered_piv.organ_normalised = filtered_piv.organ_normalised.fillna('NA')
    filtered_piv.observation_normalised = filtered_piv.observation_normalised.fillna('NA')
    filtered_piv['finding'] = filtered_piv.apply(
        lambda row: row.organ_normalised + '_' + row.observation_normalised, axis=1)
    filtered_piv = filtered_piv[['subst_id', 'finding', 'dose']]

    # Aggregate by substance and finding (as defined above), keeping the minimum dose
    # for each substance/finding instance
    group_df = filtered_piv.groupby(('subst_id', 'finding')).min().add_prefix('min_').reset_index()

    #######################################
    # Pivot so that each finding is a row #
    #######################################
    ### Quantitative
    pivotted_df = group_df.pivot_table(index='subst_id', columns='finding', values='min_dose').reset_index()
    pivotted_df['is_active'] = 'True'
    quantitative_df = pd.merge(stats_df, pivotted_df, how='left', on='subst_id',
                               left_index=False, right_index=False, sort=False)
    quantitative_df.is_active = quantitative_df.is_active.fillna('False')
    # Reorder columns
    cols = quantitative_df.columns.tolist()
    cols = cols[0:4] + [cols[-1]] + cols[4:-1]
    quantitative_df = quantitative_df[cols]

    results = {
        "data": quantitative_df.fillna(value=0)
    }

    return Response(results)

@api_view(['GET'])
def study(request):
    id = request.GET.get("id")

    res = pd.merge(study_df[study_df.study_id == id], compound_df,
                   how='left', on='subst_id', left_index=False,
                   right_index=False, sort=False)

    results = {
        'study': res.fillna(value="-")
    }
    return Response(results)

# @api_view(['GET'])
# def connectDB(request):

#     host= request.GET.get("host")
#     port= int(request.GET.get("port"))
#     sid= request.GET.get("sid")
#     user= request.GET.get("user")
#     password= request.GET.get("password")

#     # try:
#     #     dsn_tns = cx_Oracle.makedsn(host, port, sid)
#     #     conn = cx_Oracle.connect(user, password, dsn=dsn_tns)
#     # except:
#     #     #  cx_Oracle.DatabaseError as e:
#     #     conn = None
#     #     # error, = e.args
#     #     # connStatus = error.message
#     #     connStatus = 'Failed'
#     # else:
#     #     connStatus = 'Connected'

#     connStatus = 'Connected'

#     results = {
#         'connStatus': connStatus
#     }

#     send_data = FindingSerializer(results, many=False).data

#     return Response(send_data)

def create_dictionary(relations, parents):
    """
    Recursive function to create dictionary for treeviews
    """
    dict_out = {}
    for key in parents:
        dict_aux = {}
        if key in relations:
            dict_return = create_dictionary(relations, relations[key])
            dict_aux[key] = dict_return
        else:
            dict_aux[key] = {}
        dict_out = mergeDeepDict(dict_out, dict_aux)
    return dict_out

def mergeDeepDict(d1, d2):
    """
    Recursive function to merge two nested dictionaries
    """
    # update first dict with second recursively
    for k, v in d1.items():
        if k in d2:
            d2[k] = mergeDeepDict(v, d2[k])
    d1.update(d2)
    return d1

def getValuesForTree(df_filter,onto_tree_df):

    columns_name = ['child_term', 'parent_term', 'ontology']
    df_filter = df_filter[columns_name]
    tree_created = False
    last_size = -1
    while not tree_created:

        df_filter = pd.merge(df_filter, onto_tree_df, how='left', left_on='parent_term', right_on='child_term')

        df_filter = df_filter.dropna()
        df_1 = df_filter[['child_term_x', 'parent_term_x', 'ontology_x']]
        df_2 = df_filter[['child_term_y', 'parent_term_y', 'ontology_y']]

        df_1.columns = columns_name
        df_2.columns = columns_name

        df_filter = pd.concat([df_1, df_2])

        df_filter.drop_duplicates(inplace=True)

        if last_size == len(df_filter):
            tree_created = True

        last_size = len(df_filter)

    return df_filter

@api_view(['GET'])
def plot(request):
    global all_df, compound_df

    filtered_tmp = all_df[:]

    # Relevancy
    relevant = request.GET.get("treatmentRelated")
    if relevant:
        filtered_tmp = filtered_tmp[filtered_tmp.relevance == 'Treatment related']

    # Sex
    sex = request.GET.get("sex")
    if sex:
        filtered_tmp = filtered_tmp[filtered_tmp.normalised_sex == sex]

    # Exposure
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    if min_exposure and max_exposure:
        filtered_tmp = filtered_tmp[(filtered_tmp.exposure_period_days >= int(min_exposure)) &
                                    (filtered_tmp.exposure_period_days <= int(max_exposure))]

    queryDict = {}
    # Pharmacological action
    all_pharm = request.GET.getlist("pharmacological_action")
    if len(all_pharm) > 0:
        queryDict['pharmacological_action'] = 'pharmacological_action == @all_pharm'

    # Pharmacological action
    all_compound_name = request.GET.getlist("compound_name")
    if len(all_compound_name) > 0:
        queryDict['compound_name'] = 'common_name == @all_compound_name'

    # CAS number
    all_cas_number = request.GET.getlist("cas_number")
    if len(all_cas_number) > 0:
        queryDict['cas_number'] = 'cas_number == @all_cas_number'

    # Administration route
    all_routes = request.GET.getlist("routes")
    if len(all_routes) > 0:
        queryDict['routes'] = 'normalised_administration_route == @all_routes'

    # Species
    all_species = request.GET.getlist("species")
    if len(all_species) > 0:
        queryDict['species'] = 'normalised_species == @all_species'

    ##
    ## Filter organs, observations and grades by category
    ##

    # Organs
    all_organs = request.GET.getlist("organs")
    if len(all_organs) > 0:
        all_organs = all_organs[0].split(', ')
        tmp_dict = {}
        for v in all_organs:
            category, val = v.split(' | ')
            if category not in tmp_dict:
                tmp_dict[category] = [val]
            else:
                tmp_dict[category].append(val)
        queryList = []
        for category in tmp_dict:
            tmp_list = '[%s]' % (', '.join(['\'%s\'' % x.strip() for x in tmp_dict[category]]))
            queryList.append('(source == \'%s\' and organ_normalised == %s)' % (category.strip(), tmp_list))
        queryDict['organs'] = ' and '.join(list(queryList))

    # Observations
    all_observations = request.GET.getlist("observations")
    if len(all_observations) > 0:
        all_observations = all_observations[0].split(', ')
        tmp_dict = {}
        for v in all_observations:
            category, val = v.split(' | ')
            if category not in tmp_dict:
                tmp_dict[category] = [val]
            else:
                tmp_dict[category].append(val)
        queryList = []
        for category in tmp_dict:
            tmp_list = '[%s]' % (', '.join(['\'%s\'' % x.strip() for x in tmp_dict[category]]))
            queryList.append('(source == \'%s\' and observation_normalised == %s)' % (category.strip(), tmp_list))
        queryDict['observations'] = ' and '.join(list(queryList))

    # Grade
    # all_grades = request.GET.getlist("grade")
    # if len(all_grades) > 0:
    #     all_grades = all_grades[0].split(', ')
    #     tmp_dict = {}
    #     for v in all_grades:
    #         category, val = v.split(' | ')
    #         if category not in tmp_dict:
    #             tmp_dict[category] = [val]
    #         else:
    #             tmp_dict[category].append(val)
    #     queryList = []
    #     for category in tmp_dict:
    #         tmp_list = '[%s]' %(', '.join(['\'%s\'' %x.strip() for x in tmp_dict[category]]))
    #         queryList.append('(source == \'%s\' and grade == %s)' %(category.strip(), tmp_list))
    #     queryDict['grades'] = ' and '.join(list(queryList))

    #####################
    # Apply all filters #
    #####################
    query_string = ''
    if queryDict != {}:
        query_string = ' and '.join(list(queryDict.values()))
        filtered = filtered_tmp.query(query_string)
    else:
        filtered = filtered_tmp[:]

    num_studies = len(filtered.study_id.unique().tolist())
    num_structures = len(filtered.subst_id.unique().tolist())
    sources = filtered.source.dropna().unique().tolist()

    optionsDict = {}
    if not filtered.empty:
        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('organs', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['organs'] = {}
        for source in sources:
            organs = tmp_df[tmp_df.source == source].organ_normalised.dropna().unique().tolist()
            # Create nested dictionary for angular treeviews
            organs_df = organ_onto_df[organ_onto_df.child_term.isin(organs)]
            organs_df = getValuesForTree(organs_df, organ_onto_df)
            relations = organs_df.groupby(by='parent_term')['child_term'].apply(list).to_dict()
            parents = set(relations.keys()) & set(organ_onto_df[organ_onto_df.level == 1].child_term.tolist())
            optionsDict['organs'][source] = create_dictionary(relations, parents)

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('observations', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['observations'] = {}
        for source in sources:
            observations = tmp_df[tmp_df.source == source].observation_normalised.dropna().unique().tolist()
            # Create nested dictionary for angular treeviews
            observations_df = observation_onto_df[observation_onto_df.child_term.isin(observations)]
            observations_df = getValuesForTree(observations_df, observation_onto_df)
            relations = observations_df.groupby(by='parent_term')['child_term'].apply(list).to_dict()
            parents = set(relations.keys()) & set(
                observation_onto_df[observation_onto_df.level == 1].child_term.tolist())
            optionsDict['observations'][source] = create_dictionary(relations, parents)

        # tmp_dict = copy.deepcopy(queryDict)
        # tmp_dict.pop('grade', None)
        # valuesL = list(tmp_dict.values())
        # if len(valuesL) > 0:
        #     query_string = ' and '.join(valuesL)
        #     tmp_df = filtered_tmp.query(query_string)
        # else:
        #     tmp_df = filtered_tmp
        # optionsDict['grade'] = tmp_df.grade.dropna().unique().tolist()
        # optionsDict['grade'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('pharmacological_action', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['pharmacological_action'] = tmp_df.pharmacological_action.dropna().unique().tolist()
        optionsDict['pharmacological_action'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('cas_number', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['cas_number'] = tmp_df.cas_number.dropna().unique().tolist()
        optionsDict['cas_number'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('compound_name', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['compound_name'] = tmp_df.common_name.dropna().unique().tolist()
        optionsDict['compound_name'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('routes', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['routes'] = tmp_df.normalised_administration_route.dropna().unique().tolist()
        optionsDict['routes'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('species', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['species'] = tmp_df.normalised_species.dropna().unique().tolist()
        optionsDict['species'].sort()

        exposure_range = filtered.exposure_period_days.dropna().unique().tolist()
        exposure_range.sort()
        optionsDict['exposure_min'] = int(exposure_range[0])
        optionsDict['exposure_max'] = int(exposure_range[-1])

        optionsDict['sex'] = all_df.normalised_sex.dropna().unique().tolist()
        optionsDict['sex'].sort()


    plot_info = filtered.groupby(['normalised_species'])['normalised_species'].count()
    x = plot_info.index
    y = plot_info.values

    results = {
        'x': x,
        'y': y,
        'allOptions': optionsDict,
        'num_studies': num_studies,
        'num_structures': num_structures
    }

    return Response(results)
