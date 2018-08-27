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
import math
from django.views.decorators.csrf import csrf_exempt
from .serializers import FindingSerializer,Pageserializer
from API.utils import extract

def get_stats(group):
    return {'min': group.min(), 'max': group.max()}

# Get the dose range for each study
output_df= None
info_df = pd.read_pickle("API/static/data/animals_per_group_per_sex.pkl")
info_df.study_id = info_df.study_id.astype(int).astype(str)
range_df = info_df[info_df.dose > 0]
range_df = range_df.groupby(('study_id')).dose.apply(get_stats).unstack().reset_index()[['study_id','max','min']]
range_df.columns = ['study_id','dose_max','dose_min']

# Load dataframes with information on studies, compounds and findings
compound_df = pd.read_pickle("API/static/data/compound.pkl")
findings_df = pd.read_pickle("API/static/data/findings.pkl.gz", compression='gzip')
findings_df.study_id = findings_df.study_id.astype(int).astype(str)
findings_df.loc[(findings_df.relevance != 'Treatment related'), 'relevance'] = 'Not related'
study_df = pd.read_pickle("API/static/data/study.pkl")
study_df.study_id = study_df.study_id.astype(int).astype(str)
organ_onto_df = pd.read_pickle("API/static/data/organ_ontology.pkl")
observation_onto_df = pd.read_pickle("API/static/data/observation_ontology.pkl")

# Merge target/action
withAction = compound_df[compound_df.action.notnull()]
withAction.loc[:,'targetAction'] = withAction.target+' : '+withAction.action
noAction = compound_df[~compound_df.action.notnull()]
noAction.loc[:,'targetAction'] = noAction.target
target_concat_df = pd.concat([withAction, noAction], ignore_index=True)[['subst_id', 'targetAction']]
target_df = target_concat_df.groupby('subst_id').agg(lambda x : '\n'.join(x)).reset_index()
target_df.columns = ['subst_id', 'targetActionList']
target_df = pd.merge(target_df, target_concat_df, how='left', on='subst_id', 
                     left_index=False, right_index=False, sort=False)
compound_df = pd.merge(compound_df, target_df, how='left', on='subst_id', 
                     left_index=False, right_index=False, sort=False)
compound_df.drop(['target', 'action', 'pharmacological_action'], axis=1, inplace=True)

study_df = pd.merge(study_df, range_df,
                    how='left', on='study_id', left_index=False, right_index=False,
                    sort=False)
study_cmpd_df = pd.merge(study_df[['study_id', 'subst_id', 'normalised_administration_route',
                            'normalised_species', 'normalised_strain', 'source_company',
                            'exposure_period_days', 'report_number']],
                         compound_df[['subst_id', 'smiles', 'status', 'common_name',
                            'cas_number','targetActionList', 'targetAction']],
                         how='left', on='subst_id', left_index=False, right_index=False,
                         sort=False)
all_df = pd.merge(study_cmpd_df,
                  findings_df[['study_id','source', 'observation', 'parameter', 'dose',
                               'relevance', 'sex']],
                  how='left', on='study_id', left_index=False, right_index=False,
                  sort=False)

# Get the summary of data per single substance
study_count_df = all_df[['subst_id', 'normalised_species', 'study_id']].groupby(['subst_id', 'normalised_species']).study_id.nunique().reset_index()
study_count_df.columns = ['subst_id', 'normalised_species', 'study_count']
study_count_df.loc[:,'count'] = study_count_df.normalised_species + ': ' + study_count_df.study_count.astype(int).astype(str)
study_count_df = study_count_df[['subst_id', 'count']].groupby('subst_id').agg(lambda x : '\n'.join(x)).reset_index()

# Generate output dataframe
subst_info_df = pd.merge(target_df, study_count_df, 
                        how='left', on='subst_id', left_index=False, right_index=False, 
                        sort=False)
subst_info_df = subst_info_df.drop_duplicates()

@api_view(['GET'])
def initFindings(request):
    global all_df, compound_df, subst_info_df,findings_df,study_df,output_df

    num_studies = len(all_df.study_id.unique().tolist())
    num_structures = len(all_df.subst_id.unique().tolist())
    num_findings = findings_df.groupby(['dose', 'observation', 'parameter', 'relevance', 'sex', 'source', 'study_id']).ngroups

    optionsDict = {}

    optionsDict['sources'] = all_df.source.dropna().unique().tolist()
    optionsDict['sources'].sort()

    optionsDict['routes'] = all_df.normalised_administration_route.dropna().unique().tolist()
    optionsDict['routes'].remove('Unknown')
    optionsDict['routes'].remove('Unassigned')
    optionsDict['routes'].sort()

    optionsDict['sex'] = all_df.sex.dropna().unique().tolist()
    optionsDict['sex'].sort()

    optionsDict['species'] = all_df[all_df.normalised_species != 'Excluded term'].normalised_species.dropna().unique().tolist()
    optionsDict['species'].sort()

    optionsDict['pharmacological_action'] = target_df.targetAction.dropna().unique().tolist()
    optionsDict['pharmacological_action'].sort()

    optionsDict['compound_name'] = all_df.common_name.dropna().unique().tolist()
    optionsDict['compound_name'].sort()

    optionsDict['cas_number'] = all_df.cas_number.dropna().unique().tolist()
    optionsDict['cas_number'].sort()

    exposure_range = all_df.exposure_period_days.dropna().unique().tolist()
    exposure_range.sort()
    optionsDict['exposure_min'] = int(exposure_range[0])
    optionsDict['exposure_max'] = int(exposure_range[-1])

    optionsDict['parameters'] = {}
    optionsDict['observations'] = {}
    for source in optionsDict['sources']:
        organs = all_df[all_df.source.str.lower() == source.lower()].parameter.dropna().unique().tolist()
        optionsDict['parameters'][source] = organs

        observations = all_df[all_df.source.str.lower() == source.lower()].observation.dropna().unique().tolist()
        optionsDict['observations'][source] = observations

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
        num_pages = math.ceil(total / 5)

    # Range of pages to show
    previous = 0
    nexts = min([7, num_pages])
    range_pages = range(1, num_pages + 1)[previous:nexts]
    previous_page = page - 1
    next_page = page + 1

    if (next_page >= num_pages):
        next_page = 0

    output_df = pd.merge(compound_df[['subst_id', 'cas_number','common_name', 'smiles',
                                        'targetActionList']], 
                         study_count_df, 
                         how='left', on='subst_id', left_index=False, right_index=False, 
                         sort=False)

    plot_info = {}

    # Species
    normalised_species = study_df.groupby(['normalised_species'])['normalised_species'].count()
    normalised_species.sort_values(ascending=False, inplace=True)
    plot_info['normalised_species'] = [[],[]]
    sum_value = 0
    for index, value in normalised_species.iteritems():
        p = float(value)/num_studies
        if p < 0.015:
            plot_info['normalised_species'][0].append('Other')
            plot_info['normalised_species'][1].append(num_studies-sum_value)
            break
        sum_value += value
        plot_info['normalised_species'][0].append(index)
        plot_info['normalised_species'][1].append(value)

    # Relevance
    relevance = findings_df.groupby(['relevance'])['relevance'].count()
    relevance.sort_values(ascending=False, inplace=True)
    plot_info['relevance'] = [relevance.index, relevance.values]

    # Source
    source = findings_df.groupby(['source'])['source'].count()
    source.sort_values(ascending=False, inplace=True)
    plot_info['source'] = [[], []]
    sum_value = 0
    for index, value in source.iteritems():
        p = float(value)/num_findings
        if p < 0.015:
            plot_info['source'][0].append('Other')
            plot_info['source'][1].append(num_findings-sum_value)
            break
        sum_value += value
        plot_info['source'][0].append(index)
        plot_info['source'][1].append(value)

    output_df = output_df.drop_duplicates()
    output_df.common_name = output_df.common_name.str.replace(', ', '\n')

    results = {
        'data': output_df[init:end].fillna(value="-").to_dict('records'),
        'allOptions': optionsDict,
        'plotInfo': plot_info,
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page,
        'num_studies': num_studies,
        'num_structures': num_structures,
        'num_findings': num_findings
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def findings(request):

    global all_df, study_df, output_df

    #####################
    # Apply all filters #
    #####################

    filtered_tmp = all_df[:]

    # Relevancy
    relevant = request.GET.get("treatmentRelated")
    if relevant:
        filtered_tmp = filtered_tmp[filtered_tmp.relevance == 'Treatment related']

    # Sex
    sex = request.GET.get("sex")
    if sex:
        filtered_tmp = filtered_tmp[filtered_tmp.sex == sex]

    # Exposure
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    if min_exposure and max_exposure:
        filtered_tmp = filtered_tmp[(filtered_tmp.exposure_period_days >= int(min_exposure)) &
                    (filtered_tmp.exposure_period_days <= int(max_exposure))]

    ##
    ## Filter organs, observations and grades by category
    ##

    # Create a copy of the filtered df without the parameter / observation
    # filtering aplied, to use when populating the optionsDict for 
    # these two parameters
    category_filtered_tmp = filtered_tmp[:]
    # Use this df to store each parameter / observation and add them together
    # so they don't become mutually exclusive
    additive_df = pd.DataFrame(columns=filtered_tmp.columns)
    
    # Organs
    all_organs = request.GET.getlist("parameters")
    if len(all_organs) > 0:
        all_organs = all_organs[0].split('@')
        tmp_dict = {}
        for v in all_organs:
            category, val = v.split('|')
            # Expand based on the ontology
            expanded_val = list(organ_onto_df[organ_onto_df.parent_term == val].child_term)
            expanded_val = [val]+expanded_val
            if category not in tmp_dict:
                tmp_dict[category] = expanded_val
            else:
                tmp_dict[category].extend(expanded_val)
        for category in tmp_dict:
            or_df = filtered_tmp[(filtered_tmp.source == category.strip()) &
                                 (filtered_tmp.parameter.isin(tmp_dict[category]))]
            additive_df = pd.concat([additive_df, or_df])

    # Observations
    all_observations = request.GET.getlist("observations")
    if len(all_observations) > 0:
        all_observations = all_observations[0].split('@')
        tmp_dict = {}
        for v in all_observations:
            category, val = v.split('|')
            category = category.strip()
            val = val.strip()
            # Expand based on the ontology
            expanded_val = list(observation_onto_df[observation_onto_df.parent_term == val].child_term)
            expanded_val = [val]+expanded_val
            if category not in tmp_dict:
                tmp_dict[category] = expanded_val
            else:
                tmp_dict[category].extend(expanded_val)
            if category not in tmp_dict:
                tmp_dict[category] = [val]
            else:
                tmp_dict[category].append(val)
        for category in tmp_dict:
            or_df = filtered_tmp[(filtered_tmp.source == category) & 
                                 (filtered_tmp.observation.isin(tmp_dict[category]))]
            additive_df = pd.concat([additive_df, or_df])
    
    if not additive_df.empty:
        filtered_tmp = additive_df[:]

    queryDict = {}
    filtered = filtered_tmp[:]
    # Pharmacological action
    all_pharm = request.GET.getlist("pharmacological_action")
    if len(all_pharm) > 0:
        queryDict['pharmacological_action'] = 'targetAction == @all_pharm'
        filtered.query('targetAction == @all_pharm', inplace=True)
        category_filtered_tmp.query('targetAction == @all_pharm', inplace=True)

    # Compound name
    all_compound_name = request.GET.getlist("compound_name")
    if len(all_compound_name) > 0:
        # Solve issue with plus signs in compound names being converted to spaces
        plus_signs = [x.replace(' ', '+') for x in all_compound_name]
        all_compound_name = all_compound_name+plus_signs
        all_compound_name = list(set(all_compound_name))
        queryDict['compound_name'] = 'common_name == @all_compound_name'
        filtered.query('common_name == @all_compound_name', inplace=True)
        category_filtered_tmp.query('common_name == @all_compound_name', inplace=True)

    # CAS number
    all_cas_number = request.GET.getlist("cas_number")
    if len(all_cas_number) > 0:
        queryDict['cas_number'] = 'cas_number == @all_cas_number'
        filtered.query('cas_number == @all_cas_number', inplace=True)
        category_filtered_tmp.query('cas_number == @all_cas_number', inplace=True)

    # Administration route
    all_routes = request.GET.getlist("routes")
    if len(all_routes) > 0:
        queryDict['routes'] = 'normalised_administration_route == @all_routes'
        filtered.query('normalised_administration_route == @all_routes', inplace=True)
        category_filtered_tmp.query('normalised_administration_route == @all_routes', inplace=True)

    # Species
    all_species = request.GET.getlist("species")
    if len(all_species) > 0:
        queryDict['species'] = 'normalised_species == @all_species'
        filtered.query('normalised_species == @all_species', inplace=True)
        category_filtered_tmp.query('normalised_species == @all_species', inplace=True)

    ####################################
    # Generate optionsDict by applying #
    # all filters but the ones in      #
    # that category                    #
    ####################################

    optionsDict = {}
    categories = all_df.source.dropna().unique().tolist()
    if not filtered.empty:
        optionsDict['parameters'] = {}
        optionsDict['observations'] = {}
        valuesL = list(queryDict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = category_filtered_tmp.query(query_string)
        else:
            tmp_df = category_filtered_tmp[:]

        for category in categories:
            organs = tmp_df[tmp_df.source.str.lower() == category.lower()].parameter.dropna().unique().tolist()
            optionsDict['parameters'][category] = organs
            optionsDict['parameters'][category].sort()

            observations = tmp_df[tmp_df.source.str.lower() == category.lower()].observation.dropna().unique().tolist()
            optionsDict['observations'][category] = observations
            optionsDict['observations'][category].sort()
            
        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('pharmacological_action', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['pharmacological_action'] = tmp_df.targetAction.dropna().unique().tolist()
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
        try:
            optionsDict['routes'].remove('Unknown')
            optionsDict['routes'].remove('Unassigned')
        except:
            pass
        optionsDict['routes'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('species', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['species'] = tmp_df[tmp_df.normalised_species != 'Excluded term'].normalised_species.dropna().unique().tolist()
        optionsDict['species'].sort()

        exposure_range = filtered.exposure_period_days.dropna().unique().tolist()
        exposure_range.sort()
        optionsDict['exposure_min'] = int(exposure_range[0])
        optionsDict['exposure_max'] = int(exposure_range[-1])

        optionsDict['sex'] = all_df.sex.dropna().unique().tolist()
        optionsDict['sex'].sort()

    #############
    # Aggregate #
    #############

    num_studies = filtered.study_id.nunique()
    num_structures = filtered.subst_id.nunique()

    filtered_findings = filtered[['dose', 'observation', 'parameter', 'relevance', 'sex', 'source', 'study_id']].drop_duplicates().groupby(['dose', 'observation', 'parameter', 'relevance', 'sex', 'source', 'study_id'])
    num_findings = filtered_findings.ngroups
    filtered_findings = filtered_findings.count().reset_index()

    ##PLOT INFO
    plot_info = {}
    # Species
    normalised_species = filtered.groupby(['normalised_species','study_id']).count()
    normalised_species = normalised_species.reset_index().groupby(['normalised_species'])['normalised_species'].count()
    normalised_species.sort_values(ascending=False, inplace=True)
    plot_info['normalised_species'] = [[],[]]
    sum_value = 0
    for index, value in normalised_species.iteritems():
        p = float(value)/num_studies
        if p < 0.015:
            plot_info['normalised_species'][0].append('Other')
            plot_info['normalised_species'][1].append(num_studies-sum_value)
            break
        sum_value += value
        plot_info['normalised_species'][0].append(index)
        plot_info['normalised_species'][1].append(value)

    # Treatment related
    relevance = filtered_findings.groupby(['relevance'])['relevance'].count()
    relevance.sort_values(ascending=False, inplace=True)
    plot_info['relevance'] = [relevance.index, relevance.values]

    # Source
    source = filtered_findings.groupby(['source'])['source'].count()
    source.sort_values(ascending=False, inplace=True)
    plot_info['source'] = [[], []]
    sum_value = 0
    for index, value in source.iteritems():
        p = float(value)/num_findings
        if p < 0.015:
            plot_info['source'][0].append('Other')
            plot_info['source'][1].append(num_findings-sum_value)
            break
        sum_value += value
        plot_info['source'][0].append(index)
        plot_info['source'][1].append(value)

    if not filtered.empty:
        study_count_df = filtered.dropna(subset=['normalised_species'])[['subst_id', 'normalised_species', 'study_id']].groupby(['subst_id', 'normalised_species']).study_id.nunique().reset_index()
        study_count_df.columns = ['subst_id', 'normalised_species', 'study_count']
        study_count_df.loc[:,'count'] = study_count_df.normalised_species + ': ' + study_count_df.study_count.astype(int).astype(str)
        study_count_df = study_count_df[['subst_id', 'count']].groupby('subst_id').agg(lambda x : '\n'.join(x)).reset_index()

        output_df = pd.merge(filtered[['subst_id']].drop_duplicates(),
                        compound_df[['subst_id', 'cas_number', 'common_name', 'company_id', 
                                'smiles', 'status', 'targetActionList']].drop_duplicates(), 
                        how='left', on='subst_id', left_index=False, right_index=False, 
                        sort=False)
        output_df = pd.merge(output_df, study_count_df,
                        how='left', on='subst_id', left_index=False, right_index=False, 
                        sort=False)
        output_df = output_df.drop_duplicates()
        output_df.common_name = output_df.common_name.str.replace(', ', '\n')
    else:
        output_df = pd.DataFrame(columns=['subst_id', 'cas_number', 'common_name', 
                                    'company_id', 'smiles', 'status', 'targetActionList', 
                                    'count'])

    ##############
    # Pagination #
    ##############

    page = int(request.GET.get("page"))

    if page != 0:
        init = (page - 1) * 5
        end = init + 5
    else:
        init = 0
        end = len(filtered)

    num_pages = math.ceil(len(output_df) / 5)

    # Range of pages to show
    if page < 4:
        previous = 0
        nexts = min([7, num_pages])
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

    results = {
        'data': output_df[init:end].fillna(value="-").to_dict('records'),
        'allOptions': optionsDict,
        'plotInfo': plot_info,
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page,
        'num_studies': num_studies,
        'num_structures': num_structures,
        'num_findings' : num_findings
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def page(request):

    global output_df

    page = int(request.GET.get("page"))

    if page != 0:
        init = (page - 1) * 5
        end = init + 5
    else:
        init = 0
        end = len(output_df)

    num_pages = math.ceil(len(output_df) / 5)

    # Range of pages to show
    if page < 4:
        previous = 0
        nexts = min([7, num_pages])
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

    results = {
        'data': output_df[init:end].fillna(value="-").to_dict('records'),
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page
    }

    send_data = Pageserializer(results, many=False).data
    return Response(send_data)
