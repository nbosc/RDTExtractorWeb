# from snippets.models import Snippet
# from snippets.serializers import SnippetSerializer
from django.http import Http404
from rest_framework.response import Response
from rest_framework.decorators import api_view
from rest_framework import status
import pandas as pd
# Disable SettingWithCopyWarning warnings
pd.set_option('chained_assignment', None)
import json
import copy
import math
import pickle
from django.views.decorators.csrf import csrf_exempt
from .serializers import FindingSerializer,Pageserializer
from API.utils import extract

def get_stats(group):
    return {'min': group.min(), 'max': group.max()}

# Load dataframes with information on studies, compounds and findings
substance_df = pd.read_pickle("API/static/data/substance.pkl")
all_df = pd.read_pickle("API/static/data/all.pkl.gz", compression='gzip')
organ_onto_df = pd.read_pickle("API/static/data/organ_ontology.pkl")
observation_onto_df = pd.read_pickle("API/static/data/observation_ontology.pkl")

@api_view(['GET'])
def initFindings(request):

    global output_df, optionsDict

    output_df = pd.read_pickle("API/static/data/output.pkl")
    results = pickle.load(open("API/static/data/init_results.pkl", 'rb'))
    optionsDict = results['allOptions']
    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def findings(request):

    global all_df, study_df, output_df, optionsDict

    #####################
    # Apply all filters #
    #####################

    filtered = all_df[:]

    # Relevancy
    relevant = request.GET.get("treatmentRelated")
    if relevant:
        filtered = filtered[filtered.relevance == 'Treatment related']

    # Sex
    sex = request.GET.get("sex")
    if sex:
        filtered = filtered[filtered.sex == sex]

    # Exposure
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    if min_exposure and max_exposure:
        filtered = filtered[(filtered.exposure_period_days >= int(min_exposure)) &
                    (filtered.exposure_period_days <= int(max_exposure))]

    ##
    ## Filter parameters, observations and grades by category
    ##
    # Use this df to store each parameter / observation and add them together
    # so they don't become mutually exclusive
    additive_df = pd.DataFrame(columns=filtered.columns)
    
    # parameters and observations
    all_parameters = request.GET.getlist("parameters")
    all_categories = set([])
    if len(all_parameters) > 0:
        all_parameters = all_parameters[0].split('@')
        tmp_parameters_dict = {}
        for v in all_parameters:
            category, val = v.split('|')
            # Expand based on the ontology
            expanded_val = list(organ_onto_df[organ_onto_df.parent_term == val].child_term)
            expanded_val = [val]+expanded_val
            if category not in tmp_parameters_dict:
                tmp_parameters_dict[category] = expanded_val
            else:
                tmp_parameters_dict[category].extend(expanded_val)
        all_categories = set(tmp_parameters_dict.keys())

    # Observations
    all_observations = request.GET.getlist("observations")
    if len(all_observations) > 0:
        all_observations = all_observations[0].split('@')
        tmp_observations_dict = {}
        for v in all_observations:
            category, val = v.split('|')
            category = category.strip()
            val = val.strip()
            # Expand based on the ontology
            expanded_val = list(observation_onto_df[observation_onto_df.parent_term == val].child_term)
            expanded_val = [val]+expanded_val
            if category not in tmp_observations_dict:
                tmp_observations_dict[category] = expanded_val
            else:
                tmp_observations_dict[category].extend(expanded_val)
            if category not in tmp_observations_dict:
                tmp_observations_dict[category] = [val]
            else:
                tmp_observations_dict[category].append(val)
        all_categories = all_categories.union(set(tmp_parameters_dict.keys()))

    if bool(all_categories):
        # At least one observation or parameter filter has been applied
        for category in all_categories:
            if category in tmp_parameters_dict and category in tmp_observations_dict:
                or_df = filtered[(filtered.source == category.strip()) &
                                (filtered.parameter.isin(tmp_parameters_dict[category])) &
                                (filtered.observation.isin(tmp_observations_dict[category]))]
            elif category in tmp_parameters_dict:
                or_df = filtered[(filtered.source == category.strip()) &
                                (filtered.parameter.isin(tmp_parameters_dict[category]))]
            elif category in tmp_observations_dict:
                or_df = filtered[(filtered.source == category) & 
                                (filtered.observation.isin(tmp_observations_dict[category]))]
            additive_df = pd.concat([additive_df, or_df])
        filtered = additive_df[:]

    # Pharmacological action
    all_pharm = request.GET.getlist("pharmacological_action")
    if len(all_pharm) > 0:
        filtered = filtered[filtered.targetAction.isin(all_pharm)]

    # Compound name
    all_compound_name = request.GET.getlist("compound_name")
    if len(all_compound_name) > 0:
        # Solve issue with plus signs in compound names being converted to spaces
        plus_signs = [x.replace(' ', '+') for x in all_compound_name]
        all_compound_name = all_compound_name+plus_signs
        all_compound_name = list(set(all_compound_name))
        filtered = filtered[filtered.common_name.isin(all_compound_name)]

    # CAS number
    all_cas_number = request.GET.getlist("cas_number")
    if len(all_cas_number) > 0:
        filtered = filtered[filtered.cas_number.isin(all_cas_number)]

    # Administration route
    all_routes = request.GET.getlist("routes")
    if len(all_routes) > 0:
        filtered = filtered[filtered.normalised_administration_route.isin(all_routes)]

    # Species
    all_species = request.GET.getlist("species")
    if len(all_species) > 0:
        filtered = filtered[filtered.normalised_species.isin(all_species)]

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
                        substance_df[['subst_id', 'cas_number', 'common_name', 'company_id', 
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

    num_pages = math.ceil(output_df.subst_id.nunique() / 5.)

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

    if (next_page > num_pages):
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

    total = output_df.subst_id.nunique()
    if page != 0:
        init = (page - 1) * 5
        end = init + 5
    else:
        init = 0
        end = total

    num_pages = math.ceil(total / 5.)

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

    if (next_page > num_pages):
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
