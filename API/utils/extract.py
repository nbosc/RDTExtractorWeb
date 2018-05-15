##    Description    eTOX repeat-dose toxicity extraction tool
##
##    Authors:       Elisabet Gregori (elisabet.gregori@upf.edu)
##                   Ignacio Pasamontes (ignacio.pasamontes@upf.edu)
##
##    Copyright 2018 Elisabet Gregori & Ignacio Pasamontes
##
##    RDTextractor is free software: you can redistribute it 
##    and/or modify it under the terms of the GNU General Public 
##    License as published by the Free Software Foundation version 3.
##
##    RDTextractor is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

import pandas as pd


def filter_study(dict,df):

    # Exposure
    if  'min_exposure' in dict  and 'max_exposure' in dict:
        # An exposure range filter is defined
        df = df[(df.exposure_period_days >= dict['min_exposure']) &
                (df.exposure_period_days <= dict['max_exposure'])]
    elif 'min_exposure' in dict:
        # Only a.upper bound for exposure range has been set
        df = df[df.exposure_period_days >= dict['min_exposure']]
    elif 'max_exposure' in dict:
        df = df[df.exposure_period_days <= dict['max_exposure']]

    # Administration route
    if 'route' in dict:
        df = df[df.normalised_administration_route.str.lower().isin([x.lower() for x in dict['route']])]
        
    # Species
    if 'species' in dict:
        df = df[df.normalised_species.str.lower().isin([x.lower() for x in dict['species']])]
    
    return df

def expand(df,dict,onto_df):

    """
    Expand standardized observation and normalized organs based on
    the hierarchy of ontologies stored in onto_df
    """
    # Create an empty output dataframe
       
    #########
    # ORGAN #
    #########
    if 'organs' in dict:

        organs_dict = {}
        all_organs = set()
        onto_anatomy = onto_df[(onto_df.ontology == 'anatomy')]

        # get for each organ get all its childs and crate a dict
        for organ in dict['organs']:
            all_organs.union([organ])
            related_organs = onto_anatomy[(onto_anatomy.parent_term == organ)]
            all_organs = all_organs.union(related_organs.child_term)
            organs_dict.update({n: organ for n in related_organs.child_term})

        df = df[df['organ_normalised'].isin(all_organs)]
        df.loc[:, 'organ_normalised'] = df['organ_normalised'].map(organs_dict)


    ###############
    # OBSERVATION #
    ###############
    if 'observations' not in dict:

        onto_histopatho = onto_df[(onto_df['ontology'] == 'histopathology')
                                  & (onto_df['parent_term'] != "morphologic change")]

        findings_out = pd.merge(df, onto_histopatho, how='left',
                                left_on='observation_normalised',
                                right_on='child_term')

        findings_out = findings_out[['study_id', 'relevance', 'parent_term',
                                     'organ_normalised', 'dose', 'subst_id', 'report_number']]

        findings_out=findings_out.rename(index=str, columns={"parent_term": "observation_normalised"})
        findings_out.drop_duplicates(inplace=True)

    else:

        observation_dict = {}
        all_observations = set()
        onto_histopatho = onto_df[(onto_df.ontology == 'histopathology')]

        # get for each observation all its childs and crate a dict
        for observation in dict['observations']:

            all_observations.union([observation])
            related_observation = onto_histopatho[(onto_histopatho.parent_term == observation)]
            all_observations = all_observations.union(related_observation.child_term)
            observation_dict.update({n: observation for n in related_observation.child_term})

        findings_out = df[df['observation_normalised'].isin(all_observations)]
        findings_out.loc[:, 'observation_normalised'] = df['observation_normalised'].map(observation_dict)


    return findings_out

def get_stats(group):
    return {'min': group.min(), 'max': group.max()}

def run(args):

    """
    Run the data extraction based on the parsed filters and expanding
    based on the organs and morphological changes ontologies.
    """

    sys.stderr.write('\nLoading background information for version %s\n' %args.version)
    study_df, find_df = load_version(args)
    
    #################################
    # Select only relevant findings #
    #################################
    sys.stderr.write('Filtering to relevant information\n')
    relevant_studies_df = filter_study(args,study_df)
    relevant_find = find_df[find_df.study_id.isin(relevant_studies_df.study_id)]
    relevant_find = pd.merge(relevant_find, study_df[['study_id', 'subst_id']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)
    
    ###################################
    # Get stats for relevant findings #
    ###################################
    # Get the number of studies per substance
    count_df = relevant_find.groupby(('subst_id')).study_id.nunique().to_frame().reset_index()
    # Get the global dose range per substance
    range_df = relevant_find[relevant_find.dose > 0]
    range_df = range_df.groupby(('subst_id')).dose.apply(get_stats).unstack().reset_index()
    # Get all stats into a single dataframe
    stats_df = pd.merge(count_df, range_df, how='inner', on='subst_id', 
                        left_index=False, right_index=False, sort=False)
    stats_df.columns = ['subst_id', 'study_count', 'dose_max', 'dose_min']

    ###################################################################
    # Expand based on anatomical and morphological changes ontologies #
    ###################################################################
    # Expand organs and histopathological findings according to the ontologies 
    # and filter by finding-based arguments
    sys.stderr.write('Expand based on anatomic and morphological change ontologies\n')
    filtered_find = expand(relevant_find,args)

    if filtered_find.empty:
        raise Exception('Filtered out all rows, so the dataframe is empty.')

    ######################################
    # Aggragate by substance and finding #
    ######################################
    # Define finding as organ+observation
    filtered_find.organ_normalised = filtered_find.organ_normalised.fillna('NA')
    filtered_find.observation_normalised = filtered_find.observation_normalised.fillna('NA')
    filtered_find['finding'] = filtered_find.apply(lambda row: row.organ_normalised+'_'+row.observation_normalised, axis=1)
    filtered_find = filtered_find[['subst_id', 'finding', 'dose']]

    # Aggregate by substance and finding (as defined above), keeping the minimum dose 
    # for each substance/finding instance
    group_df = filtered_find.groupby(('subst_id', 'finding')).min().add_prefix('min_').reset_index()
    
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
    cols = cols[0:4]+[cols[-1]]+cols[4:-1]
    quantitative_df = quantitative_df[cols]
    
    ### Qualitative
    group_df.loc[:,'min_dose'] = 1
    pivotted_df = group_df.pivot_table(index='subst_id', columns='finding', values='min_dose').reset_index()
    pivotted_df['is_active'] = 'True'
    qualitative_df = pd.merge(stats_df, pivotted_df, how='left', on='subst_id',
                                left_index=False, right_index=False, sort=False)
    qualitative_df.is_active = qualitative_df.is_active.fillna('False')
    # Reorder columns
    cols = qualitative_df.columns.tolist()
    cols = cols[0:4]+[cols[-1]]+cols[4:-1]
    qualitative_df = qualitative_df[cols]

    ####################
    # Save the results #
    ####################
    quantitative_df.to_csv(args.output_basename+'_quant.tsv', 
                            sep='\t', index=False)
    qualitative_df.to_csv(args.output_basename+'_qual.tsv', 
                            sep='\t', index=False)
