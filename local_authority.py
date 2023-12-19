'''Investigate the link between local authority deprivation and voting in the U

Local authority boundary data:
    Source: Office for National Statistics licensed under the Open Government Licence v.3.0
    Contains OS data © Crown copyright and database right 2023

Local authority deprivation data:
    2019 data
    Source: Office for National Statistics licensed under the Open Government Licence v.3.0

Election boundary data
    Source: Office for National Statistics licensed under the Open Government Licence v.3.0
    Contains OS data © Crown copyright and database right 2023

Election results data:
    2019 data
    Source: Office for National Statistics licensed under the Open Government Licence v.3.0
'''

import os
import logging
from pathlib import Path
import zipfile
import pandas as pd
import geopandas
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
import requests

logging.basicConfig(level=logging.INFO)

class LocalAuthorityData:
    '''Contains the data pertaining to Local Authorities, such as size and 
    location, and deprivation rates'''

    def __init__(self):
        self.local_authority_shapes_raw = """raw/local_auth/Local_Authority_
                                           Districts__December_2019__Boundaries
                                           _UK_BFC.shp"""
        self.local_authority_deprivation_raw = """raw/local_auth/
                                                localincomedeprivationdata.xlsx"""
        self.local_authority_shapes_gdf = geopandas.GeoDataFrame()
        self.local_authority_deprivation_df = pd.DataFrame()
        self.local_authority_gdf = geopandas.GeoDataFrame()
        self.output_file = 'combined/local_authority_deprivation'
        logging.info('LAD: Vars init for LocalAuthorityData')
        self.initialise_data()


    def initialise_data(self):
        '''Check to see if the dataframes have already been combined, if
        not create the relevant folders and combine them'''

        logging.info('LAD: Initialise_data function begin')
        # Check to see if the file exists, if it is, load it in as
        # self.local_authority_gdf. If not, run through the download, and load
        # functions
        if Path(self.output_file).is_file():
            logging.info(f'LAD: Combined dataset as {self.output_file} already exists')
            logging.info('LAD: Merged database file exists, loading in now')
            self.local_authority_gdf = geopandas.read_file(self.output_file)
            logging.info('LAD: Merged database loaded')
        else:
            logging.info('LAD: Merged data does not exists, downloading')

            # Create folders for separating data if they don't exist
            if not os.path.exists('raw'):
                os.makedirs('raw')
            if not os.path.exists('combined'):
                os.makedirs('combined')
            if not os.path.exists('raw/local_auth'):
                os.makedirs('raw/local_auth')
            logging.info('LAD: Folders created')

            self.download_files()  # Downloads the shape and data files if needed
            self.load_files()  # Loads the files into the program
            self.combine_dataframes()  # Creates a combined file of shape + data

    def download_files(self):
        '''Download the local authority boundaries and deprication data 
        if it's not downloaded already. Stores in the 'raw' directory
        If already downloaded, does not download anything'''

        deprivation_url = """https://www.ons.gov.uk/file?uri=/peoplepopulation
                            andcommunity/personalandhouseholdfinances/income
                            andwealth/datasets/mappingincomedeprivationatalocal
                            authoritylevel/2019/localincomedeprivationdata.xlsx"""

        authority_boundaries_url = """https://stg-arcgisazurecdataprod1.az
                                    .arcgis.com/exportfiles-1559-17928/
                                    LAD_Dec_2019_Boundaries_UK_BFC_2022_
                                    2942927368901013076.zip?sv=2018-03-28
                                    &sr=b&sig=6lfIyORiFMeWFv0NGIX645Kqsk
                                    YRvd210sr4jOQzXcg%3D&se=2023-09-09T14%
                                    3A55%3A37Z&sp=r"""

        logging.info('LAD: Checking if deprivation data exists')
        # Check to see if the deprivation data has already been downloaded
        if Path(self.local_authority_deprivation_raw).is_file():
            logging.info(f"""LAD: File {self.local_authority_deprivation_raw}
                          already exists""")
        else:
            # Try to download the deprivation data
            logging.info('LAD: Deprivation data needs downloading')
            response_depriv = requests.get(deprivation_url, stream=True,
                                           timeout=5000)
            if response_depriv.status_code == 200:
                # Write the download data to the file specified
                with open(self.local_authority_deprivation_raw, 'wb') as file:
                    for chunk in response_depriv.iter_content(chunk_size=1024):
                        if chunk:
                            file.write(chunk)
                logging.info(f"""LAD: File {self.local_authority_deprivation_raw}
                             downloaded successfully.""")
            else:
                logging.info(f"""LAD: Failed to download the file.
                             Status code: {response_depriv.status_code}""")

        logging.info('LAD: Checking is local authority shape data exists')
        # Check to see if the authority boundary data has already been downloaded
        if Path(self.local_authority_shapes_raw).is_file():
            logging.info(f'LAD: File {self.local_authority_shapes_raw} already exists')
        else:
            # Try to download the authority data as a zip archive
            logging.info('LAD: Downloading local authority shape data')
            local_auth_shape_zip = 'raw/local_auth/local_auth_shape_zip.zip'
            response_auth = requests.get(authority_boundaries_url, stream=True,
                                         timeout=5000)
            if response_auth.status_code == 200:
                # Write the download data to the file specified
                with open(local_auth_shape_zip, 'wb') as file:
                    for chunk in response_auth.iter_content(chunk_size=1024):
                        if chunk:
                            file.write(chunk)
                logging.info(f"""LAD: File {self.local_authority_shapes_raw}
                              zip archive downloaded successfully.""")
            else:
                logging.info(f"""LAD: Failed to download the file.
                              Status code: {response_auth.status_code}""")

            # Destination for the unzipped files is in the raw folder
            destination_directory = 'raw/local_auth'

            logging.info('LAD: Unzipping the shape date files')
            # Unzip the file
            with zipfile.ZipFile(local_auth_shape_zip, 'r') as zip_ref:
                zip_ref.extractall(destination_directory)
        logging.info('LAD: File preparation complete, files in raw directory')

    def load_files(self):
        '''Loads the shape file and the deprivation data into separate dataframes
        Set the index as the Local Authority Code (varies by year, so hardcoded
        for 2019)'''

        # Load in the raw shape data as a geopandas dataframe
        logging.info('LAD: Begin loading files')
        self.local_authority_shapes_gdf = geopandas.read_file(self.local_authority_shapes_raw)
        self.local_authority_shapes_gdf = self.local_authority_shapes_gdf.set_index('lad19cd')

        # Load in the correct sheet from the excel document with all the data in
        # Load in as a pandas df
        logging.info('LAD: Loaded shape gdf, loading deprivation df')
        self.local_authority_deprivation_df = pd.read_excel(
                                                self.local_authority_deprivation_raw,
                                                'Rankings for all indicators',
                                                header=1)
        self.local_authority_deprivation_df = (self.local_authority_deprivation_df
                                               .set_index
                                               ('Local Authority District code (2019)'))

    def combine_dataframes(self):
        '''Combine the shape and deprivation dataframes into one dataframe
        Keeps all shape data, and adds deprivation data to the correct
        boundary.
        Saves the combined geopandas dataframe as a GeoJSON file for quicker
        loading next time the program is run'''

        # Merge the data dataframe into the geopandas dataframe using the
        # local authority code as the key to match up the data
        # Throws away any deprivation data that doesn't have a matching shape
        logging.info('LAD: Merging dataframes')
        self.local_authority_gdf = (self.local_authority_shapes_gdf
                                    .merge(self.local_authority_deprivation_df,
                                           left_on='lad19cd',
                                           right_on='Local Authority District code (2019)',
                                           how='left'))

        # Fill any missing data with 0, so it will still be plotted
        self.local_authority_gdf = self.local_authority_gdf.fillna(0)
        logging.info('LAD: Dataframes merged')

        # Save the file as a GeoJSON
        self.local_authority_gdf.to_file(self.output_file, driver='GeoJSON')


class ElectoralResults:
    '''Contains data related to electoral boundaries, and the voting
    results of the 2019 UK General Election'''

    def __init__(self):
        self.electoral_shapes_raw = 'raw/electoral/PCON_DEC_2019_UK_BFC.shp'
        self.electoral_results_raw = 'raw/electoral/HoC-GE2019-results-by-constituency-xlsx.xlsx'
        self.electoral_shapes_gdf = geopandas.GeoDataFrame()
        self.electoral_results_df = pd.DataFrame()
        self.electoral_gdf = geopandas.GeoDataFrame()
        self.output_file = 'combined/electoral_results'
        self.party_colour_dict = {}
        logging.info('ER: Vars init for ElectoralResults')
        self.initialise_data()


    def initialise_data(self):
        '''Check to see if the dataframes have already been combined,
        if not create the relevant folders and combine them'''

        logging.info('ER: Initialise_data function begin')
        # Check to see if the file exists, if it is, load it in as
        # self.electoral_gdf. If not, run through the download, and load functions
        if Path(self.output_file).is_file():
            logging.info(f'ER: Combined dataset as {self.output_file} already exists')
            logging.info('ER: Merged database file exists, loading in now')
            self.electoral_gdf = geopandas.read_file(self.output_file)
            logging.info('ER: Merged database loaded')
        else:
            logging.info('ER: Merged data does not exists, downloading')

            # Create folders for separating data if they don't exist
            if not os.path.exists('raw'):
                os.makedirs('raw')
            if not os.path.exists('combined'):
                os.makedirs('combined')
            if not os.path.exists('raw/electoral'):
                os.makedirs('raw/electoral')
            logging.info('ER: Folders created')

            self.download_files()  # Downloads the shape and data files if needed
            self.load_files()  # Loads the files into the program
            self.combine_dataframes()  # Creates a combined file of shape + data

        logging.info('ER: Adding new data columns')
        self.party_colours()  # Create a column containing party colour hex
        self.voter_turnout()  # Calculate the voter turnout

    def download_files(self):
        '''Download the electoral boundaries and voting rresults data 
        if it's not downloaded already. Stores in the 'raw' directory
        If already downloaded, does not download anything'''

        electoral_shape_download_url = """https://stg-arcgisazurecdataprod1
                                        .az.arcgis.com/exportfiles-1559-16027/
                                        WPC_Dec_2019_Boundaries_UK_BFC_V2_2022
                                        _-7467643007870890518.zip?sv=2018-03-28
                                        &sr=b&sig=TyQpZecp00jsj0qDPBI438B5nT6B
                                        %2B1XHhmRCI3KOu0U%3D&se=2023-09-09T15%
                                        3A34%3A42Z&sp=r"""

        electoral_results_url = """https://researchbriefings.files.parliament
                                 .uk/documents/CBP-8749/HoC-GE2019-results-by-
                                 constituency-xlsx.xlsx"""

        logging.info('ER: Checking if results data exists')
        # Check to see if the results data has already been downloaded
        if Path(self.electoral_results_raw).is_file():
            logging.info(f'ER: File {self.electoral_results_raw} already exists')
        else:
            # Try to download the deprivation data
            logging.info('ER: Deprivation data needs downloading')
            response_depriv = requests.get(electoral_results_url, stream=True,
                                           timeout=5000)
            if response_depriv.status_code == 200:
                # Write the download data to the file specified
                with open(self.electoral_results_raw, 'wb') as file:
                    for chunk in response_depriv.iter_content(chunk_size=1024):
                        if chunk:
                            file.write(chunk)
                logging.info(f'ER: File {self.electoral_results_raw} downloaded successfully.')
            else:
                logging.info(f"""ER: Failed to download the file.
                              Status code: {response_depriv.status_code}""")

        logging.info('ER: Checking is electoral boundary shape data exists')
        # Check to see if the electoral boundary data has already been downloaded
        if Path(self.electoral_shapes_raw).is_file():
            logging.info(f'ER: File {self.electoral_shapes_raw} already exists')
        else:
            # Try to download the authority data as a zip archive
            logging.info('ER: Downloading local authority shape data')
            elec_shape_zip = 'raw/electoral/elec_shape_zip.zip'
            response_auth = requests.get(electoral_shape_download_url, stream=True,
                                         timeout=5000)
            if response_auth.status_code == 200:
                # Write the download data to the file specified
                with open(elec_shape_zip, 'wb') as file:
                    for chunk in response_auth.iter_content(chunk_size=1024):
                        if chunk:
                            file.write(chunk)
                logging.info(f"""ER: File {self.electoral_shapes_raw} zip
                              archive downloaded successfully.""")
            else:
                logging.info(f"""ER: Failed to download the file.
                             Status code: {response_auth.status_code}""")

            # Destination for the unzipped files is in the raw folder
            destination_directory = 'raw/electoral'

            logging.info('ER: Unzipping the shape date files')
            # Unzip the file
            with zipfile.ZipFile(elec_shape_zip, 'r') as zip_ref:
                zip_ref.extractall(destination_directory)
        logging.info('ER: File preparation complete, files in raw directory')

    def load_files(self):
        '''Loads the shape file and the results data into separate dataframes
        Set the index as the Electoral code'''

        # Load in the raw shape data as a geopandas dataframe
        logging.info('ER: Begin loading files')
        self.electoral_shapes_gdf = geopandas.read_file(self.electoral_shapes_raw)
        self.electoral_shapes_gdf = self.electoral_shapes_gdf.set_index('PCON19CD')

        # Load in the excel document of results as a pandas df
        logging.info('ER: Loaded shape gdf, loading results df')
        self.electoral_results_df = pd.read_excel(self.electoral_results_raw, header=0)
        self.electoral_results_df = self.electoral_results_df.set_index('ons_id')

    def combine_dataframes(self):
        '''Combine the shape and results dataframes into one dataframe
        Keeps all shape data, and adds results data to the correct#
        boundary
        Saves the combined geopandas dataframe as a GeoJSON file for quicker
        loading next time the program is run'''

        # Merge the data dataframe into the geopandas dataframe using the
        # electoral area as the key to match up the data
        # Throws away any deprivation data that doesn't have a matching shape
        logging.info('ER: Merging dataframes')
        self.electoral_gdf = (self.electoral_shapes_gdf
                              .merge(self.electoral_results_df,
                                     left_on='PCON19CD', right_on='ons_id',
                                     how='left'))
        self.electoral_gdf = self.electoral_gdf.fillna(0)
        logging.info('ER: Dataframes merged')
        self.electoral_gdf.to_file(self.output_file, driver='GeoJSON')

    def party_colours(self):
        '''Add a column containing the colour of the winning party in each
        constituency'''

        # Choose the colum to use for the colour and make a colour dictionary
        # that repesents each party's colour
        logging.info('ER: Adding party colour values column')
        vote_result_column = 'first_party'
        self.party_colour_dict= {'Con': '#0087DC', 'Lab': '#E4003B','LD': '#FAA61A',
                            'Green': '#02A59B', 'Spk': '#772464','DUP': '#D46A4C',
                            'SF': '#326760', 'SDLP': '#2AA82C', 'Alliance': '#772464',
                            'SNP': '#FDF38E', 'PC': '#005B54'}        

        # Create a new colum that has the representative colour for the party
        # that got into power in that seat
        self.electoral_gdf['party_colour'] = (self.electoral_gdf[vote_result_column]
                                              .map(self.party_colour_dict))

    def voter_turnout(self):
        '''Calculate the turnout in each constituency by dividing the sum of
         the valid and invalid votes  by the electorate (total number of 
         peopel who can vote)'''

        #Make a new column that contains the turnout as a decimal
        logging.info('ER: Calculating voter turnout')
        self.electoral_gdf['turnout'] = ((self.electoral_gdf['valid_votes'] +
                                         self.electoral_gdf['invalid_votes']) /
                                         self.electoral_gdf['electorate'])


class MultiplePlots:
    '''Display one plot showing deprivation, voter turnout, and GE reuslts
    Voter turnout and deprivation are overlaid and transparent to see trend'''

    def __init__(self, electoral, deprivation):
        self.electoral = electoral
        self.deprivation = deprivation
        logging.info('MP: Initialised dataframes')

    def plot(self):
        '''Create graphs showing the different data sets
        Figure 1: Show the combined turnout and deprivation and GE results'''

        # Set the coordinate systems of the two data sets to the same
        self.electoral = self.electoral.to_crs(epsg=3857)
        self.deprivation = self.deprivation.to_crs(epsg=3857)
        logging.info('MP: Set ploy coord systems to epsg 3857')


        logging.info('MP: Getting Deprivation data and colourmap')
        deprivation_data_columnm = self.depriv_plot()

        # Create a normalised column for each dataframe, to overlay easier
        elect_normalised = self.normalise_results(self.electoral, 'turnout')
        depriv_normalised = self.normalise_results(self.deprivation, deprivation_data_columnm)

        # Create a new figure wit the overlay and GE election results
        fig_comparison, ax_comparison = plt.subplots(1, 2, figsize=(12, 6))  # Create two subplots

        # Overlay the deprivation and voter turnout maps
        logging.info('MP: Plotting turnout')
        self.electoral.plot(ax=ax_comparison[0], column=self.electoral[elect_normalised],
                            cmap='Reds_r', alpha=0.5, edgecolor='black', linewidth=0.1)
        logging.info('MP: Plotting deprivation')
        self.deprivation.plot(ax=ax_comparison[0], column=depriv_normalised,
                            cmap='Blues', alpha=0.5, edgecolor='black', linewidth=0.1)

        # Create a colorbar for voter turnout
        voter_divider = make_axes_locatable(ax_comparison[0])
        voter_cax = voter_divider.append_axes("right", size="5%", pad=0.05)
        voter_sm = plt.cm.ScalarMappable(cmap='Reds_r', norm=plt.Normalize(
                                        vmin=min(self.electoral[elect_normalised]),
                                        vmax=max(self.electoral[elect_normalised])))
        voter_sm._A = []
        voter_cbar = plt.colorbar(voter_sm, cax=voter_cax)
        voter_cbar.set_label('Normalised Voter Turnout')
        logging.info('MP: Created voter turnout colourbar')

        # Create a colorbar for the deprivation
        deprivation_divider = make_axes_locatable(ax_comparison[0])
        deprivation_cax = deprivation_divider.append_axes("bottom", size="5%", pad=0.1)
        deprivation_sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(
                                        vmin=min(self.deprivation[depriv_normalised]),
                                        vmax=max(self.deprivation[depriv_normalised])))
        deprivation_sm._A = []
        deprivation_cbar = plt.colorbar(deprivation_sm, cax=deprivation_cax,
                                        orientation='horizontal')
        deprivation_cbar.set_label('Normalised Deprivation')
        logging.info('MP: Created deprivation colourbar')

        logging.info('MP: Plotting GE results data')
        self.electoral.plot(ax=ax_comparison[1], color=self.electoral['party_colour'],
                            legend=False, edgecolor='black', linewidth=0.1)

        # Get unique values from 'first_party' column and create a legend
        unique_parties = self.electoral['first_party'].unique()
        legend_labels = {party: party for party in unique_parties}
        ax_comparison[1].legend(
                        handles=[plt.Line2D([0], [0],
                            color=self.electoral['party_colour']
                            .loc[self.electoral['first_party'] == party].values[0],
                            label=legend_labels[party]) for party in unique_parties],
                            bbox_to_anchor=(1.02, 1), loc='upper left')
        logging.info('MP: Created GE results legend')

        # Remove x and y axes for both subplots
        ax_comparison[0].axis('off')
        ax_comparison[1].axis('off')
        logging.info('MP: Removed axis features')

        # Set titles for the subplots
        ax_comparison[0].set_title('Overlay of Normalised 2019 \n GE Voter Turnout and'
                                + ' Income Deprivation Rate \n in Local Authorities')
        ax_comparison[1].set_title('2019 UK General Election Results')
        logging.info('MP: Added figure titles')

        return fig_comparison


    def depriv_plot(self):
        '''Plot the local authorities, using the deprivation rate to colour
        in each boundary. Authorities with no data are coloured grey'''

        # Choose the colum to use for the colour and get a colourmap
        deprivation_data_column = 'Income deprivation rate'

        # Normalise the scale to the highest and lowest values in the column
        norm = Normalize(vmin=self.deprivation[deprivation_data_column].min(),
                         vmax=self.deprivation[deprivation_data_column].max())

        return deprivation_data_column

    def normalise_results(self, dataframe, column):
        '''Takes a dataframe and column name, and creates a new column in the
        dataframe for normalised data (0-1). This makes it easier to overlay
        plots with colorbars'''

        normalised_column = column + "_norm"
        dataframe[normalised_column] = ((dataframe[column] - dataframe[column].min())
                                      / (dataframe[column].max() - dataframe[column].min()))

        return normalised_column


if __name__ == "__main__":
    local_auth = LocalAuthorityData()
    #auth_plot = local_auth.plot()

    elec = ElectoralResults()
    #elec.party_colours()
    #elec.voter_turnout()
    #elec.plot_turnout()
    #elec_plot = elec.plot_turnout()

    graphs = MultiplePlots(elec.electoral_gdf, local_auth.local_authority_gdf)
    #fig_indiv, fig_overlay, fig_comparison = graphs.plot()
    fig_comparison = graphs.plot()
    plt.tight_layout()
    plt.show()
