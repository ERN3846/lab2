import requests
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import streamlit as st
import numpy as np

def retrieve_ppi_string(target_protein):
    """
    Retrieve protein-protein interaction data from STRING DB
    """
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606  # Human species ID
    }
    
    try:
        response = requests.get(string_url, params=params)
        response.raise_for_status()
        
        if not response.json():
            st.warning(f"No interactions found for protein {target_protein}")
            return pd.DataFrame()
        
       # Convert JSON response to DataFrame            
        network_df = pd.json_normalize(response.json())
        
        if network_df.empty:
            st.warning("Retrieved data is empty")
            return pd.DataFrame()
            
        return network_df
        
    except requests.exceptions.RequestException as e:
        st.error(f"Error retrieving data from STRING: {str(e)}")
        return pd.DataFrame()

def retrieve_ppi_biogrid(target_protein):
    """
    Retrieve protein-protein interaction data from BioGRID
    """
    bio_url = "https://webservice.thebiogrid.org/interactions/"
    params = {
        "searchNames": "true",
        "geneList": target_protein,
        "includeInteractors": "true",
        "taxId": 9606,   # Human species ID
        "format": "json",
        "accesskey": "0fb13a467d24525f9a83d5d1efb3f211"
    }
    
    try:
        response = requests.get(bio_url, params=params)
        response.raise_for_status()
 
         # Check if response contains data       
        data = response.json()
        if not data:
            st.warning(f"No interactions found for protein {target_protein} in BioGRID.")
            return pd.DataFrame()
        
        # Convert JSON response to DataFrame  
        network_df = pd.DataFrame.from_dict(data, orient='index')
    
        # Check for required columns and return DataFrame    
        if 'OFFICIAL_SYMBOL_A' in network_df.columns and 'OFFICIAL_SYMBOL_B' in network_df.columns:
            return network_df
        else:# BioGRID
            st.warning("BioGRID data does not contain expected columns.")
            print("Available columns:", network_df.columns)
            return pd.DataFrame()

    except requests.exceptions.RequestException as e:
        st.error(f"Error retrieving data from BioGRID: {str(e)}")
        return pd.DataFrame()
    except ValueError:
        st.error("Error: BioGRID API response could not be parsed as JSON.")
        return pd.DataFrame()

def process_network_data(df, database):
    """
    Process the network data based on the database source
    """
    if df.empty:
        return pd.DataFrame()
        
    if database == "STRING":
        if 'preferredName_A' in df.columns and 'preferredName_B' in df.columns:
            return df[['preferredName_A', 'preferredName_B']]
    else:  
        if 'OFFICIAL_SYMBOL_A' in df.columns and 'OFFICIAL_SYMBOL_B' in df.columns:
            return df[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
            
    st.error(f"Required columns not found in {database} data")
    return pd.DataFrame()

def generate_network(df, database):
    """
    Generate network graph from processed data
    """
    if df.empty:
        return nx.Graph()
        
    source_col = 'preferredName_A' if database == "STRING" else 'OFFICIAL_SYMBOL_A'
    target_col = 'preferredName_B' if database == "STRING" else 'OFFICIAL_SYMBOL_B'
    
    return nx.from_pandas_edgelist(df, source_col, target_col)

def get_top_proteins(G, n=5):
    """
    Get top N proteins by degree centrality
    """
    degree_centrality = nx.degree_centrality(G)
    return sorted(degree_centrality.items(), key=lambda x: -x[1])[:n]

def plot_network(G, top_proteins=None):
    """
    Plot network graph with highlighted top proteins
    """
    plt.figure(figsize=(12, 8))
    
    # Generate layout    
    layout = nx.spring_layout(G, k=1/np.sqrt(G.number_of_nodes()), seed=42)
    
    # Draw base network
    nx.draw(G, layout,
            node_color='lightblue',
            node_size=1000,
            with_labels=True,
            font_size=8,
            font_weight='bold',
            edge_color='gray',
            alpha=0.7)
    
    # Highlight top proteins if provided
    if top_proteins:
        top_nodes = [node for node, _ in top_proteins]
        nx.draw_networkx_nodes(G, layout,
                             nodelist=top_nodes,
                             node_color='orange',
                             node_size=1000)
    
    return plt.gcf()

# Streamlit UI
st.title("Protein-Protein Interaction Network Analysis")

# User inputs
target_protein = st.text_input("Enter Protein ID:", 
                             help="Enter a protein identifier (e.g., TP53, BRCA1)")
database_option = st.selectbox("Select Database:", ["STRING", "BioGRID"])
top_n = st.slider("Number of top proteins to highlight:", min_value=1, max_value=10, value=5)

if st.button("Analyze Network"):
    if not target_protein:
        st.error("Please enter a protein ID")
    else:
        with st.spinner(f"Retrieving PPI data from {database_option}..."):
            # Retrieve and process data
            if database_option == "STRING":
                df = retrieve_ppi_string(target_protein)
            else:
                df = retrieve_ppi_biogrid(target_protein)
                
            processed_df = process_network_data(df, database_option)
            
            if not processed_df.empty:
                # Create network and analyze
                G = generate_network(processed_df, database_option)
                
                # Display basic network statistics    
                st.subheader("Network Statistics")
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Number of Nodes", G.number_of_nodes())
                with col2:
                    st.metric("Number of Edges", G.number_of_edges())
                with col3:
                    st.metric("Network Density", round(nx.density(G), 4))
                
                # Get and display top proteins
                top_proteins = get_top_proteins(G, top_n)
                
                st.subheader("Top Proteins by Centrality")
                top_proteins_df = pd.DataFrame(top_proteins, 
                                             columns=['Protein', 'Centrality Score'])
                st.dataframe(top_proteins_df)
                
                # Plot network
                st.subheader("Network Visualization")
                fig = plot_network(G, top_proteins)
                st.pyplot(fig)
                plt.close()
                
                # Option to download the network image 
                if st.button("Download Network Image"):
                    plt.figure(figsize=(12, 8))
                    fig = plot_network(G, top_proteins)
                    plt.savefig(f"{target_protein}_network.png", 
                              bbox_inches='tight', 
                              dpi=300)
                    st.success(f"Network image saved as {target_protein}_network.png")
                    plt.close()
