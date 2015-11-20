Physical protein-protein interaction (PPI) network 

Interaction data are derived from the BioGRID database (version 3.2.121, downloaded on 02.02.2015).

The following filtering steps were performed on the entire dataset of interactions from the file BIOGRID-ORGANISM-3.2.121.mitab/BIOGRID-ORGANISM-Homo_sapiens-3.2.121.mitab.txt.  
1. Only interactions between human proteins were considered (taxid:9606).
2. Only direct interactions and physical associations were considered (psi-mi:MI:0407, psi-mi:MI:0915). 
3. Interactions were counted symmetrically (A-B equals B-A).
4. Unique PMIDs (PubMed references), BioGRID:IDs and interaction types were counted for every interaction. 
5. Alternative aliases for each interacting protein were discarded (updated 17.02.15). 

Two files with interactions are provided.
- PPI_network_BioGrid_ALL.txt contains all protein-protein interactions (n=138,722)
- PPI_network_BioGrid_HC.txt contains all high-confidence protein-protein interactions with at least two independent studies (PMIDs) (n=21,038). 

Contact: Juri Reimand (Juri.Reimand@utoronto.ca)