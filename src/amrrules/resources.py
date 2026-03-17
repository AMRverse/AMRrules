"""Resource management for AMRFinderPlus data files."""

import csv
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import obonet
import tempfile
import tarfile

class ResourceManager:
    """Manages external resource files required for assigning and annotating rules."""
    
    def __init__(self):
        """Initialize resource manager with default resource directory."""
        self.dir = Path(__file__).parent / "resources"
        self._amrfp_card_convert_cache: Optional[list] = None
        self._amrfp_db_version: Optional[str] = None
        self._refseq_nodes_cache: Optional[dict] = None
        self._card_drug_map: Optional[dict] = None
    
    def setup_all_resources(self):
        """
        Download and set up all required external resources, from AMRFP and CARD databases.
        """
        amrfp_success = self.download_amrfp_resources()
        card_success = self.download_card_archives()

        if amrfp_success and card_success:
            print("All resources have been successfully set up.")
            return True
        else:
            return False

    # Functions for downloading AMRFP and CARD files
    def download_amrfp_resources(self):
        """
        Download AMRFinderPlus Ref Gene Hierarchy and the database version number into the resource directory.
        """
        # URLs for the AMRFinderPlus resources
        amrfp_nodes_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt'
        amrfp_version_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/version.txt'
        
        # Files to download
        file_urls = {
            'ReferenceGeneHierarchy.txt': amrfp_nodes_url,
            'version.txt': amrfp_version_url
        }
        
        success = True
        for filename, url in file_urls.items():
            target_path = self.dir / filename
            print(f"Downloading {filename} from {url}...")
            
            try:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    
                with open(target_path, 'wb') as f:
                    f.write(content)
                
                print(f"Successfully downloaded {filename}")
                
                # Get the database version
                if filename == 'version.txt':
                    self._amrfp_db_version = content.decode('utf-8').strip()
                    print(f"AMRFinderPlus database version: {self._amrfp_db_version}")
                    
            except Exception as e:
                print(f"Error downloading {filename}: {e}")
                success = False
        
        if success:
            print("AMRFinderPlus database files downloaded successfully.")
            # Make sure to update the amrfp version number
            self.get_amrfp_db_version()
        else:
            print("Warning: The AMRFinderPlus resources could not be downloaded.")
        
        return success

    def download_card_archives(self):
        """
        Download and extract CARD ontology and data files into the resource directory.
        """
        # URLs for the CARD database files
        card_ontology_url = "https://card.mcmaster.ca/download/5/ontology-v4.0.1.tar.bz2"
        card_data_url = "https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2"

        # Files to extract from each archive
        card_ontology_files = ["aro.obo", "ncbi_taxonomy.tsv"]
        card_data_files = ["aro_categories.tsv"]

        # Download and extract files from both archives
        ontology_success = self._download_and_extract(card_ontology_url, card_ontology_files)
        data_success = self._download_and_extract(card_data_url, card_data_files)
        
        if ontology_success and data_success:
            print("CARD archives downloaded and extracted successfully.")
            
        # Verify that the files exist
        missing_files = []
        for file_name in card_ontology_files + card_data_files:
            file_path = self.dir / file_name
            if not file_path.exists():
                missing_files.append(file_name)
                
        if missing_files:
            print(f"Warning: The following files are missing: {', '.join(missing_files)}")
            return False
        else:
            print("All required files are present in the resources directory.")
            return ontology_success and data_success
    
    def _download_and_extract(self, url, files_to_extract):
        """Function to help download CARD archives and extract specific files, saving them into the resources directory."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            # Download the archive
            print(f"Downloading {url}...")
            try:
                with urllib.request.urlopen(url) as response:
                    temp_file.write(response.read())
                    temp_file.flush()
            except Exception as e:
                print(f"Error downloading {url}: {e}")
                return False

            temp_path = temp_file.name
        
        try:
            # Extract specific files
            with tarfile.open(temp_path, "r:bz2") as tar:
                all_members = tar.getmembers()
                # For each file we want to extract
                for target_file in files_to_extract:
                    found = False
                    # Look for exact match or file within subdirectory
                    for member in all_members:
                        basename = Path(member.name).name
                        if basename == target_file:
                            print(f"Extracting {member.name} as {target_file}...")
                            # Extract but rename to the target filename
                            member_obj = tar.extractfile(member)
                            if member_obj:
                                with open(self.dir / target_file, 'wb') as f:
                                    f.write(member_obj.read())
                                found = True
                                break
                    
                    if not found:
                        print(f"Warning: Could not find {target_file} in the archive")
                        
            return True
        except Exception as e:
            print(f"Error extracting files: {e}")
            return False
        finally:
            # Clean up the temp file
            try:
                Path(temp_path).unlink(missing_ok=True)
            except:
                pass

    # Functions for parsing AMRFP and CARD resources into data structures used elsewhere
    def refseq_nodes(self) -> dict:

        if self._refseq_nodes_cache is None:
            refseq_file = self.dir / "ReferenceGeneHierarchy.txt"
            if refseq_file.exists():
                self._refseq_nodes_cache = self._load_refseq_nodes(str(refseq_file))
            else:
                # Return empty dict if file doesn't exist yet
                self._refseq_nodes_cache = {}
        return self._refseq_nodes_cache

    def _load_refseq_nodes(self, node_file: str):
        """Load RefSeq nodes from the given file."""
        self._refseq_nodes_cache = {}

        refseq_hierarchy = csv.DictReader(open(node_file, 'r'), delimiter='\t')
        for row in refseq_hierarchy:
            node_id = row.get('node_id')
            parent_node = row.get('parent_node_id')
            self._refseq_nodes_cache[node_id] = parent_node

        return self._refseq_nodes_cache

    def get_amrfp_card_conversion(self) -> dict:

        if self._amrfp_card_convert_cache is None:
            card_file = self.dir / "amrfp_to_card_drugs_classes.txt"
            if card_file.exists():
                self._amrfp_card_convert_cache = self._load_amrfp_card_conversion(str(card_file))
            else:
                # Return empty dict if file doesn't exist yet
                self._amrfp_card_convert_cache = {}
        return self._amrfp_card_convert_cache

    def _load_amrfp_card_conversion(self, card_file: str):
        """
        Load the dictionary that converts AMRFP Subclasses to the CARD drug and drug class ontology.
        Also create a dict that just converts CARD drugs to their classes.
        """
        self._amrfp_card_convert_cache = {}

        with open(card_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                amrfp_subclass = row.get('AFP_Subclass')
                card_drug = row.get('CARD drug')
                card_class = row.get('CARD drug class')
                self._amrfp_card_convert_cache[amrfp_subclass] = {
                    'drug': card_drug,
                    'class': card_class
                }

        return self._amrfp_card_convert_cache

    def get_amrfp_db_version(self) -> str:
        """
        Get the AMRFinderPlus database version from the downloaded version.txt file.
        
        Returns:
            str: The version string or "Unknown" if the file is not available
        """
        version_file = self.dir / "version.txt"
        if version_file.exists():
            try:
                with open(version_file, 'r') as f:
                    self._amrfp_db_version = f.read().strip()
                return self._amrfp_db_version
            except Exception as e:
                print(f"Error reading AMRFinderPlus version file: {e}")
                return "Unknown"
        else:
            print("AMRFinderPlus version file not found.")
            return "Unknown"
    
    def _parse_obo_for_descendants(self, obo_file_path: str, term_id: str) -> List[str]:
        """Find all descendants of a given term in the OBO file."""
        children = []
        with open(obo_file_path, 'r') as file:
            current_term = None
            for line in file:
                line = line.strip()
                if line.startswith("[Term]"):
                    current_term = None
                elif line.startswith("id: "):
                    current_term = line.split("id: ")[1]
                elif line.startswith("is_a: ") and current_term:
                    parent = line.split("is_a: ")[1].split(" ! ")[0]
                    if parent == term_id:
                        children.append(current_term)
        return children
    
    def _extract_card_drugs(self, obo_file_path: str, categories_file_path: str) -> List[Tuple[str, str, str]]:
        """Extract drug names and classes from CARD ontology files."""
        # Load the OBO file
        card_ontology = obonet.read_obo(obo_file_path)
        
        # Load the categories TSV
        with open(categories_file_path, 'r', newline='') as file:
            card_categories = csv.DictReader(file, delimiter='\t')
            
            # Get drug classes from categories
            drug_classes: Dict[str, str] = {}
            for row in card_categories:
                if row['ARO Category'] == 'Drug Class':
                    drug_classes[row['ARO Name']] = row['ARO Accession']
        
        # Get term names
        id_to_name = {id_: data.get("name") for id_, data in card_ontology.nodes(data=True)}
        
        # Collect results
        output_dict = {}
        
        # Extract drugs for each drug class
        for drug_class, aro_accession in drug_classes.items():
            try:
                children = self._parse_obo_for_descendants(obo_file_path, aro_accession)
                for child in children:
                    if child in id_to_name:
                        child_name = id_to_name[child]
                        output_dict[child_name] = drug_class
            except Exception as e:
                print(f"Error processing {aro_accession} for {drug_class}: {e}")
                continue
        
        # Additional specific betalactam and other drug classes
        betalac_aros = [
            'ARO:3009105', 'ARO:3009106', 'ARO:3009107', 'ARO:3009108', 
            'ARO:3009109', 'ARO:3009123', 'ARO:3009124', 'ARO:3009125',
            'ARO:3000035', 'ARO:3007783', 'ARO:0000022', 'ARO:3007629',
            'ARO:3000707'
        ]
        
        for aro_accession in betalac_aros:
            if aro_accession in id_to_name:
                drug_class = id_to_name[aro_accession]
                children = self._parse_obo_for_descendants(obo_file_path, aro_accession)
                for child in children:
                    if child in id_to_name:
                        child_name = id_to_name[child]
                        output_dict[child_name] = drug_class

        return output_dict

    def get_card_drug_class_map(self):
        if self._card_drug_map is None:
            obo_file = self.dir / "aro.obo"
            categories_file = self.dir / "aro_categories.tsv"
            self._card_drug_map = self._extract_card_drugs(str(obo_file), str(categories_file))
        return self._card_drug_map
