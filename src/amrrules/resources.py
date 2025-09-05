"""Resource management for AMRFinderPlus data files."""

import csv
import urllib.request
from pathlib import Path
from typing import Optional

class ResourceManager:
    """Manages cached resource files for validation."""
    
    def __init__(self):
        """Initialize resource manager with default resource directory."""
        self.dir = Path(__file__).parent / "resources"
        # Create the resources directory if it doesn't exist
        #self.dir.mkdir(parents=True, exist_ok=True)
        self._amrfp_card_convert_cache: Optional[list] = None
        self._amrfp_db_version: Optional[str] = None
        self._refseq_nodes_cache: Optional[dict] = None
    
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
    
    def download_amrfp_resources(self):
        """
        Download AMRFinderPlus reference files into the resource directory.
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
            print("AMRFinderPlus resources downloaded successfully.")
            # Make sure to update the cached version
            self.get_amrfp_db_version()
        else:
            print("Warning: Some AMRFinderPlus resources could not be downloaded.")
        
        return success
    
    def setup_all_resources(self):
        """
        Download and set up all required resources.
        """
        amrfp_success = self.download_amrfp_resources()

        if amrfp_success:
            print("All resources have been successfully set up.")
            return True
        else:
            print("Warning: Some resources could not be set up properly.")
            return False