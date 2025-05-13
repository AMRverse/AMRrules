import csv
import urllib.request
import io

def download_and_parse_reference_gene_hierarchy(url):
    print(f"Downloading Reference Gene Hierarchy from {url}...")
    
    # Use urllib to download the file
    with urllib.request.urlopen(url) as response:
        content = response.read().decode('utf-8')  # Decode the response as UTF-8
    
    # Parse the downloaded content as a TSV file
    content_io = io.StringIO(content)
    reader = csv.DictReader(content_io, delimiter='\t')
    
    amrfp_nodes = {}
    for row in reader:
        node_id = row.get('node_id')
        parent_node = row.get('parent_node_id')
        amrfp_nodes[node_id] = parent_node

    print(f"Downloaded and parsed {len(amrfp_nodes)} rows from the Reference Gene Hierarchy.")
    return amrfp_nodes

def parse_rules_file(rules_file, tool):
    rules = []
    #if tool == 'amrfp':
    #    key_value = 'nodeID'
    with open(rules_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rules.append(row)
    return rules