import os
import re
from collections import defaultdict

# ==========================================
# CONFIGURATION
# ==========================================
# Path to the directory containing your GFF/GFF3 files
GFF_DIR = "./gffs"

# Path to the directory containing OrthoFinder MSA files (e.g., .fa or .mfa)
MSA_DIR = "OrthoFinder/Results_Jun25_1/MultipleSequenceAlignments"

# Path to the directory where formatted Malin files will be saved
OUTPUT_DIR = "./malin_inputs"

# Optional: Manual mapping of GFF filename prefix to Malin organism tag
# If a file is named 'Homo_sapiens.gff', it defaults to organism 'Homo_sapiens'
# You can override it here if you need short codes (e.g., {'Homo_sapiens': 'Hsap'})
ORGANISM_MAP = {}
# ==========================================

def clean_id(identifier):
    """
    Removes all non-alphanumeric characters and converts to lowercase 
    to prevent punctuation mismatches (e.g., 'WGS:JAIUWE' vs 'WGS_JAIUWE').
    """
    if not identifier:
        return ""
    # Strip any leading 'rna-' or 'cds-' prefixes if they are inconsistent
    identifier = re.sub(r'^(rna|cds)-', '', identifier, flags=re.IGNORECASE)
    return re.sub(r'[^a-zA-Z0-9]', '', identifier).lower()

def parse_gff_introns(gff_path, organism_tag):
    """
    Parses GFF and indexes introns using normalized, punctuation-free ID keys.
    """
    cds_dict = defaultdict(list)
    print(f"Processing GFF for organism: {organism_tag}...")
    
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            # Extract potential IDs
            id_choices = []
            for tag in ['Parent=', 'orig_transcript_id=', 'protein_id=', 'ID=']:
                match = re.search(f'{tag}([^;]+)', attributes)
                if match:
                    id_choices.append(match.group(1).strip())
            
            if id_choices:
                # Use the normalized version of the first ID as primary dictionary key
                primary_key = clean_id(id_choices[0])
                cds_dict[primary_key].append((start, end, strand))
                
                # Link all other normalized aliases to the same data
                for alias in id_choices[1:]:
                    norm_alias = clean_id(alias)
                    if norm_alias and norm_alias != primary_key:
                        cds_dict[norm_alias] = cds_dict[primary_key]

    intron_map = {}
    for tracking_key, cds_list in cds_dict.items():
        if not cds_list or tracking_key in intron_map:
            continue
            
        strand = cds_list[0][2]
        if strand == '+':
            cds_list.sort(key=lambda x: x[0])
        else:
            cds_list.sort(key=lambda x: x[0], reverse=True)
            
        introns = []
        cumulative_length = 0
        for i in range(len(cds_list) - 1):
            cds_len = cds_list[i][1] - cds_list[i][0] + 1
            cumulative_length += cds_len
            introns.append(str(cumulative_length))
            
        intron_map[tracking_key] = {
            'organism': organism_tag,
            'introns': introns
        }
        
    return intron_map

def find_matching_data(header, all_intron_data):
    """
    Normalizes the FASTA header and checks if any GFF key 
    is a normalized substring of it.
    """
    # Clean the full FASTA header line
    clean_fasta_header = clean_id(header.split()[0].lstrip('>'))
    
    # Check for direct match or substring match using normalized keys
    for gff_norm_key in all_intron_data:
        if gff_norm_key in clean_fasta_header:
            return all_intron_data[gff_norm_key]
            
    return None
    
def process_msas(msa_dir, output_dir, all_intron_data):
    """
    Reads MSA files, updates headers with Malin tags, and saves them.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    msa_files = [f for f in os.listdir(msa_dir) if f.endswith(('.fa', '.fasta', '.mfa'))]
    
    for msa_file in msa_files:
        input_path = os.path.join(msa_dir, msa_file)
        output_path = os.path.join(output_dir, msa_file)
        
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            current_header = ""
            sequence_lines = []
            
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous record
                    if current_header:
                        write_malin_record(current_header, sequence_lines, outfile, all_intron_data)
                    current_header = line
                    sequence_lines = []
                else:
                    sequence_lines.append(line)
                    
            # Process the final record in the file
            if current_header:
                write_malin_record(current_header, sequence_lines, outfile, all_intron_data)

def write_malin_record(header, sequence_lines, outfile, all_intron_data):
    """
    Appends the required organism and intron tags to the header and writes to file.
    """
    match = find_matching_data(header, all_intron_data)
    
    if match:
        org_tag = match['organism']
        intron_list = match['introns']
        intron_str = ",".join(intron_list) if intron_list else ""
        
        # Build Malin compliant header
        new_header = f"{header} /organism={org_tag} {{i {intron_str} i}}"
    else:
        # Fallback if gene wasn't found in GFF maps
        print(f"Warning: Could not find structural data for sequence: {header}")
        new_header = f"{header} /organism=Unknown {{i i}}"
        
    outfile.write(f"{new_header}\n")
    outfile.write(f"\n".join(sequence_lines) + "\n")

def main():
    # Step 1: Parse all GFF files
    all_intron_data = {}
    if not os.path.exists(GFF_DIR):
        print(f"Error: GFF directory '{GFF_DIR}' not found. Please create it and add your GFF files.")
        return
        
    gff_files = [f for f in os.listdir(GFF_DIR) if f.endswith(('.gff', '.gff3'))]
    
    for gff_file in gff_files:
        base_name = os.path.splitext(gff_file)[0]
        # Use short map tag if defined, otherwise use the filename base
        org_tag = ORGANISM_MAP.get(base_name, base_name)
        
        gff_path = os.path.join(GFF_DIR, gff_file)
        organism_introns = parse_gff_introns(gff_path, org_tag)
        all_intron_data.update(organism_introns)
        
    # Step 2: Format MSA files
    if not os.path.exists(MSA_DIR):
        print(f"Error: MSA directory '{MSA_DIR}' not found. Please verify your OrthoFinder output path.")
        return
        
    print("\nFormatting MSA files for Malin...")
    process_msas(MSA_DIR, OUTPUT_DIR, all_intron_data)
    print(f"Done! Formatted files are located in: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
