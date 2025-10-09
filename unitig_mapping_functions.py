def unitig_samples(unitig):
    '''Locate which camples contain the untig
    This is useful for later analytics'''
    import pandas as pd
    replicating_features = [unitig]
    replicating_features =  [int(feature) for feature in replicating_features]
    all_features = set(range(0,5458393))
    rows_to_skip = all_features - set(replicating_features) 
    rows_to_skip = rows_to_skip 
    rows_to_skip = list(rows_to_skip)
    rows_to_skip = [x+1 for x in rows_to_skip]            
    replicating_features_frame = pd.read_csv("/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/DBGWAS/bugwas_input.unique_rows.binary", sep=' ',  skiprows = rows_to_skip, index_col= "ps")
    #replicating_features_frame = 1 - replicating_features_frame
    #save for plotting on tree
    replicating_features_frame.T.to_csv("/Users/kvk22/Desktop/Scripts/Philogeny_aware_ML/proofs/aztreonam_markers.txt")
    unitig_samples = (replicating_features_frame.columns[(replicating_features_frame == 1).any(axis=0)])
    print(f"Located unitig with frequency {replicating_features_frame.T.mean().values[0]}")
    return unitig_samples, replicating_features_frame


def write_unitig_string_fasta(unitig_string):
    '''Saves the unitig as fasta file and returns the path to it'''
    fasta_path_to_unitig = "/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/genomes/unitig_string.fasta"
    with open(fasta_path_to_unitig, "w") as f:
        f.write(">unitig_string\n")
        f.write(unitig_string + "\n")
    return fasta_path_to_unitig

def check_presence(unitig_string, unitig_samples, unitig_string_fasta):
    import os
    import subprocess
    import pandas as pd

    all_annotatos = []
    for unitig_sample in unitig_samples:
        #file paths
        bwa_path = "/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/genomes/bwa/bwa"
        sample_fasta = f"/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/genomes/{unitig_sample}"
        unitig_fasta = unitig_string_fasta

        #Check if all files exist
        for file in [bwa_path, sample_fasta, unitig_fasta]:
            if not os.path.exists(file):
                print(f"Sample {sample_fasta} missing")
                continue

        #command
        cmd = [bwa_path, "fastmap", "-w", "10", "-l", str(len(unitig_string)), sample_fasta, unitig_fasta]

        # Run command 
        result = subprocess.run(cmd, capture_output=True, text=True)
        all_annotatos.append(result.stdout)
        
    #filter strings to clear only the needed parts
    processed = [
        entry.split('\nEM', 1)[-1].split('\n', 1)[0] if '\nEM' in entry else entry.split('\n', 1)[0]
        for entry in all_annotatos
    ]
    split_entries = [entry.split('\t') for entry in processed]
    split_entries = [entry.split('\t') for entry in processed]
    df = pd.DataFrame(split_entries)
    df = df.drop(columns=[0,1])
    df = df.iloc[:, :3]
    df = df.dropna()
    df = df.rename(columns={df.columns[0]: "unitig length", df.columns[1]: "unitig count", df.columns[2]: f"occurance_1"})
    df = df.reset_index(drop=True)
    return df

def extract_gene(matching_annotations):
    import re
    import os
    import gffpandas.gffpandas as gffpd
    import warnings
    import pandas as pd
    warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
    contexts_across_samples = {}
    locations_across_samples = {}

    for sample_itter in matching_annotations["occurance_1"]:
        #format the strings for reading
        s = sample_itter
        match = re.match(r"^(.*?)([+-])(\d+)$", s)
        if match:
            contig= match.group(1).strip(":")
            location_on_contig = match.group(3)     
            if ".con." in contig:
                sample_fna = contig.split(".con.")[0]

        #read the GFF file and keep only the annotations
        path = "/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/Arthemis/annotations/"
        gff_file = path + sample_fna + ".fna.gff"
        if os.path.exists(gff_file):
            annotation = gffpd.read_gff3(gff_file)
        else:
            print(f"File not found: {gff_file}")
            continue
        annotation = annotation.df
        annotation = annotation.dropna()

        #extract the gene of interest
        unitig_location = annotation[annotation["seq_id"] == contig]
        location_on_contig = int(location_on_contig)
        location_gene = unitig_location[
            (unitig_location["start"] <= location_on_contig) & (unitig_location["end"] >= location_on_contig)
        ]

        #extract genetic context
        matching_indices = location_gene.index.tolist()
        # Collect indices for context: 3 before and 3 after each match
        context_indices = []
        for idx in matching_indices:
            start = max(0, idx - 3)
            end = idx + 4 
            context_indices.extend(range(start, end))
        # Remove duplicates and keep only valid indices
        context_indices = sorted(set(i for i in context_indices if i in unitig_location.index))
        context = unitig_location.loc[context_indices]

        #rename the unitig indexes
        context = context.rename(index={5134: "unitig_index"})
        location_gene = location_gene.rename(index={5134: "unitig_index"})

        #add to dictionary
        if len(context)>0:
            contexts_across_samples.update({sample_fna: context})
            locations_across_samples.update({sample_fna: location_gene})
    return contexts_across_samples, locations_across_samples


def extract_gene_seqeunces(locations_across_samples, unitig_string):
    '''Get the DNA and protein sequences for this unitig'''
    from Bio import SeqIO
    from Bio.Seq import Seq

    genes = {}

    for sample in locations_across_samples:
        location = locations_across_samples[sample]

        #read genome
        genome_path = "/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/genomes/"
        fasta_path = genome_path + sample + ".fna"
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id == location["seq_id"].values[0]:
                seqeunce = record.seq
                break

        #confirm contig is in sequence otherwise reverse complement
        if unitig_string.lower() in str(seqeunce):
            seqeunce = str(seqeunce)
        else:
            seqeunce_reversed = str(seqeunce.reverse_complement())
            #id thill cnnot locate raise error
            if unitig_string.lower() not in seqeunce_reversed:
                raise ValueError("unitig_string is not present in the sequence")
        
        #get gene string 
        start_value = int(location["start"].values[0])-1
        end_value = int(location["end"].values[0])
        gene_string = str(seqeunce[start_value:end_value])

        #translate to protein
        protein_seq = str(Seq(gene_string).translate())
        if protein_seq[0] != "M":
            raise ValueError(f"Something wrong with translation. Protein sequence: {protein_seq}")

        #add to dictionary
        genes[sample] = [gene_string, protein_seq]
        
    return genes

def clustal_alignment(extracted_genes_proteins):
    import os
    import tempfile
    from Bio import SeqIO
    import pandas as pd
    import subprocess
    import shutil

    #align genes
    #Create a temporary directory for FASTA files
    temp_dir = tempfile.mkdtemp()
    #Write each sequence to a FASTA file in the temp directory
    fasta_path = os.path.join(temp_dir, "msa_input.fasta")
    with open(fasta_path, "w") as fasta_file:
        for key in extracted_genes_proteins:
            fasta_file.write(f">{key}\n{extracted_genes_proteins[key][0]}\n")

    #Run Clustal Omega to generate the alignment (output in FASTA format)
    msa_output_path = os.path.join(temp_dir, "msa_output.fasta")
    subprocess.run([
        "clustalo", "-i", fasta_path, "-o", msa_output_path, "--force", "--outfmt=fasta"
    ], check=True)

    #Read the aligned sequences from the output FASTA
    aligned_seqs = {}
    for record in SeqIO.parse(msa_output_path, "fasta"):
        aligned_seqs[record.id] = str(record.seq)

    #Convert the alignment to a pandas DataFrame
    alignment_genes = pd.DataFrame(
        [list(seq) for seq in aligned_seqs.values()],
        index=list(aligned_seqs.keys())
    ).T  # Transpose for SNP analysis

    #Clean up: delete the temporary directory and files
    shutil.rmtree(temp_dir)


    #align proteins
    temp_dir = tempfile.mkdtemp()
    #Write each sequence to a FASTA file in the temp directory
    fasta_path = os.path.join(temp_dir, "msa_input.fasta")
    with open(fasta_path, "w") as fasta_file:
        for key in extracted_genes_proteins:
            fasta_file.write(f">{key}\n{extracted_genes_proteins[key][1]}\n")

    #Run Clustal Omega to generate the alignment (output in FASTA format)
    msa_output_path = os.path.join(temp_dir, "msa_output.fasta")
    subprocess.run([
        "clustalo", "-i", fasta_path, "-o", msa_output_path, "--force", "--outfmt=fasta"
    ], check=True)

    #Read the aligned sequences from the output FASTA
    aligned_seqs = {}
    for record in SeqIO.parse(msa_output_path, "fasta"):
        aligned_seqs[record.id] = str(record.seq)

    #Convert the alignment to a pandas DataFrame
    alignmed_proteins = pd.DataFrame(
        [list(seq) for seq in aligned_seqs.values()],
        index=list(aligned_seqs.keys())
    ).T  # Transpose for SNP analysis
    #Clean up: delete the temporary directory and files
    shutil.rmtree(temp_dir)


    return alignment_genes, alignmed_proteins