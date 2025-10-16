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

<<<<<<< HEAD
def extract_gene(matching_annotations, unitig_string = "CGGCAGCGTCAGATGTGTATAAGAGACAGTA"):
=======
def extract_gene(matching_annotations):
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480
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
<<<<<<< HEAD
        if len(location_gene) == 0: print(f"Unitig outside gene in sample {sample_fna}")
        if len(location_gene) > 0:
            location_gene.index = ["unitig_index"]
        
        #extract genetic context
        unitig_location = unitig_location[~((unitig_location["start"] <= location_on_contig) & (unitig_location["end"] >= location_on_contig))]
        downstream = unitig_location[unitig_location["start"] < location_on_contig].head(3)
        upstream = unitig_location[unitig_location["start"] > location_on_contig].head(3)
        context = pd.concat([downstream,location_gene, upstream])
=======

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
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480

        #add to dictionary
        if len(context)>0:
            contexts_across_samples.update({sample_fna: context})
            locations_across_samples.update({sample_fna: location_gene})
    return contexts_across_samples, locations_across_samples


<<<<<<< HEAD

def extract_gene_seqeunces_unitig(alternative_proteins_locations):
    from tqdm import tqdm
    import os
    from Bio import SeqIO
    from Bio.Seq import Seq

    extracted_genes_proteins = {}

    for sample in tqdm(alternative_proteins_locations, desc="Extracting genes from all samples"):
        genes_in_samples = []
        for entry in alternative_proteins_locations[sample]:
            start = int(entry.split("\t")[3])-1
            end = int(entry.split("\t")[4])
            contig = entry.split("\t")[0]

            path_to_annotations = "/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/Arthemis/annotations/"
            gff_path = path_to_annotations + sample + ".gff"
            
            #check if file exists
            if os.path.exists(gff_path):
                gff_lines = []
                fasta_lines = []
                in_fasta = False
                with open(gff_path, "r") as f:
                    for line in f:
                        if line.startswith("##FASTA"):
                            in_fasta = True
                            continue
                        if in_fasta:
                            fasta_lines.append(line.rstrip())
                        else:
                            gff_lines.append(line.rstrip())
                # gff_lines contains GFF annotation lines
                dna_sequence_str = "".join(fasta_lines)
                split_sequences = dna_sequence_str.split(">")
                split_sequences = [seq for seq in split_sequences if seq]

                #extract the gene from the correct contig
                matching_entries = [seq for seq in split_sequences if contig in seq]
                cleaned_entries = [seq.replace(contig, "") for seq in matching_entries]
                gene_string = cleaned_entries[0][start:end]

                #translate to protein
                protein_seq = str(Seq(gene_string).translate())
                if protein_seq[0] != "M":
                    gene_string = str(Seq(gene_string).reverse_complement())
                    protein_seq = str(Seq(gene_string).translate())
                    if protein_seq[0] != "M":
                        raise ValueError(f"Something wrong with translation. Protein sequence: {protein_seq}")
                
                genes_in_samples.append([gene_string, protein_seq])
                del gene_string
                del protein_seq

            extracted_genes_proteins[sample] = genes_in_samples
    return extracted_genes_proteins


def unitigs_seqneces_not_in_genes(unitig_string, matching_annotations, genes):

    import re
    import os
    from Bio.Seq import Seq
    from Bio import SeqIO


    unitig_string = "CGGCAGCGTCAGATGTGTATAAGAGACAGTA"


    extracted_genes_proteins = {}

    contigs_of_interest = matching_annotations["occurance_1"].str.split(":", expand=True)
    contigs_of_interest.columns = ["col1", "col2"]

    for sample_itter in matching_annotations["occurance_1"]:
        #format the strings for reading
        s = sample_itter
        match = re.match(r"^(.*?)([+-])(\d+)$", s)
        if match:
            contig= match.group(1).strip(":")
            location_on_contig = match.group(3)     
            if ".con." in contig:
                sample_fna = contig.split(".con.")[0]
        
        if sample_fna+".fna" not in genes.keys():
            sample_fasta = f"/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/genomes/{sample_fna}.fna"
            

            try:
                with open(sample_fasta, "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        # process record as needed
                        if record.id in contigs_of_interest["col1"].values:
                            seqeunce = str(record.seq)
                            unitig_seq = unitig_string.lower()
                            seq_str = str(seqeunce).lower()
                            match_start = seq_str.find(unitig_seq)

                            if match_start != -1:
                                # Calculate start and end positions with padding, ensuring bounds
                                start = max(0, match_start - 10)
                                end = min(len(seq_str), match_start + len(unitig_seq) + 10)
                                context_seq = seq_str[start:end]
                                #add to the retrun dictionary 
                                extracted_genes_proteins[sample_fna+".fna"] = [[context_seq.upper()]]


            except (FileNotFoundError, IOError):
                continue
    
    merged_dict = {**extracted_genes_proteins, **genes}
    return merged_dict



=======
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
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480

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
<<<<<<< HEAD
    msa_output_path_gen = os.path.join(temp_dir, "msa_output_gen.fasta")
    subprocess.run([
        "clustalo", "-i", fasta_path, "-o", msa_output_path_gen, "--force", "--outfmt=fasta"
=======
    msa_output_path = os.path.join(temp_dir, "msa_output.fasta")
    subprocess.run([
        "clustalo", "-i", fasta_path, "-o", msa_output_path, "--force", "--outfmt=fasta"
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480
    ], check=True)

    #Read the aligned sequences from the output FASTA
    aligned_seqs = {}
<<<<<<< HEAD
    for record in SeqIO.parse(msa_output_path_gen, "fasta"):
=======
    for record in SeqIO.parse(msa_output_path, "fasta"):
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480
        aligned_seqs[record.id] = str(record.seq)

    #Convert the alignment to a pandas DataFrame
    alignment_genes = pd.DataFrame(
        [list(seq) for seq in aligned_seqs.values()],
        index=list(aligned_seqs.keys())
    ).T  # Transpose for SNP analysis

<<<<<<< HEAD
=======
    #Clean up: delete the temporary directory and files
    shutil.rmtree(temp_dir)


>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480
    #align proteins
    temp_dir = tempfile.mkdtemp()
    #Write each sequence to a FASTA file in the temp directory
    fasta_path = os.path.join(temp_dir, "msa_input.fasta")
    with open(fasta_path, "w") as fasta_file:
        for key in extracted_genes_proteins:
            fasta_file.write(f">{key}\n{extracted_genes_proteins[key][1]}\n")

    #Run Clustal Omega to generate the alignment (output in FASTA format)
<<<<<<< HEAD
    msa_output_path_prot = os.path.join(temp_dir, "msa_output_prot.fasta")
    subprocess.run([
        "clustalo", "-i", fasta_path, "-o", msa_output_path_prot, "--force", "--outfmt=fasta"
=======
    msa_output_path = os.path.join(temp_dir, "msa_output.fasta")
    subprocess.run([
        "clustalo", "-i", fasta_path, "-o", msa_output_path, "--force", "--outfmt=fasta"
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480
    ], check=True)

    #Read the aligned sequences from the output FASTA
    aligned_seqs = {}
<<<<<<< HEAD
    for record in SeqIO.parse(msa_output_path_prot, "fasta"):
=======
    for record in SeqIO.parse(msa_output_path, "fasta"):
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480
        aligned_seqs[record.id] = str(record.seq)

    #Convert the alignment to a pandas DataFrame
    alignmed_proteins = pd.DataFrame(
        [list(seq) for seq in aligned_seqs.values()],
        index=list(aligned_seqs.keys())
    ).T  # Transpose for SNP analysis
    #Clean up: delete the temporary directory and files
<<<<<<< HEAD


    return alignment_genes, alignmed_proteins, msa_output_path_gen, msa_output_path_prot




def extract_protein_locations_all_matches_by_gene_panaroo_name(unitig_locations_across_samples):
    '''The function takes the protein location of the proteins whcih contain the unitig and matches it to the 
    panarro annotations to see which group the genes match to if possible.
    The function returns a dictionary with all samples where the protein was douns and the respective locations 
    on all genomes. This inlcudes the genomes which include the unitig'''
    import pandas as pd
    import gffpandas.gffpandas as gffpd
    import os
    import re

    gene_panaroo = []
    for key in unitig_locations_across_samples:
        if len(unitig_locations_across_samples[key])>0:
            column = unitig_locations_across_samples[key]["attributes"].str.split(";")
            get_ID = column.values[0][0]
            gene_panaroo.append(get_ID.split("=")[1])


    #find the feature in other samples 
    gene_presence_absence = pd.read_csv("/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/genomes/gene_presence_absence.csv")

    indices = []
    for gene_name in gene_panaroo:
        indices.append(gene_presence_absence.index[gene_presence_absence.isin([gene_name]).any(axis=1)].values[0])

    # Extract the rows from gene_presence_absence
    gene_presence_absence = gene_presence_absence.drop(columns = ["Gene", "Non-unique Gene name", "Annotation"])
    selected_rows = gene_presence_absence.loc[indices]

    #get gene locations for all samples
    selected_rows = selected_rows.dropna(axis=1, how='all')
    alternative_proteins = {}
    for columns in selected_rows.columns:
        #read each gf file 
        path_to_annotations = "/Users/kvk22/Library/CloudStorage/OneDrive-UniversityofCambridge/Desktop/Arthemis/annotations/"
        gff_file = path_to_annotations + columns + ".gff"
        #check if file exists
        if os.path.exists(gff_file):
            values = selected_rows[columns].dropna().astype(str).tolist()
            pattern = re.compile(r'(' + '|'.join(re.escape(val) for val in values) + r')')
            matching_lines = []
            with open(gff_file, "r") as f:
                for line in f:
                    if pattern.search(line):
                        matching_lines.append(line.strip())
                        print(pattern)
                        print(matching_lines)
                        if len(matching_lines) == len(values):  # Stop after first match for each value
                            break
            alternative_proteins[columns] = matching_lines
    return alternative_proteins


def find_correlated_regions(extracted_genes_proteins, unitig_string):
    from collections import defaultdict


    all_kmers = set()
    for sample in extracted_genes_proteins:
        sequences = extracted_genes_proteins[sample]
        for sequence in sequences:
            k = 31
            substrings = [sequence[0][i:i+k] for i in range(len(sequence[0]) - k + 1)]
            all_kmers.update(substrings)


    # Create De Bruijn graph from k-mers
    de_bruijn_graph = defaultdict(list)
    for kmer in all_kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        de_bruijn_graph[prefix].append(suffix)

    #extend to unitigs
    unitigs = []
    visited = set()

    for node in de_bruijn_graph:
        out_edges = de_bruijn_graph[node]
        in_edges = [n for n, edges in de_bruijn_graph.items() if node in edges]
        if len(out_edges) == 1 and len(in_edges) == 1:
            continue  # skip nodes in the middle of a linear path

        for next_node in out_edges:
            if (node, next_node) in visited:
                continue
            path = [node]
            current = next_node
            while (current in de_bruijn_graph and
                len(de_bruijn_graph[current]) == 1 and
                sum([current in edges for edges in de_bruijn_graph.values()]) == 1):
                path.append(current)
                visited.add((path[-2], current))
                current = de_bruijn_graph[current][0]
            path.append(current)
            unitig = path[0] + ''.join(n[-1] for n in path[1:])
            unitigs.append(unitig)
    unitigs.append(unitig_string)
    #get unique unitigs 
    unitigs_list = []
    for u in unitigs:
        if u not in unitigs_list:
            unitigs_list.append(u)


    #create presence/absence matrix
    import pandas as pd

    samples = list(extracted_genes_proteins.keys())
    unitigs_list = list(unitigs)

    # Initialize matrix: rows=unitigs, columns=samples
    presence_absence = pd.DataFrame(0, index=unitigs_list, columns=samples)

    for sample in samples:
        seq = extracted_genes_proteins[sample][0][0]  # sequence string for the sample
        for unitig in unitigs_list:
            if unitig in seq:
                presence_absence.at[unitig, sample] = 1


    #reindex
    unitigs = pd.DataFrame(unitigs)
    target_index = list(unitigs[unitigs[0] == unitig_string].index)[0]
    unitigs = unitigs.rename(index={target_index: "target_unitig"})
    presence_absence.index = unitigs.index

    #correlate to target unitig
    presence_absence = presence_absence.T
    correlations = presence_absence.corrwith(presence_absence["target_unitig"])




    #plot the correlations
    import matplotlib.pyplot as plt
    import numpy as np

    # Convert to absolute values
    abs_correlations = correlations.abs()

    # Find the numerical index of "target_unitig"
    target_idx = list(abs_correlations.index).index("target_unitig")

    # Reindex: replace "target_unitig" with its numerical index
    new_index = list(abs_correlations.index)
    new_index[target_idx] = target_idx
    abs_correlations.index = new_index

    # Note: target_idx is the position of the target unitig
    print("Numerical index of target_unitig:", target_idx)

    # Plot
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 6))
    plt.plot(abs_correlations.index, abs_correlations.values, label="Absolute Correlation")
    plt.xlabel("Unitig Index")
    plt.ylabel("Absolute Correlation")
    plt.title("Absolute Correlation to Target Unitig")
    plt.axvline(x=target_idx, color='red', linestyle='--', label="target_unitig")
    plt.legend()
    plt.show()

    #prind which unitigs are hgihly correlated:
    unitigs_correlations = pd.merge(pd.DataFrame(correlations), unitigs, right_index=True, left_index=True)
    unitigs_correlations.columns = ["correlation", "unitig"]
    return unitigs_correlations


def run_ncbi_annotations(genes, number_of_hits_to_run = 10):
    from Bio.Blast import NCBIWWW
    import pandas as pd
    from io import StringIO
    from Bio.Blast import NCBIXML

    all_hits = pd.DataFrame()  

    first_30_genes = dict(list(genes.items())[:number_of_hits_to_run])
    for key in first_30_genes:
        protein = genes[key][0][1]

        # Example: run BLASTp for each protein sequence in list_of_entries
        blast_results = []
        # Run BLASTp (this queries NCBI's online BLAST server)
        result_handle = NCBIWWW.qblast("blastp", "nr", protein)
        blast_results.append(result_handle.read())
        result_handle.close()


        # Parse the BLAST result (assuming blast_results[0] contains the XML string)
        blast_record = NCBIXML.read(StringIO(blast_results[0]))

        # Extract top 5 hits with coverage and identity
        hits = []
        for alignment in blast_record.alignments[:5]:
            hsp = alignment.hsps[0]
            coverage = hsp.align_length / alignment.length
            identity = (hsp.identities / hsp.align_length) * 100
            hits.append({
                "hit_id": alignment.hit_id,
                "hit_def": alignment.hit_def,
                "length": alignment.length,
                "e_value": hsp.expect,
                "score": hsp.score,
                "coverage": coverage,
                "identity": identity
            })

        # Put results in a DataFrame
        df_hits = pd.DataFrame(hits)
        all_hits = pd.concat([all_hits, df_hits], ignore_index=True)
    return all_hits
=======
    shutil.rmtree(temp_dir)


    return alignment_genes, alignmed_proteins
>>>>>>> 15d7af32ab701fe8e54d73b9f5dcf17bac926480
