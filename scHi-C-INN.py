import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
import argparse

def get_chromosomes(genome_type):
    if genome_type == 'hg19' or genome_type == 'hg38':
        return ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
    elif genome_type == 'mm9' or genome_type == 'mm10':
        return ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY"]
    else:
        raise ValueError("Unsupported genome type. Choose from 'hg19', 'hg38', 'mm9', 'mm10'.")

def process_chromosome(ch, base_dir, outputs, inputs, correlation_dir, chrlen_dic, cell_num, num_neighbors):
    output_folder = os.path.join(base_dir, outputs, ch)
    os.makedirs(output_folder, exist_ok=True)
    hic_who = pd.read_csv(os.path.join(correlation_dir, f'cell_inter_{ch}.csv'))
    
    for i in range(1, cell_num+1, 1):
        size = chrlen_dic[ch]
        f = np.loadtxt(os.path.join(base_dir, f'{inputs}/{ch}/{ch}_cell{i}_{inputs}.txt'))
        f = f.reshape((size,size))
        xs, ys = np.where(f == 0)
        res = hic_who[hic_who.cell1 == i]['cell2'].values[:num_neighbors]

        if len(res) > 0:
            neighbor_matrices = [np.loadtxt(os.path.join(base_dir, f'{inputs}/{ch}/{ch}_cell{res[j]}_{inputs}.txt')).reshape((size, size)) for j in range(len(res))]

            for x, y in zip(xs, ys):
                if np.abs(x - y) < 10:
                    f[x, y] = np.mean([matrix[x, y] for matrix in neighbor_matrices])
        
        np.savetxt(os.path.join(output_folder, f'{ch}_cell{i}_{outputs}.txt'), f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process chromosome data for different genomes")
    parser.add_argument("--base_dir", required=True, help="Base directory")
    parser.add_argument("--cell_num", type=int, required=True, help="Cell number")
    parser.add_argument("--genome_type", required=True, help="Genome type (hg19, hg38, mm9, mm10)")
    parser.add_argument("--outputs", required=True, help="Output type")
    parser.add_argument("--inputs", required=True, help="Input type")
    parser.add_argument("--correlation_dir", required=True, help="Directory with correlation coefficient results")
    parser.add_argument("--num_neighbors", type=int, required=True, help="Number of neighbors for imputation")

    args = parser.parse_args()

    chrom = get_chromosomes(args.genome_type)
    
    # Define chromosome lengths based on the genome type
    # For demonstration, using simplified lengths. Update with accurate lengths for each genome type.
    if args.genome_type in ['hg19', 'hg38']:
        chr_len = [250, 244, 199, 192, 181, 172, 160, 147, 142, 136, 136, 134, 116, 108, 103, 91, 82, 79, 60, 64, 49, 52, 156]
    elif args.genome_type in ['mm9', 'mm10']:
        chr_len = [195, 182, 160, 157, 152, 149, 145, 129, 124, 130, 122, 120, 108, 125, 103, 98, 91, 88, 61, 48, 51]  # Example lengths
    else:
        raise ValueError("Unsupported genome type.")

    chrlen_dic = {chrom[i]: chr_len[i] for i in range(len(chrom))}

    with Pool() as pool:
        pool.map(partial(process_chromosome, base_dir=args.base_dir, outputs=args.outputs, inputs=args.inputs, correlation_dir=args.correlation_dir, chrlen_dic=chrlen_dic, cell_num=args.cell_num, num_neighbors=args.num_neighbors), chrom)
