import os, json
import pandas as pd
import numpy as np
from multiprocessing import Pool
import argparse

def process_chromosome(c, cell_num, label_df, input_dir, output_dir):
    input_file = os.path.join(input_dir, f'{c}.txt')
    df = pd.read_csv(input_file, sep=',')
    all_tgt = pd.DataFrame()

    for i in range(1, cell_num+1):
        c_list = []
        j = i + 1
        if i == 1:
            while j <= cell_num:
                c_list.append(str(i)+'-'+str(j))
                j += 1
        else:
            for m in range(1, i):
                c_list.append(str(m)+'-'+str(i))
            while j <= cell_num:
                c_list.append(str(i)+'-'+str(j))
                j += 1
        selected_df = df.loc[df.index.isin(c_list)] # rows for specific cell
        cell1 = label_df[label_df.cell==i]['group'].values[0]
        cell2 = label_df[label_df.group==cell1]['cell'].values
        sameG_idx1 = map(lambda x,y: str(x)+'-'+str(y), [i]*len(cell2), cell2)
        sameG_idx2 = map(lambda x,y: str(x)+'-'+str(y), cell2, [i]*len(cell2))
        sameG_list1 = list(sameG_idx1)
        sameG_list2 = list(sameG_idx2)
        sameG_list = sameG_list1 + sameG_list2
        sameG_df = selected_df.loc[selected_df.index.isin(sameG_list)]
        # sort and get the max five cells (absolute value)
        maxs = np.sort(sameG_df['weighted_pearson'])
        maxs = maxs[~np.isnan(maxs)]
        maxs_slt = maxs[-5:len(maxs)]
        # concatenate them to a new dataframe
        final_df = pd.DataFrame()
        for k in maxs_slt:
            final_df = pd.concat([final_df, sameG_df[sameG_df['weighted_pearson']==k]])
        ind1 = list([i]*len(final_df))
        ind2 = []
        for d in final_df.index:
            two = d.split('-')
            if int(two[1])==i:
                ind2.append(int(two[0]))
            else:
                ind2.append(int(two[1]))
        try:
            sep_df = pd.DataFrame({'cell1':ind1, 'cell2':ind2, 'hicrep':final_df['weighted_pearson']})
            all_tgt = pd.concat([all_tgt, sep_df])
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f'cell_inter_{c}.csv')
            all_tgt.to_csv(output_file, index=False)
        except:
            pass

def get_chromosomes(genome_type):
    if genome_type == 'hg19' or genome_type == 'hg38':
        return ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
    elif genome_type == 'mm9' or genome_type == 'mm10':
        return ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY"]
    else:
        raise ValueError("Unsupported genome type. Choose from 'hg19', 'hg38', 'mm9', 'mm10'.")

def main():
    parser = argparse.ArgumentParser(description="Process chromosome data for different genomes")
    parser.add_argument("--input_dir", required=True, help="Input directory")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--cell_num", type=int, required=True, help="Cell number")
    parser.add_argument("--label_dir", required=True, help="Cell type label directory")
    parser.add_argument("--genome_type", required=True, help="Genome type (hg19, hg38, mm9, mm10)")

    args = parser.parse_args()

    chrom = get_chromosomes(args.genome_type)

    f = open(os.path.join(args.label_dir, 'label_info.json'))
    label = json.load(f)    
    label = label['cell type']
    cells = [i for i in range(1, args.cell_num + 1)]
    label_df = pd.DataFrame({'cell': cells, 'group': label})

    with Pool() as pool:
        pool.starmap(process_chromosome, [(c, args.cell_num, label_df, args.input_dir, args.output_dir) for c in chrom])

if __name__ == "__main__":
    main()
