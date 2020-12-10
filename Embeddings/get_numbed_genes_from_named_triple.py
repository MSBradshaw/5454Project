file_path = '/Users/michael/PycharmProjects/VariantRankingPheKnowLator/Data/PheKnowLator_Subclass_InverseRelations_NotClosed_NoOWL_Triples_Identifiers_11_24_2020.txt'
name_to_num = {}
count = 0
for line in open(file_path,'r'):
    row = line.strip().split('\t')
    for item in row:
        if item not in name_to_num:
            name_to_num[item] = count
            count += 1

gene_nums = [name_to_num[key] for key in name_to_num.keys() if 'https://www.ncbi.nlm.nih.gov/gene/' in key ]
with open('gene_numbers.txt','w') as out:
    for num in gene_nums:
        out.write(str(num)+ '\n')

with open('genes_only_embedding.emb','w') as out:
    for line in open('Embeddings/PheKnowLator_Subclass_InverseRelations_NotClosed_NoOWL_Integers-NOT-TRIPLE_node2vec_Embeddings.emb','r'):
        row = line.strip().split(' ')
        if int(row[0]) in gene_nums:
            print('Num!')
            out.write(line)