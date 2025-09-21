#!/usr/bin/env python 

import re
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import sys

# sequência de referência do promotor (com ambiguidades IUPAC)
promoter_ref = "TSSCAATKTKACTSKCATAWGTTTWWCGCCAC"

# sequência de aminoácidos de referência
protein_ref = "MSEQHRYPTWDFVGALNKIASALAEGMEVAGANARVLYWSPQMNRVCAVSKLKLHVFDSL"

# standard genetic code: https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
genetic_code = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# converter ambiguidades IUPAC para expressão regular
def iupac_to_regex(iupac_seq):
    iupac_dict = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]',
        'W': '[AT]', 'K': '[GT]', 'M': '[AC]',
        'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]',
        'V': '[ACG]', 'N': '[ACGT]'
    }
    regex = ''
    for base in iupac_seq.upper():
        regex += iupac_dict.get(base, base)
    return regex

# traduzir sequência de DNA em proteína
def translate_dna(sequence):
    protein = ""
    for i in range(0, len(sequence) - (len(sequence) % 3), 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:  # verifica se o códon está completo
            aa = genetic_code.get(codon, 'X')
            if aa == '*':  # stop codon
                break
            protein += aa
    return protein

# calcula similaridade entre duas sequências (proteínas)
def calculate_similarity(seq1, seq2):
    min_length = min(len(seq1), len(seq2))
    if min_length == 0:
        return 0.0
    
    # conta aminoácidos idênticos nas posições correspondentes
    matches = sum(1 for i in range(min_length) if seq1[i] == seq2[i])
    
    # retorna a proporção de correspondências
    similarity = matches / min_length
    return similarity

# encontra o início do gene (sequência codificante)
def find_coding_sequence(genome, protein_ref):
    best_start_pos = -1
    best_similarity = 0.0
    
    # procurar apenas os start codons (ATG)
    start_codon = "ATG"
    start_positions = [m.start() for m in re.finditer(start_codon, genome)]
    
    for start_pos in start_positions:
        # extraím a sequência codificante a partir do start codon
        coding_seq = genome[start_pos:]
        
        # traduz a sequência para proteína
        translated = translate_dna(coding_seq)
        
        # calcula a similaridade com a proteína de referência
        # OBS: aqui consideramos apenas o comprimento da proteína traduzida para comparação
        max_ref_len = min(len(translated), len(protein_ref))
        similarity = calculate_similarity(translated[:max_ref_len], protein_ref[:max_ref_len])
        
        # se a similaridade for maior que um limiar e tiver comprimento mínimo significativo
        # (pelo menos 10 aminoácidos para evitar falsos positivos)
        if similarity > best_similarity and len(translated) >= 10:
            best_similarity = similarity
            best_start_pos = start_pos
            
    # se encontrarmos uma sequência com boa similaridade (>70%)
    if best_similarity > 0.7:
        ##print(f"\nMelhor correspondência:")
        ##print(f"Posição inicial: {best_start_pos}")
        ##print(f"Similaridade: {best_similarity:.2f}")
        return best_start_pos
 
    return -1  # return -1 caso não encontre o início do gene

# determina a atividade transcricional com base na distância
# usamos como referência a tabela de distância do desafio
def classify_transcription_activity(distance):
    if distance is None or distance == "N/A":
        return "desconhecida"
    distance = int(distance)
    if distance > 401:
        return "baixa"
    elif 201 <= distance <= 400:
        return "moderada"
    elif 101 <= distance <= 200:
        return "alta"
    else:  # distance <= 100
        return "nula"

# analisar as cepas
def analyze_strains(sequences):
    results = []
    promoter_pattern = re.compile(iupac_to_regex(promoter_ref))
    
    total_sequences = len(sequences)
    print(f"Analisando {total_sequences} sequências...")
    
    for i, (seq_id, genome) in enumerate(sequences.items(), 1):
        if i % 50 == 0 or i == total_sequences:  # mostrar progresso a cada 50 cepas
            print(f"Processando cepa {i} de {total_sequences} ({i/total_sequences*100:.1f}%)...")
        
        # procurar pelo promotor
        promoter_match = promoter_pattern.search(genome)
        promoter_position = promoter_match.start() if promoter_match else -1
        
        # procurar pelo início do gene
        gene_start = find_coding_sequence(genome, protein_ref)
        
        if promoter_position >= 0 and gene_start >= 0:
            # calcular a distância entre o fim do promotor e o início do gene
            promoter_end = promoter_position + len(promoter_ref)
            distance = gene_start - promoter_end
            
            # classificar a atividade transcricional
            activity = classify_transcription_activity(distance)
            
            results.append({
                "Strain": seq_id,
                "Promoter_Position": promoter_position,
                "Gene_Start": gene_start,
                "Distance": distance,
                "Transcription_Activity": activity
            })
        else:
            # caso não encontrar o promotor ou o gene
            results.append({
                "Strain": seq_id,
                "Promoter_Position": promoter_position,
                "Gene_Start": gene_start,
                "Distance": "N/A",
                "Transcription_Activity": "Desconhecida"
            })
    
    return pd.DataFrame(results)

# ler arquivo FASTA usando Biopython
def read_fasta_file(fasta_file):
    sequences = {}
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = str(record.seq)
        return sequences
    except Exception as e:
        print(f"Erro ao ler o arquivo FASTA: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) < 2:
        print("Usage: python screening.py dataset_exercise.fasta")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    
    print(f"Lendo sequências do arquivo: {fasta_file}")
    
    sequences = read_fasta_file(fasta_file)
    print(f"Total de sequências lidas: {len(sequences)}")
    
    print("Analisando atividade transcricional...")
    results = analyze_strains(sequences)
    
    # exibir resultados
    print("\nResultados da análise (primeiras 5 linhas):")
    print(results.head())
    
    # exibir estatísticas de atividade transcricional
    activity_counts = results['Transcription_Activity'].value_counts()
    print("\nDistribuição de atividade transcricional:")
    print(activity_counts)
    
    # calcular estatísticas adicionais
    valid_distances = [d for d in results['Distance'] if d != "N/A"]
    if valid_distances:
        print(f"\nEstatísticas de distância promotor-gene:")
        print(f"Média: {sum(valid_distances)/len(valid_distances):.2f} nucleotídeos")
        print(f"Mínima: {min(valid_distances)} nucleotídeos")
        print(f"Máxima: {max(valid_distances)} nucleotídeos")
    
    # save results 
    output_file = "transcription_activity_results.csv"
    results.to_csv(output_file, index=False)
    print(f"\nResultados salvos em '{output_file}'")

    # salvando um relatório resumido
    # apenas para curiosidade
    with open("transcription_report.txt", "w") as f:
        f.write("RELATÓRIO DE ATIVIDADE TRANSCRICIONAL\n")
        f.write("===================================\n\n")
        f.write(f"Total de cepas analisadas: {len(results)}\n\n")
        f.write("Distribuição por atividade transcricional:\n")
        for activity, count in activity_counts.items():
            percentage = (count / len(results)) * 100
            f.write(f"- {activity}: {count} cepas ({percentage:.1f}%)\n")
        
        if valid_distances:
            f.write("\nEstatísticas de distância promotor-gene:\n")
            f.write(f"- Média: {sum(valid_distances)/len(valid_distances):.2f} nucleotídeos\n")
            f.write(f"- Mínima: {min(valid_distances)} nucleotídeos\n")
            f.write(f"- Máxima: {max(valid_distances)} nucleotídeos\n")
    
    print(f"Relatório resumido salvo em 'transcription_report.txt'")

if __name__ == "__main__":
    main()
