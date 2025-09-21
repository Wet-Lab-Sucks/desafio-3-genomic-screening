# Desafio 3 - Screening de Atividade Transcricional em E. coli

## Descrição do Desafio

Neste desafio de bioinformática genômica, a tarefa era ajudar pesquisadores a identificar a melhor cepa de `E. coli` para a produção de biocombustíveis.  

>[!NOTE]
>No contexo deste desafio, genes exógenos foram inseridos em cepas de `E. coli` para otimizar a produção de biocombustíveis.  

Os pesquisadores observaram que a transcrição de uma proteína endógena `"X"` afeta **negativamente** o desempenho da produção de biocombustível. A expressão dessa proteína é, por sua vez, regulada pela `distância entre a sua região promotora e o gene codificador`.

O desafio consistiu em auxiliar os pesquisadores a realizar um screening de perfil transcricional para `500 cepas`, a fim de identificar a mais promissora para a produção de biocombustível. 

O objetivo final era classificar a atividade transcricional de cada cepa (`baixa`, `moderada`, `alta` ou `nula`) com base na distância de nucleotídeos entre o promotor e o gene da proteína `"X"`.

## Dados

Para esta análise, foram fornecidos fragmentos genômicos de 500 cepas.  

A tarefa foi realizada utilizando as seguintes sequências de referência e critérios de classificação:
* Sequência de Referência do Promotor: `TSSCAATKTKACTSKCATAWGTTTWWCGCCAC`
* Sequência de Referência da Proteína `"X"`: `MSEQHRYPTWDFVGALNKIASALAEGMEVAGANARVLYWSPQMNRVCAVSKLKLHVFDSL`
* Código Genético: [Standard](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables)

A classificação da atividade transcricional foi baseada nas seguintes faixas de distância:

| **Distância do Promotor ao Gene**	| **Nível de Atividade Transcricional** |
| ----------------------------- | --------------------------------- |
| Acima de 401 nucleotídeos	    | Baixa                             |
| 201 a 400 nucleotídeos	    | Moderada                          |
| 101 a 200 nucleotídeos	    | Alta                              |
| Até 100 nucleotídeos	        | Nula                              |


## Nossa estratégia

Nossa abordagem para este problema de genômica comparativa foi desenvolver um [pipeline](src/screening.py) que consistiu em: 
* Processar os 500 fragmentos genômicos, procurando pelas sequências de referência do promotor e da proteína `"X"` para calcular a distância entre elas;
* Com base nessa distância e na tabela de classificação fornecida pela LBB, nós classificamos a atividade transcricional de cada cepa como `baixa`, `moderada`, `alta` ou `nula`;
* Geramos um relatório final para a submissão:

```bash
Sequência   Nível
seq1    alta
seq2    baixa
seq3    baixa
seq4    nula
seq5    moderada
```

## Resultados

Essa abordagem se mostrou precisa e eficiente, atingindo um excelente resultado no desafio.  
Tivemos a pontuação máxima neste desafio: 10/10!
