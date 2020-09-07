% Pre-process gene names
function [processed_gene_name] = gene_name_preprocess(gene_name, ncbi_name)

for row = 1:size(gene_name, 1)
    if contains(ncbi_name(row, :), 'deleted_in_esophageal_cancer_1') == 1
        gene_name{row, :} = 'DEC1';
    elseif contains(ncbi_name(row, :), 'Homo_sapiens_15_kDa_selenoprotein__SEP15__transcript_variant_2__mRNA._Source_RefSeq_mRNA_Acc_NM_203341') == 1
        gene_name{row, :} = 'SEP15';
    elseif strcmp(ncbi_name{row, :}, 'septin_1') == 1
        gene_name{row, :} = 'SEPT1';
    elseif strcmp(ncbi_name{row, :}, 'septin_2') == 1
        gene_name{row, :} = 'SEPT2';
    elseif strcmp(ncbi_name{row, :}, 'septin_3') == 1
        gene_name{row, :} = 'SEPT3';
    elseif strcmp(ncbi_name{row, :}, 'septin_4') == 1
        gene_name{row, :} = 'SEPT4';       
    elseif strcmp(ncbi_name{row, :}, 'septin_5') == 1
        gene_name{row, :} = 'SEPT5';        
    elseif strcmp(ncbi_name{row, :}, 'septin_6') == 1
        gene_name{row, :} = 'SEPT6';    
    elseif strcmp(ncbi_name{row, :}, 'septin_7') == 1
        gene_name{row, :} = 'SEPT7';
    elseif strcmp(ncbi_name{row, :}, 'septin_8') == 1
        gene_name{row, :} = 'SEPT8';        
    elseif strcmp(ncbi_name{row, :}, 'septin_9') == 1
        gene_name{row, :} = 'SEPT9';        
    elseif strcmp(ncbi_name{row, :}, 'septin_10') == 1
        gene_name{row, :} = 'SEPT10';
    elseif strcmp(ncbi_name{row, :}, 'septin_11') == 1
        gene_name{row, :} = 'SEPT11';        
    elseif strcmp(ncbi_name{row, :}, 'septin_12') == 1
        gene_name{row, :} = 'SEPT12';      
    elseif strcmp(ncbi_name{row, :}, 'septin_13') == 1
        gene_name{row, :} = 'SEPT13';
    elseif strcmp(ncbi_name{row, :}, 'septin_14') == 1
        gene_name{row, :} = 'SEPT14';
    elseif contains(ncbi_name(row, :), 'mitochondrial_amidoxime_reducing_component_1') == 1
        gene_name{row, :} = 'MARC1';
    elseif contains(ncbi_name(row, :), 'mitochondrial_amidoxime_reducing_component_2') == 1
        gene_name{row, :} = 'MARC2';    
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__1__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH1';
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__2__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH2';
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__3__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH3';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__4__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH4';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__5') == 1
        gene_name{row, :} = 'MARCH5';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__6__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH6';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__7__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH7';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__8__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH8';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__9__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH9';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__10__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH10';        
    elseif contains(ncbi_name(row, :), 'membrane-associated_ring_finger__C3HC4__11__E3_ubiquitin_protein_ligase') == 1
        gene_name{row, :} = 'MARCH11';
    end
end

processed_gene_name = cellstr(string(gene_name));

end