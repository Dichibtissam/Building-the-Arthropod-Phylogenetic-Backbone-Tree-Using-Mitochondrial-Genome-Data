#--------------------------------#
#            Étape 1             #
#   Transfert du fichier local   #
#--------------------------------#

# Transférer le fichier ncbi_dataset.tsv du PC local vers le serveur distant asellus
scp /home/wsl/home/wsl/stage/essey_3180182025/finale/ncbi_dataset.tsv idich@umr5023-asellusc1.univ-lyon1.fr:/home/idich/internship/ibtissam/metadata_study/


#--------------------------------#
#            Étape 2             #
#  Vérification du fichier       #
#--------------------------------#

# Afficher les premières lignes du fichier pour s'assurer qu'il a bien été transféré
head ncbi_dataset.tsv

#--------------------------------#
#            Étape 3             #
#  Extraction des TaxIDs         #
#--------------------------------#

# Extraire la 6e colonne (TaxID) du fichier et enregistrer les résultats dans taxids.txt                                         
awk -F'\t' 'NR>1 {print $6}' ncbi_dataset.tsv > taxids.txt

# Afficher les premières lignes du fichier taxid_list.txt pour vérifier le résultat                                                  
head taxids.txt

# Compte le nombre de lignes dans taxid_list.txt                                                     
wc -l taxids.txt

#--------------------------------#
#            Étape 4             #
#  Installation de TaxonKit      #
#--------------------------------#

                                                                                                                                     
# Vérification si le fichier taxids.txt existe
ls -l taxids.txt
                                                                                                                                     
# Télécharger la base de données taxonomique de TaxonKit
wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# Décompresser le fichier téléchargé
tar -xvzf taxdump.tar.gz


#--------------------------------#
#            Étape 5             #
#   Attribution des Taxonomies   #
#--------------------------------#

# Obtenir les lignées taxonomiques pour les TaxIDs en utilisant TaxonKit
cat taxids.txt | taxonkit reformat -I 1 -f "{p};{c};{o};{f};{g};{s}" -F > out_reformat_taxids.txt

#--------------------------------#
#            Étape 6             #
#   Traitement des données       #
#--------------------------------#

# Charger les bibliothèques nécessaires
python
import pandas as pd

# Charger le fichier contenant les taxonomies
df = pd.read_csv("out_reformat_taxids.txt", sep="\t", header=None, names=["taxid", "lineage"])

# Séparer la colonne 'lineage' en plusieurs colonnes
df_split = df["lineage"].str.split(";", expand=True)

# Renommer les colonnes 
df_split.columns = ["phylum", "class", "order", "family", "genus", "species"]

# Ajouter la colonne 'taxid' au dataframe final
df_final = pd.concat([df["taxid"], df_split], axis=1)

# Sauvegarder le fichier sous forme de CSV avec des en-têtes claires
df_final.to_csv("out_reformat_taxids_with_headers.csv", index=False)

# Afficher les premières lignes pour vérification
print(df_final.head())

# Quitter l'exécution du script Python et retourner au terminal
quit()

# Si l'on souhaite compter uniquement les lignes de données (sans l'en-tête)
tail -n +2 out_reformat_taxids_with_headers.csv | wc -l



#--------------------------------#
#            Etape 7             #
#    Choix des out-groups        #
#--------------------------------#


# Charger les bibliothèques nécessaires
library(dplyr)
library(readr)

# Lire les fichiers
outgroup_df <- read.csv("out_reformat_taxids_with_headers.csv", header = TRUE, stringsAsFactors = FALSE)
ncbi_data <- read_tsv("ncbi_dataset.tsv", col_types = cols())

# Convertir 'taxid' en numérique pour une meilleure gestion
outgroup_df$taxid <- as.numeric(outgroup_df$taxid)

# Fusionner avec les accession numbers (colonne 'Organism taxid' dans ncbi_data)
ncbi_data$`Organism taxid` <- as.numeric(ncbi_data$`Organism taxid`)  # Assurer que 'taxid' est numérique
merged_df <- left_join(outgroup_df, ncbi_data, by = c("taxid" = "Organism taxid"), multiple = "all")

# Fonction mise à jour pour sélectionner les outgroups avec les accession numbers et taxids
select_outgroups <- function(df) {
  outgroup_list <- list()
  
  # Sélectionner les ordres ayant au moins 5 individus
  order_counts <- df %>% 
    group_by(order) %>% 
    tally() %>% 
    filter(n >= 5) %>% 
    pull(order)
  
  df_filtered <- df %>% filter(order %in% order_counts)
  
  for (order_name in unique(df_filtered$order)) {
    
    # Identifier les classes de l'ingroup
    ingroup_class <- df_filtered %>% 
      filter(order == order_name & !is.na(class) & class != "") %>% 
      pull(class) %>% 
      unique()
    
    if (length(ingroup_class) == 0) next
    
    # Sélectionner les outgroups pour les classes de l'ingroup
    outgroup_candidates <- df_filtered %>% 
      filter(class %in% ingroup_class, order != order_name)
    
    if (nrow(outgroup_candidates) == 0) next
    
    # Sélectionner aléatoirement un ordre parmi les candidats
    selected_order <- sample(unique(outgroup_candidates$order), 1)
    
    # Sélectionner les génomes des outgroups
    selected_outgroup_genomes <- outgroup_candidates %>% 
      filter(order == selected_order)
    
    # Vérifier si le nombre de génomes est suffisant
    if (nrow(selected_outgroup_genomes) < 5) next
    
    # Sélectionner aléatoirement 5 génomes
    selected_outgroup_genomes <- selected_outgroup_genomes %>% sample_n(5)
    
    Sys.sleep(1)  # Pause d'une seconde entre chaque sélection
    
    # Compter le nombre de génomes dans l'ingroup
    ingroup_genomes <- df_filtered %>% 
      filter(order == order_name)
    
    ingroup_genome_count <- nrow(ingroup_genomes)
    
    # Accession numbers et taxids de l'ingroup
    ingroup_accessions <- paste(ingroup_genomes$`GenBank accession`, collapse = ", ")
    ingroup_taxids <- paste(ingroup_genomes$taxid, collapse = ", ")
    
    # Accession numbers et taxids des outgroups
    outgroup_accessions <- paste(selected_outgroup_genomes$`GenBank accession`, collapse = ", ")
    outgroup_taxids <- paste(selected_outgroup_genomes$taxid, collapse = ", ")
    
    # Ajouter à la liste
    outgroup_list[[order_name]] <- list(
      order = order_name, 
      class = ingroup_class, 
      outgroup_order = selected_order, 
      outgroup_accessions = outgroup_accessions,
      ingroup_accessions = ingroup_accessions,
      ingroup_taxids = ingroup_taxids,
      outgroup_taxids = outgroup_taxids,
      ingroup_genome_count = ingroup_genome_count
    )
  }
  
  # Créer un DataFrame à partir de la liste
  outgroup_df <- bind_rows(lapply(outgroup_list, function(x) {
    data.frame(order = x$order,
               class = x$class,
               outgroup_order = x$outgroup_order,
               outgroup_accessions = x$outgroup_accessions,
               ingroup_accessions = x$ingroup_accessions,
./script_alignment_Macse.sh               ingroup_taxids = x$ingroup_taxids,
               outgroup_taxids = x$outgroup_taxids,
               ingroup_genome_count = x$ingroup_genome_count)
  }))
  
  return(outgroup_df)
}

# Application de la fonction
outgroup_selection <- select_outgroups(merged_df)

# Export des résultats
write.csv(outgroup_selection, "outgroup_selection_with_accessions_and_taxids.csv", row.names = FALSE)

# Affichage des résultats
print(outgroup_selection)

    
 
