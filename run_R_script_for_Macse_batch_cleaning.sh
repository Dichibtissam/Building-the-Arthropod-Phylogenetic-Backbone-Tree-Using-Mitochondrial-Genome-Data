#!/bin/bash

# Définir le chemin vers le script R
R_SCRIPT_PATH="/home/idich/internship/ibtissam/alignement/Scripts/Diplostraca/R_script_for_Macse_batch_cleaning.R"  # Remplacez "your_r_script.R" par le nom exact de votre script R

# Exécutez le script R avec Rscript
Rscript "$R_SCRIPT_PATH"

# Vérifiez si le script R a bien été exécuté
if [ $? -eq 0 ]; then
    echo "✔️ Script R exécuté avec succès."
else
    echo "❌ L'exécution du script R a échoué."
fi

