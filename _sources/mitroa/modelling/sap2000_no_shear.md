
# Επίλυση πλαισίου με το SAP2000 (αγνοώντας το έργο των τεμνουσών)

Θα εξεταστεί το ίδιο παράδειγμα, με μόνη διαφορά ότι θα αγνοηθεί το έργο των τεμνουσών.

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_foreas.png
---
width: 500px
align: center
---
```

## Διατομές δομικών στοιχείων

### Modification factors

Για να γίνει αυτό, η μόνη διαφοροποίηση στην προσομοίωση του παραδείγματος σχετίζεται με τα modification factors στις ιδιότητες των δομικών στοιχείων.

Μηδενίζονται τα πεδία που σχετίζονται με τις τέμνουσες δυνάμες (Shear Area). Αυτό πρέπει να γίνει και στις διατομές όλων των δομικών στοιχείων (δοκού και στύλων).

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap2_mod_factors.png
---
width: 400px
align: center
name: pi_sap2_mod_factors
---
```

## Αποτελέσματα ανάλυσης

Τρέχοντας την ανάλυση παρατηρείται ότι τα διαγράμματα των εντατικών μεγεθών αλλάζουν και πλέον ταυτίζονται όμως με τα αποτελέσματα της επίλυσης με την κλασική στατική (με το χέρι).

### Διάγραμμα ροπών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap2_results_moments.png
---
width: 650px
align: center
name: pi_sap2_results_moments
---
```

### Διάγραμμα τεμνουσών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap2_results_shear.png
---
width: 650px
align: center
name: pi_sap2_results_shear
---
```

### Διάγραμμα αξονικών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap2_results_axial.png
---
width: 650px
align: center
name: pi_sap2_results_axial
---
```

### Παρατηρήσεις

Παρατηρείται ότι τα παραπάνω διαγράμματα έχουν τιμές που **ταυτίζονται** με τα αποτελέσματα της επίλυσης με την κλασική στατική (με το χέρι).