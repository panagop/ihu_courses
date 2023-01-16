
# Επίλυση πλαισίου με το SAP2000 (default)

## Γεωμετρία φορέα

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_foreas.png
---
width: 500px
align: center
name: pi_foreas
---
```


### Δημιουργία νέου φορέα

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_newmodel.png
---
width: 450px
align: center
name: pi_sap_newmodel
---
```

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_2dframe.png
---
width: 450px
align: center
name: pi_sap_2dframe
---
```

### Υλικά

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_material.png
---
width: 350px
align: center
name: pi_sap_material
---
```

### Συνθήκες στήριξης

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_resrtaints.png
---
width: 600px
align: center
name: pi_sap_resrtaints
---
```

## Διατομές δομικών στοιχείων

### Δοκός

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_beam.png
---
width: 400px
align: center
name: pi_sap_beam
---
```

### Στύλοι

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_col.png
---
width: 400px
align: center
name: pi_sap_col
---
```

Προσοχή ώστε να οριστούν σωστά οι πλευρές (μικρή/μεγάλη) του υποστυλώματος (συσχέτιση τοπικών και καθολικών αξόνων)

### Modification factors

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_mod_factors.png
---
width: 400px
align: center
name: pi_sap_mod_factors
---
```

Στα modification factors των διατομών αφήνουμε παντού μονάδες, εκτός από τη μάζα και το βάρος που δε θέλουμε να υπολογιστεί αυτόματα από τις διαστάσεις των στοιχείων (άλλωστε είχαμε ορίσει μηδενικό ειδικό βάρος στο υλικό)

### Ανάθεση των διατομών στις ράβδους του φορέα

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_frame_sections.png
---
width: 500px
align: center
name: pi_sap_frame_sections
---
```

## Φορτία

### Load patterns/cases

Δίνεται μόνο ένα load pattern/case καθώς θα γίνει απλή επαλληλία της οριζόντιας φόρτισης και του ομοιόμορφου φορτίου της δοκούς, χωρίς κάποιο πρόσθετο συνδυασμό (με διαφορετικούς συντελεστές) μεταξύ τους.

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_load_patterns.png
---
width: 500px
align: center
name: pi_sap_load_patterns
---
```

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_load_cases.png
---
width: 500px
align: center
name: pi_sap_load_cases
---
```

### Εισαγωγή επικόμβιων φορτίων

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_nodal_loads.png
---
width: 500px
align: center
name: pi_sap_nodal_loads
---
```

### Εισαγωγή ομοιόμορφων φορτίων δοκού

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_frame_loads.png
---
width: 600px
align: center
name: pi_sap_frame_loads
---
```

### Φορέας με την εξωτερική φόρτιση

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_loads.png
---
width: 500px
align: center
name: pi_sap_loads
---
```

## Αποτελέσματα ανάλυσης

### Διάγραμμα ροπών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_results_moments.png
---
width: 650px
align: center
name: pi_sap_results_moments
---
```

### Διάγραμμα τεμνουσών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_results_shear.png
---
width: 650px
align: center
name: pi_sap_results_shear
---
```

### Διάγραμμα αξονικών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_results_axial.png
---
width: 650px
align: center
name: pi_sap_results_axial
---
```

### Παρατηρήσεις

Παρατηρείται ότι τα παραπάνω διαγράμματα έχουν τιμές που είναι κοντά, δεν ταυτίζονται όμως με τα αποτελέσματα της επίλυσης με την κλασική στατική (με το χέρι).