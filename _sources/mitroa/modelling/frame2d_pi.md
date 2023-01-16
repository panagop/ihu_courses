# Επίλυση πλαισίου με το Frame2d

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_foreas.png
---
width: 500px
align: center
---
```

To βασιμό μενού για την εισαγωγή των δεδομένων στο Frame2d είναι το `MODEL`

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_model_menu.png
---
width: 300px
align: center
name: pi_frame2d_model_menu
---
```

## Μονάδες

Επιλέγεται το σύστημα των μονάδων. Στο παράδειγμα θα χρησιμοποιηθεί το `Default Metric`

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_units.png
---
width: 600px
align: center
name: pi_frame2d_units
---
```

## Ιδιότητες δομικών στοιχείων

### Υλικό

Εισάγεται το νέο υλικό και παρατηρείται ότι η μόνη ιδιότητα που μας ενδιαφέρει είναι το μέτρο ελαστικότητας (δεν υπάρχει το μέτρο διάτμησης, ο λόγος του Poisson κτλ.).

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_materials.png
---
width: 450px
align: center
name: pi_frame2d_materials
---
```

### Διατομές

#### Δοκός

Υπάρχει η δυνατότητα για απ'ευθείας εισαγωγή των τιμώς της ροπής αδράνειας $I$ και του εμβαδού $A$ της διατομής, είτε να χρησιμοποιηθεί το `Calculator`, οπότε να γίνει η εισαγωγή των διαστάσεων της διατομής. 
Προσοχή στις μονάδες και στη συσχέτιση των τοπικών αξόνων αναφοράς με το καθολικό σύστημα συντεταγμένων του φορέα. Στο παράδειγμα $a$ και $b$ είναι το πλάτος και το ύψος της δοκού, αντίστοιχα, οπότε πρέπει να γίνει η επιλογή `Aplly y-y`.

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_beam_section.png
---
width: 650px
align: center
name: pi_frame2d_beam_section
---
```

#### Στύλοι

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_col_section.png
---
width: 650px
align: center
name: pi_frame2d_col_section
---
```

## Γεωμετρία φορέα

### Κόμβοι

Οι κόμβοι δίνονται μέσω των συντεταγμένων τους και όχι με γραφικό τρόπο.
Οι συνθήκες στήριξης ορίζονται παγιώνοντας τους σχετικούς βαθμούς ελευθερίας. Στην περίπτωση που εξετάζεται υπάρχουν πακτώσεις στη βάση του φορέα, οπότε παγιώνονται όλοι οι βαθμοί ελευθερίας στους κόμβους 1(Α) και 4(D).

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_nodes.png
---
width: 550px
align: center
name: pi_frame2d_nodes
---
```

Στο σχήμα που ακολουθεί φαίνεται η τοποθέτηση και η ονομασία των κόμβων του φορέα.

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_nodes2.png
---
width: 450px
align: center
name: pi_frame2d_nodes2
---
```

### Δομικά στοιχεία

Τα δομικά στοιχεία εισάγονται ενώνοντας τους κόμβους, επιλέγοντας το κατάλληλο υλικό και τη διατομή τους. Στην περίπτωση που υπάρχουν εσωτερικές αρθρώσεις στην αρχή ή στο τέλος ενός στοιχείου επιλέγεται το σχετικό άκρο στα `Hinges` (δεν υπάρχει άρθρωση στο φορέα που εξετάζεται).

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_elements.png
---
width: 600px
align: center
name: pi_frame2d_elements
---
```

## Φορτία

### Εισαγωγή επικόμβιων φορτίων

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_nodal_loads.png
---
width: 550px
align: center
name: pi_frame2d_nodal_loads
---
```

### Εισαγωγή ομοιόμορφων φορτίων δοκού

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_element_loads.png
---
width: 550px
align: center
name: pi_frame2d_element_loads
---
```

## Αποτελέσματα ανάλυσης

### Διάγραμμα ροπών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_results_moments.png
---
width: 600px
align: center
name: pi_frame2d_results_moments
---
```

### Διάγραμμα τεμνουσών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_results_shears.png
---
width: 600px
align: center
name: pi_frame2d_results_shears
---
```

### Διάγραμμα αξονικών

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_results_axials.png
---
width: 600px
align: center
name: pi_frame2d_results_axials
---
```

Αν τα διαγράμματα εμφανίζονται με περίεργο scaling (π.χ. είναι τεράστια), τότε με τα `+` και `-` κουμπιά, πάνω αριστερά, διορθώνεται το μέγεθος ώστε να εμφανίζονται λογικά.

### Παρατηρήσεις

Τα αποτελέσματα είναι όμοια με αυτά που προκύπτουν με την *κλασική στατική*, άρα και με την επίλυση με το SAP2000 όταν αγνοείται το έργο των τεμνουσών δυνάμεων.
